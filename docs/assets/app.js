/* Interactive supplementary tables. No dependencies.
   Data comes from docs/data/*.json, written by scripts/build_supp_tables.py. */
(function () {
  "use strict";

  const CONFIG = {
    paperUrl: "https://www.medrxiv.org/content/10.64898/2025.12.23.25342890",
    repoUrl: "https://github.com/VilteBaltra/genetic-architecture-of-brain-age-gap",
    pageSize: 100,
    sigThreshold: 0.05,
  };

  const $ = (sel, el = document) => el.querySelector(sel);
  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  const NA = new Set(["NA", "N/A", "NaN", "–", "-"]);
  const isBlank = (v) => v === null || v === undefined || (typeof v === "string" && NA.has(v));

  let manifest = null;
  const cache = new Map();      // id -> table data
  const state = new Map();      // id -> view state
  let current = null;
  let searchIndex = null;

  $("#link-paper").href = CONFIG.paperUrl;
  $("#link-repo").href = CONFIG.repoUrl;

  /* ---------- formatting ---------- */
  function fmtNum(v, col) {
    if (typeof v !== "number") return esc(v);
    if (col.kind === "id" || Number.isInteger(v)) return String(v);
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-3 || a >= 1e6)) {
      const [m, e] = v.toExponential(2).split("e");
      return `${m}×10<sup>${Number(e)}</sup>`;
    }
    return String(Number(v.toPrecision(4)));
  }

  function plainNum(v, col) {
    if (typeof v !== "number") return String(v);
    if (col.kind === "id" || Number.isInteger(v)) return String(v);
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-3 || a >= 1e6)) return v.toExponential(2);
    return String(Number(v.toPrecision(4)));
  }

  function linkFor(v, col) {
    const s = String(v);
    switch (col.kind) {
      case "snp": return /^rs\d+$/.test(s) ? `https://www.ncbi.nlm.nih.gov/snp/${s}` : null;
      case "ensg": return /^ENSG\d+$/.test(s) ? `https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${s}` : null;
      case "gene": {
        const g = s.replace(/[()]/g, "").trim();
        return /^[A-Za-z0-9.\-]+$/.test(g) ? `https://www.genecards.org/cgi-bin/carddisp.pl?gene=${encodeURIComponent(g)}` : null;
      }
      case "pmid": return /^\d+$/.test(s) ? `https://pubmed.ncbi.nlm.nih.gov/${s}/` : null;
      case "url": return /^https?:\/\//.test(s) ? s : /^www\./.test(s) ? `https://${s}` : null;
    }
    return null;
  }

  function cellHTML(v, col) {
    if (v === null || v === undefined) return "";
    if (typeof v === "string" && NA.has(v)) return `<span class="na">${esc(v)}</span>`;
    const href = linkFor(v, col);
    if (href) return `<a href="${esc(href)}" target="_blank" rel="noopener" data-stop>${esc(v)}</a>`;
    if (typeof v === "number") return fmtNum(v, col);
    let s = esc(v);
    // Link "PMID: 12345" and bare URLs inside free text.
    // One pass, so a link made for a PMID is never re-linked as a URL.
    s = s.replace(/(PMID:?\s*(\d{6,9}))|(https?:\/\/[^\s<&]+[^\s<&.,;)])/g, (m, pm, id, url) => pm
      ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${id}/" target="_blank" rel="noopener" data-stop>${pm}</a>`
      : `<a href="${url}" target="_blank" rel="noopener" data-stop>${url}</a>`);
    return s;
  }

  /* ---------- filtering ---------- */
  const NUMRE = "[-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?";
  function makeFilter(text, col) {
    const t = text.trim();
    if (!t) return null;
    if (col.type === "num") {
      let m = t.match(new RegExp(`^(<=|>=|<|>|=)?\\s*(${NUMRE})$`));
      if (m) {
        const op = m[1] || "=", x = parseFloat(m[2]);
        return (v) => typeof v === "number" && (op === "<" ? v < x : op === ">" ? v > x : op === "<=" ? v <= x : op === ">=" ? v >= x : v === x);
      }
      m = t.match(new RegExp(`^(${NUMRE})\\s*(?:\\.\\.|to|–|—)\\s*(${NUMRE})$`));
      if (m) {
        const lo = parseFloat(m[1]), hi = parseFloat(m[2]);
        return (v) => typeof v === "number" && v >= lo && v <= hi;
      }
    }
    const q = t.toLowerCase();
    if (q.startsWith("=")) {
      const exact = q.slice(1).trim();
      return (v) => v !== null && String(v).toLowerCase() === exact;
    }
    return (v) => v !== null && String(v).toLowerCase().includes(q);
  }

  /* ---------- data loading ---------- */
  async function getJSON(path) {
    const r = await fetch(path);
    if (!r.ok) throw new Error(`${path}: HTTP ${r.status}`);
    return r.json();
  }

  async function loadTable(id) {
    if (cache.has(id)) return cache.get(id);
    const t = await getJSON(`data/${id}.json`);
    t.rows.forEach((r, i) => { r._i = i; r._text = r.map((v) => (v === null ? "" : String(v))).join("\u0001").toLowerCase(); });
    t.pCols = t.columns.map((c, i) => (c.kind === "p" ? i : -1)).filter((i) => i >= 0);
    // Main p-value column: the plain "pval"/"p"/"P.uncorrected"-style one, else the first p column.
    const MAIN_P = ["pval", "p", "p-value", "pvalue", "p.uncorrected", "top p"];
    t.mainP = t.pCols.find((i) => MAIN_P.includes(t.columns[i].name.toLowerCase()));
    if (t.mainP === undefined) t.mainP = t.pCols.length ? t.pCols[0] : -1;
    t.stickyCol = t.columns.findIndex((c) => c.kind !== "section");
    if (t.heatmap) {
      let lo = Infinity, hi = -Infinity;
      t.heatCols = new Set();
      t.columns.forEach((c, i) => { if (c.type === "num" && c.kind !== "id") t.heatCols.add(i); });
      for (const r of t.rows) for (const i of t.heatCols) if (typeof r[i] === "number") { lo = Math.min(lo, r[i]); hi = Math.max(hi, r[i]); }
      t.heatRange = [lo, hi];
    }
    cache.set(id, t);
    return t;
  }

  function getState(id) {
    if (!state.has(id)) state.set(id, { q: "", filters: {}, sort: null, page: 0, hidden: new Set(), sigOnly: false, sigCol: undefined });
    return state.get(id);
  }

  /* ---------- sidebar ---------- */
  function buildTOC() {
    const toc = $("#toc"), sel = $("#table-select");
    let html = "", opts = "";
    for (const g of manifest.groups) {
      const items = manifest.tables.filter((t) => t.group === g.key);
      if (!items.length) continue;
      html += `<section><h2>${esc(g.label)}</h2><ul>`;
      opts += `<optgroup label="${esc(g.label)}">`;
      for (const t of items) {
        html += `<li><a href="#${t.id}" data-id="${t.id}" title="${esc(t.caption)}"><span class="tid">${t.id}</span><span class="cap">${esc(t.caption)}</span></a></li>`;
        opts += `<option value="${t.id}">${t.id} · ${esc(t.caption)}</option>`;
      }
      html += `</ul></section>`;
      opts += `</optgroup>`;
    }
    toc.innerHTML = html;
    sel.innerHTML = opts;
    sel.addEventListener("change", () => { location.hash = sel.value; });
  }

  function markCurrent(id) {
    document.querySelectorAll("#toc a").forEach((a) => a.setAttribute("aria-current", a.dataset.id === id ? "true" : "false"));
    $("#table-select").value = id;
  }

  /* ---------- table view ---------- */
  async function show(id, opts = {}) {
    const meta = manifest.tables.find((t) => t.id === id) || manifest.tables[0];
    id = meta.id;
    markCurrent(id);
    const main = $("#main");
    if (!cache.has(id)) main.innerHTML = `<div class="loading">Loading ${esc(id)} (${meta.rows.toLocaleString()} rows)…</div>`;
    let t;
    try { t = await loadTable(id); } catch (e) {
      main.innerHTML = `<div class="empty">Could not load ${esc(id)}. ${esc(e.message)}. If you opened index.html straight from disk, serve the folder instead (for example <code>python -m http.server</code>).</div>`;
      return;
    }
    current = id;
    const st = getState(id);
    if (opts.q !== undefined) { st.q = opts.q; st.page = 0; }
    renderShell(t, st);
    renderBody(t, st);
    if (opts.scroll !== false) main.scrollIntoView({ block: "start" });
    document.title = `${id} · Brain Age Gap Supplementary Tables`;
  }

  // The section column is shown as divider rows instead, except while sorted.
  function visibleCols(t, st) {
    const sec = t.columns.findIndex((c) => c.kind === "section");
    return t.columns.map((c, i) => i).filter((i) => !st.hidden.has(i) && !(i === sec && !st.sort));
  }

  function isActive(st) {
    return !!(st.q.trim() || Object.keys(st.filters).length || st.sort || st.sigOnly || st.hidden.size);
  }

  function renderShell(t, st) {
    const main = $("#main");
    const hasP = t.pCols.length > 0;
    if (st.sigCol === undefined) st.sigCol = t.mainP;
    const pLabel = (i) => (t.columns[i].group ? `${t.columns[i].group} · ` : "") + t.columns[i].name;
    const wide = t.columns.length > 10;
    main.innerHTML = `
      <div class="t-head">
        <span class="t-id">Table ${t.id}</span>
        <h2>${esc(t.caption)}</h2>
        <span class="t-meta">${t.rows.length.toLocaleString()} rows · ${t.columns.length} columns · click a row to see all its values</span>
      </div>
      <div class="toolbar">
        <input id="tq" class="grow" type="search" placeholder="Search this table" value="${esc(st.q)}" autocomplete="off" spellcheck="false" aria-label="Search this table">
        <button id="filt" class="btn ${st.showFilters ? "on" : ""}" type="button" aria-pressed="${!!st.showFilters}">Filter columns</button>
        ${hasP ? `<span class="sigctl">
          <button id="sig" class="btn ${st.sigOnly ? "on" : ""}" type="button" aria-pressed="${st.sigOnly}" title="Keep rows where the selected p-value column is below ${CONFIG.sigThreshold}">p &lt; ${CONFIG.sigThreshold}</button>
          ${t.pCols.length > 1
            ? `<label class="sigcol"><span>in</span><select id="sigcol" aria-label="p-value column used for the p &lt; ${CONFIG.sigThreshold} filter">${t.pCols.map((i) => `<option value="${i}" ${i === st.sigCol ? "selected" : ""}>${esc(pLabel(i))}</option>`).join("")}</select></label>`
            : `<span class="sigcol"><span>in ${esc(pLabel(st.sigCol))}</span></span>`}
        </span>` : ""}
        ${wide ? `<div class="colpicker"><button id="colbtn" class="btn" type="button" aria-expanded="false">Columns</button><div id="colmenu" class="colmenu" hidden></div></div>` : ""}
        <span class="quiet">
          <button id="reset" class="btn ghost" type="button" hidden>Clear all</button>
          <button id="csv" class="btn ghost" type="button">Download CSV</button>
          <button id="copylink" class="btn ghost" type="button">Copy link</button>
        </span>
      </div>
      <p class="hint" id="hint" ${st.showFilters ? "" : "hidden"}>Type in the boxes under each heading. Numbers accept <code>&lt;0.05</code>, <code>&gt;=10</code> or a range like <code>1e-8..1e-5</code>. Start with <code>=</code> for an exact text match.</p>
      <div class="scroller" id="scroller"><table class="data" id="grid"></table></div>
      <div class="pager"><span class="status" id="status"></span><span class="btns" id="pagebtns"></span></div>
      <div class="legend" id="legend"></div>
      <div class="notes" id="notes"></div>`;

    let timer;
    $("#tq").addEventListener("input", (e) => { clearTimeout(timer); timer = setTimeout(() => { st.q = e.target.value; st.page = 0; renderBody(t, st); }, 150); });
    $("#filt").addEventListener("click", (e) => {
      st.showFilters = !st.showFilters;
      e.currentTarget.classList.toggle("on", st.showFilters); e.currentTarget.setAttribute("aria-pressed", st.showFilters);
      $("#hint").hidden = !st.showFilters;
      renderBody(t, st);
      if (st.showFilters) { const f = $("#grid [data-filter]"); if (f) f.focus(); }
    });
    if (hasP && $("#sigcol") && $("#sigcol").tagName === "SELECT") $("#sigcol").addEventListener("change", (e) => {
      st.sigCol = +e.target.value; st.page = 0;
      const lg = $("#sig-legend"); if (lg) lg.textContent = `${pLabel(st.sigCol)} below ${CONFIG.sigThreshold}`;
      renderBody(t, st);
    });
    if (hasP) $("#sig").addEventListener("click", (e) => {
      st.sigOnly = !st.sigOnly; st.page = 0;
      e.currentTarget.classList.toggle("on", st.sigOnly); e.currentTarget.setAttribute("aria-pressed", st.sigOnly);
      renderBody(t, st);
    });
    $("#reset").addEventListener("click", () => { state.delete(t.id); show(t.id, { scroll: false }); });
    $("#csv").addEventListener("click", () => downloadCSV(t, st));
    $("#copylink").addEventListener("click", () => copyLink(t.id));
    if (wide) setupColumnMenu(t, st);

    // Legend and notes
    const leg = [];
    if (hasP) leg.push(`<span><span class="sig">●</span> <span id="sig-legend">${esc(pLabel(st.sigCol))} below ${CONFIG.sigThreshold}</span></span>`);
    if (t.heatmap) leg.push(`<span><span class="sw" style="background:linear-gradient(90deg, transparent, color-mix(in srgb, var(--heat) 60%, transparent))"></span>${fmtNum(t.heatRange[0], {})} → ${fmtNum(t.heatRange[1], {})}</span>`);
    $("#legend").innerHTML = leg.join("");
    let notes = t.notes.map((n) => `<p>${cellHTML(n, {})}</p>`).join("");
    if (t.glossary.length) notes += `<h3>Column definitions</h3><dl class="glossary">${t.glossary.map(([k, v]) => `<dt>${esc(k)}</dt><dd>${esc(v)}</dd>`).join("")}</dl>`;
    $("#notes").innerHTML = notes;
  }

  function setupColumnMenu(t, st) {
    const btn = $("#colbtn"), menu = $("#colmenu");
    const groups = new Map();
    t.columns.forEach((c, i) => { const g = c.group || "Columns"; if (!groups.has(g)) groups.set(g, []); groups.get(g).push(i); });
    const draw = () => {
      menu.innerHTML = `<div class="row"><button class="btn" type="button" data-all="1">Show all</button></div>` +
        [...groups].map(([g, idx], gi) => `<fieldset><legend><label><input type="checkbox" data-group="${gi}" ${idx.every((i) => !st.hidden.has(i)) ? "checked" : ""}> ${esc(g)}</label></legend>` +
          idx.map((i) => `<label><input type="checkbox" data-col="${i}" ${st.hidden.has(i) ? "" : "checked"} ${i === t.stickyCol ? "disabled" : ""}> ${esc(t.columns[i].name)}</label>`).join("") + `</fieldset>`).join("");
    };
    draw();
    const gList = [...groups.values()];
    btn.addEventListener("click", () => { menu.hidden = !menu.hidden; btn.setAttribute("aria-expanded", String(!menu.hidden)); });
    document.addEventListener("click", (e) => { if (!menu.hidden && !e.target.closest(".colpicker")) { menu.hidden = true; btn.setAttribute("aria-expanded", "false"); } });
    menu.addEventListener("change", (e) => {
      const el = e.target;
      if (el.dataset.col) { const i = +el.dataset.col; el.checked ? st.hidden.delete(i) : st.hidden.add(i); }
      if (el.dataset.group) { for (const i of gList[+el.dataset.group]) if (i !== t.stickyCol) el.checked ? st.hidden.delete(i) : st.hidden.add(i); }
      draw(); renderBody(t, st);
    });
    menu.addEventListener("click", (e) => { if (e.target.dataset.all) { st.hidden.clear(); draw(); renderBody(t, st); } });
  }

  function filteredRows(t, st) {
    const q = st.q.trim().toLowerCase();
    const fs = Object.entries(st.filters).map(([i, txt]) => [+i, makeFilter(txt, t.columns[+i])]).filter(([, f]) => f);
    let rows = t.rows.filter((r) => {
      if (q && !r._text.includes(q)) return false;
      for (const [i, f] of fs) if (!f(r[i])) return false;
      if (st.sigOnly && !(typeof r[st.sigCol] === "number" && r[st.sigCol] < CONFIG.sigThreshold)) return false;
      return true;
    });
    if (st.sort) {
      const { col, dir } = st.sort;
      const num = t.columns[col].type === "num";
      rows = rows.slice().sort((a, b) => {
        const x = a[col], y = b[col];
        const bx = isBlank(x), by = isBlank(y);
        if (bx || by) return bx === by ? a._i - b._i : bx ? 1 : -1;   // blanks last
        let c;
        if (num && typeof x === "number" && typeof y === "number") c = x - y;
        else if (typeof x === "number" && typeof y !== "number") c = -1;
        else if (typeof y === "number" && typeof x !== "number") c = 1;
        else c = String(x).localeCompare(String(y), undefined, { numeric: true, sensitivity: "base" });
        return c === 0 ? a._i - b._i : c * dir;
      });
    }
    return rows;
  }

  function cellClass(v, col, t, i, st) {
    const cls = [];
    if (col.kind === "section") cls.push("section");
    else if (col.type === "num") cls.push("num");
    else if (["snp", "ensg", "id"].includes(col.kind)) cls.push("mono");
    else if (typeof v === "string" && v.length > 40) cls.push("trunc");
    if (i === st.sigCol && typeof v === "number" && v < CONFIG.sigThreshold) cls.push("p-sig");
    if (i === t.stickyCol) cls.push("sticky");
    return cls.join(" ");
  }

  function heatStyle(t, i, v) {
    if (!t.heatmap || !t.heatCols.has(i) || typeof v !== "number") return "";
    const [lo, hi] = t.heatRange;
    const f = hi > lo ? (v - lo) / (hi - lo) : 0;
    return ` style="background:color-mix(in srgb, var(--heat) ${Math.round(f * 60)}%, var(--surface))"`;
  }

  function renderHead(t, st, cols) {
    const hasGroups = cols.some((i) => t.columns[i].group);
    let groupRow = "";
    if (hasGroups) {
      let k = 0;
      while (k < cols.length) {
        const g = t.columns[cols[k]].group;
        let span = 1;
        while (k + span < cols.length && t.columns[cols[k + span]].group === g) span++;
        const sticky = cols[k] === t.stickyCol && span === 1 ? " sticky" : "";
        groupRow += g ? `<th colspan="${span}" class="g${sticky}"><span>${esc(g)}</span></th>` : `<th colspan="${span}" class="${sticky.trim()}"></th>`;
        k += span;
      }
      groupRow = `<tr class="groups">${groupRow}</tr>`;
    }
    const names = cols.map((i) => {
      const c = t.columns[i];
      const s = st.sort && st.sort.col === i ? `<span class="arrow" aria-hidden="true">${st.sort.dir > 0 ? "↑" : "↓"}</span>` : "";
      const aria = st.sort && st.sort.col === i ? (st.sort.dir > 0 ? "ascending" : "descending") : "none";
      return `<th class="sortable ${c.type === "num" ? "num" : ""} ${i === t.stickyCol ? "sticky" : ""}" data-col="${i}" aria-sort="${aria}" tabindex="0" title="Sort by ${esc(c.name)}">${esc(c.name)}${s}</th>`;
    }).join("");
    const filters = cols.map((i) => {
      const c = t.columns[i];
      const ph = c.type === "num" ? (c.kind === "p" ? "<0.05" : "e.g. >0") : "filter";
      return `<th class="${i === t.stickyCol ? "sticky" : ""}"><input type="text" data-filter="${i}" id="f-${t.id}-${i}" value="${esc(st.filters[i] || "")}" placeholder="${ph}" aria-label="Filter ${esc(c.name)}" autocomplete="off" spellcheck="false"></th>`;
    }).join("");
    const showF = st.showFilters || Object.keys(st.filters).length > 0;
    return `<thead>${groupRow}<tr class="names">${names}</tr>${showF ? `<tr class="filters">${filters}</tr>` : ""}</thead>`;
  }

  function renderBody(t, st) {
    const cols = visibleCols(t, st);
    const rows = filteredRows(t, st);
    const pages = Math.max(1, Math.ceil(rows.length / CONFIG.pageSize));
    st.page = Math.min(st.page, pages - 1);
    const start = st.page * CONFIG.pageSize;
    const slice = rows.slice(start, start + CONFIG.pageSize);
    const secCol = t.columns.findIndex((c) => c.kind === "section");

    let body = "";
    let prevSec;
    for (const r of slice) {
      if (secCol >= 0 && !st.sort && r[secCol] !== prevSec) {
        body += `<tr class="sec-row"><td colspan="${cols.length}"><span>${esc(r[secCol] || "")}</span></td></tr>`;
      }
      prevSec = secCol >= 0 ? r[secCol] : null;
      body += `<tr data-i="${r._i}">` + cols.map((i) => {
        const v = r[i], c = t.columns[i];
        const inner = cellHTML(v, c);
        const title = typeof v === "number" && !Number.isInteger(v) ? ` title="${v}"`
          : typeof v === "string" && v.length > 40 ? ` title="${esc(v)}"` : "";
        return `<td class="${cellClass(v, c, t, i, st)}"${heatStyle(t, i, v)}${title}>${inner}</td>`;
      }).join("") + `</tr>`;
    }
    if (!slice.length) body = `<tr><td colspan="${cols.length}" class="empty">No rows match. Clear the search or filters to see all ${t.rows.length.toLocaleString()} rows.</td></tr>`;

    // Keep focus in a filter input across re-renders.
    const active = document.activeElement && document.activeElement.dataset ? document.activeElement.dataset.filter : undefined;
    const caret = active !== undefined ? document.activeElement.selectionStart : null;

    const grid = $("#grid");
    grid.innerHTML = renderHead(t, st, cols) + `<tbody>${body}</tbody>`;

    if (active !== undefined) {
      const el = grid.querySelector(`[data-filter="${active}"]`);
      if (el) { el.focus(); try { el.setSelectionRange(caret, caret); } catch (e) { /* ignore */ } }
    }

    // Header interactions
    grid.querySelectorAll("th.sortable").forEach((th) => {
      const go = () => {
        const i = +th.dataset.col;
        if (!st.sort || st.sort.col !== i) st.sort = { col: i, dir: 1 };
        else if (st.sort.dir === 1) st.sort.dir = -1;
        else st.sort = null;
        st.page = 0; renderBody(t, st);
      };
      th.addEventListener("click", go);
      th.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); go(); } });
    });
    let timer;
    grid.querySelectorAll("[data-filter]").forEach((inp) => inp.addEventListener("input", () => {
      clearTimeout(timer);
      timer = setTimeout(() => { const i = inp.dataset.filter; if (inp.value.trim()) st.filters[i] = inp.value; else delete st.filters[i]; st.page = 0; renderBody(t, st); }, 200);
    }));
    grid.querySelector("tbody").addEventListener("click", (e) => {
      if (e.target.closest("[data-stop]")) return;
      const tr = e.target.closest("tr[data-i]");
      if (tr) openDrawer(t, t.rows[+tr.dataset.i]);
    });

    const reset = $("#reset"); if (reset) reset.hidden = !isActive(st);

    // Status and pager
    const total = t.rows.length;
    const filtered = rows.length !== total;
    $("#status").textContent = rows.length
      ? `${(start + 1).toLocaleString()}–${(start + slice.length).toLocaleString()} of ${rows.length.toLocaleString()} rows${filtered ? ` (filtered from ${total.toLocaleString()})` : ""}`
      : `0 of ${total.toLocaleString()} rows`;
    const pb = $("#pagebtns");
    pb.innerHTML = pages > 1
      ? `<button class="btn" type="button" data-p="first" ${st.page === 0 ? "disabled" : ""}>First</button>
         <button class="btn" type="button" data-p="prev" ${st.page === 0 ? "disabled" : ""}>Previous</button>
         <span class="status">Page ${st.page + 1} of ${pages}</span>
         <button class="btn" type="button" data-p="next" ${st.page >= pages - 1 ? "disabled" : ""}>Next</button>
         <button class="btn" type="button" data-p="last" ${st.page >= pages - 1 ? "disabled" : ""}>Last</button>`
      : "";
    pb.onclick = (e) => {
      const p = e.target.dataset.p; if (!p) return;
      st.page = p === "first" ? 0 : p === "prev" ? st.page - 1 : p === "next" ? st.page + 1 : pages - 1;
      renderBody(t, st); $("#scroller").scrollTop = 0;
    };
  }

  /* ---------- row drawer ---------- */
  function openDrawer(t, r) {
    const d = $("#drawer");
    const key = r[t.stickyCol];
    $("#drawer-title").textContent = `${t.id} · ${key === null ? "Row details" : key}`;
    $("#drawer-body").innerHTML = t.columns.map((c, i) => {
      const v = r[i];
      if (v === null) return "";
      const val = typeof v === "number" ? (Number.isInteger(v) ? v : `${fmtNum(v, c)} <span class="muted">(${v})</span>`) : cellHTML(v, c);
      return `<dt>${c.group ? `<span class="grp">${esc(c.group)} · </span>` : ""}${esc(c.name)}</dt><dd class="${typeof v === "number" ? "num" : ""}">${val}</dd>`;
    }).join("");
    d.hidden = false;
    $("#drawer-close").focus();
  }
  $("#drawer-close").addEventListener("click", () => { $("#drawer").hidden = true; });
  document.addEventListener("keydown", (e) => { if (e.key === "Escape") $("#drawer").hidden = true; });

  /* ---------- export ---------- */
  function downloadCSV(t, st) {
    const cols = visibleCols(t, st);
    const rows = filteredRows(t, st);
    const q = (v) => { if (v === null || v === undefined) return ""; const s = String(v); return /[",\n\r]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; };
    const head = cols.map((i) => q(t.columns[i].group ? `${t.columns[i].group} | ${t.columns[i].name}` : t.columns[i].name)).join(",");
    const lines = rows.map((r) => cols.map((i) => q(r[i])).join(","));
    const blob = new Blob(["﻿" + [head, ...lines].join("\r\n")], { type: "text/csv;charset=utf-8" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `Table_${t.id}.csv`;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
    toast(`Saved ${rows.length.toLocaleString()} rows as Table_${t.id}.csv`);
  }

  function copyLink(id) {
    const url = location.href.split("#")[0] + "#" + id;
    const done = () => toast("Link copied");
    try { navigator.clipboard.writeText(url).then(done, () => prompt(url)); } catch (e) { prompt(url); }
    function prompt(u) { toast(u); }
  }

  let toastTimer;
  function toast(msg) {
    let el = $(".toast");
    if (!el) { el = document.createElement("div"); el.className = "toast"; el.setAttribute("role", "status"); document.body.appendChild(el); }
    el.textContent = msg; el.hidden = false;
    clearTimeout(toastTimer); toastTimer = setTimeout(() => { el.hidden = true; }, 2600);
  }

  /* ---------- cross-table search ---------- */
  function setupGlobalSearch() {
    const input = $("#global-q"), out = $("#global-results");
    let timer;
    input.addEventListener("input", () => {
      clearTimeout(timer);
      timer = setTimeout(async () => {
        const q = input.value.trim().toLowerCase();
        if (q.length < 2) { out.hidden = true; return; }
        if (!searchIndex) {
          out.hidden = false; out.innerHTML = `<span class="muted">Loading search index…</span>`;
          try { searchIndex = await getJSON("data/search.json"); } catch (e) { out.innerHTML = `<span class="muted">Search index unavailable.</span>`; return; }
        }
        const hits = [];
        for (const key in searchIndex) {
          const k = key.toLowerCase();
          if (k.includes(q)) hits.push([key, k === q ? 0 : k.startsWith(q) ? 1 : 2]);
        }
        hits.sort((a, b) => a[1] - b[1] || a[0].length - b[0].length || a[0].localeCompare(b[0]));
        const shown = hits.slice(0, 12);
        out.hidden = false;
        if (!shown.length) { out.innerHTML = `<span class="muted">No table contains “${esc(input.value.trim())}”.</span>`; return; }
        out.innerHTML = shown.map(([key]) => {
          const tabs = Object.entries(searchIndex[key]).sort((a, b) => parseInt(a[0].slice(1)) - parseInt(b[0].slice(1)) || a[0].localeCompare(b[0]));
          return `<div class="hit"><span class="val">${esc(key)}</span><span class="chips">${tabs.map(([id, n]) => `<button class="chip" type="button" data-id="${id}" data-q="${esc(key)}" title="${n} row${n > 1 ? "s" : ""} in ${id}">${id}${n > 1 ? ` ×${n}` : ""}</button>`).join("")}</span></div>`;
        }).join("") + (hits.length > shown.length ? `<span class="muted">${hits.length - shown.length} more values match. Type more to narrow.</span>` : "");
      }, 180);
    });
    out.addEventListener("click", (e) => {
      const b = e.target.closest(".chip"); if (!b) return;
      const id = b.dataset.id, q = b.dataset.q;
      getState(id).filters = {};
      if (location.hash.slice(1) === id) show(id, { q }); else { pendingQuery = q; location.hash = id; }
    });
  }
  let pendingQuery;

  /* ---------- routing ---------- */
  function route() {
    const id = decodeURIComponent((location.hash || "").slice(1)).toLowerCase();
    const hit = manifest.tables.find((t) => t.id.toLowerCase() === id);
    const target = hit ? hit.id : manifest.tables[0].id;
    const q = pendingQuery; pendingQuery = undefined;
    show(target, { q, scroll: !!id });
  }

  (async function init() {
    try { manifest = await getJSON("data/manifest.json"); } catch (e) {
      $("#main").innerHTML = `<div class="empty">Could not load the table list (${esc(e.message)}). If you opened index.html straight from disk, serve the folder instead, for example with <code>python -m http.server</code> inside <code>docs/</code>.</div>`;
      return;
    }
    buildTOC();
    setupGlobalSearch();
    window.addEventListener("hashchange", route);
    route();
  })();
})();
