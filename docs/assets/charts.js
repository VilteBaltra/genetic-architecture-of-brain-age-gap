/* Optional chart views for selected supplementary tables.
   To turn charts off, delete the <script src="assets/charts.js"> line in index.html.
   Each entry in window.SuppCharts adds a Table / Chart switch to that table. */
(function () {
  "use strict";

  /* Colours: three method hues validated for colour-blind separation and contrast
     against the light and dark surfaces. Shapes repeat the encoding. */
  const css = `
  :root { --m1: #2a6fa3; --m2: #c4651f; --m3: #2e8a68; }
  @media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --m1: #4a9ad6; --m2: #cf7a2c; --m3: #339c76; } }
  :root[data-theme="dark"] { --m1: #4a9ad6; --m2: #cf7a2c; --m3: #339c76; }
  .chartwrap { position: relative; background: var(--surface); border: 1px solid var(--line); border-radius: var(--radius); padding: 1rem 1rem .5rem; }
  .fp-head { display: flex; flex-wrap: wrap; gap: .5rem 1.25rem; align-items: center; margin-bottom: .5rem; font-size: var(--fs-sm); }
  .fp-title { font-weight: 600; }
  .fp-controls { display: flex; flex-wrap: wrap; gap: .5rem; }
  .fp-controls label { display: inline-flex; align-items: center; gap: .35rem; color: var(--muted); font-size: var(--fs-sm); }
  .fp-controls select { width: auto; padding: .35rem .5rem; }
  .fp-legend { display: flex; flex-wrap: wrap; gap: .3rem 1rem; color: var(--muted); font-size: var(--fs-xs); margin-bottom: .5rem; }
  .fp-legend span { display: inline-flex; align-items: center; gap: .35rem; }
  .fp-legend svg { overflow: visible; }
  .fp-caption { color: var(--muted); font-size: var(--fs-xs); margin: .25rem 0 0; max-width: 80ch; }
  .fp-svg { display: block; width: 100%; height: auto; overflow: visible; }
  .fp-svg text { font-family: var(--font-body); }
  .fp-svg .lab { fill: var(--ink); font-size: 12.5px; }
  .fp-svg .tick { fill: var(--muted); font-size: 11px; font-family: var(--font-data); }
  .fp-svg .axt { fill: var(--muted); font-size: 11.5px; }
  .fp-svg .grid { stroke: var(--line); stroke-width: 1; }
  .fp-svg .zero { stroke: var(--muted); stroke-width: 1; stroke-dasharray: 3 3; }
  .fp-svg .band { fill: var(--surface-2); }
  .fp-svg .ci { stroke-width: 2; stroke-linecap: round; }
  .fp-svg .mk { stroke-width: 2; }
  .fp-svg .hollow { fill: var(--surface); }
  .fp-svg .s1 { stroke: var(--m1); } .fp-svg .f1 { fill: var(--m1); }
  .fp-svg .s2 { stroke: var(--m2); } .fp-svg .f2 { fill: var(--m2); }
  .fp-svg .s3 { stroke: var(--m3); } .fp-svg .f3 { fill: var(--m3); }
  .fp-svg .hit { fill: transparent; cursor: pointer; }
  .fp-svg g[data-g] { pointer-events: none; }
  .fp-svg .hit:hover + g .mk, .fp-svg g.on .mk { stroke-width: 3; }
  .fp-tip { position: absolute; z-index: 5; pointer-events: none; background: var(--ink); color: var(--bg); border-radius: 4px; padding: .4rem .55rem; font-size: var(--fs-xs); line-height: 1.4; max-width: 18rem; box-shadow: var(--shadow); }
  .fp-tip b { font-weight: 600; }
  .fp-empty { padding: 2rem; text-align: center; color: var(--muted); }
  `;
  const style = document.createElement("style");
  style.textContent = css;
  document.head.appendChild(style);

  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  const num = (v) => (typeof v === "number" && isFinite(v) ? v : null);
  const fmt = (v) => {
    if (v === null || v === undefined) return "–";
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-3 || a >= 1e5)) return v.toExponential(2).replace("e", "×10^");
    return String(Number(v.toPrecision(3)));
  };
  const col = (t, ...names) => t.columns.findIndex((c) => names.includes(c.name.toLowerCase()));

  function niceTicks(a, b, n) {
    const span = b - a || 1;
    const raw = span / n;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const e = raw / mag;
    const step = (e >= 7.5 ? 10 : e >= 3.5 ? 5 : e >= 1.5 ? 2 : 1) * mag;
    const out = [];
    for (let v = Math.ceil(a / step) * step; v <= b + step * 1e-9; v += step) out.push(Number(v.toFixed(12)));
    return out;
  }

  function marker(shape, x, y, r, cls) {
    if (shape === "diamond") {
      const d = r * 1.3;
      return `<path class="${cls}" d="M${x},${y - d}L${x + d},${y}L${x},${y + d}L${x - d},${y}Z"/>`;
    }
    if (shape === "square") return `<rect class="${cls}" x="${x - r * 0.9}" y="${y - r * 0.9}" width="${r * 1.8}" height="${r * 1.8}" rx="1"/>`;
    return `<circle class="${cls}" cx="${x}" cy="${y}" r="${r}"/>`;
  }

  function legendMark(shape, n, hollow) {
    const cls = `mk s${n} ${hollow ? "hollow" : "f" + n}`;
    return `<svg class="fp-svg" width="14" height="14" viewBox="-7 -7 14 14" aria-hidden="true" style="width:14px;display:inline">${marker(shape, 0, 0, 4.5, cls)}</svg>`;
  }

  /* Generic forest plot.
     groups: [{ label, items: [{ est, lo, hi, s (series index 1-3), sig, row, tip }] }]
     series: [{ label, shape }]  (index + 1 = colour slot)
     core:   series indexes used to set the x range (others are clipped with an arrow) */
  function forest(el, o, helpers) {
    const groups = o.groups.filter((g) => g.items.length);
    if (!groups.length) {
      el.insertAdjacentHTML("beforeend", `<div class="fp-empty">No rows to chart. Clear the search or filters, or pick another option above.</div>`);
      return;
    }
    const width = Math.max(320, el.clientWidth - 32);
    const labelW = Math.min(250, Math.round(width * 0.34));
    const right = 16, axisH = 26;
    const plotL = labelW + 12, plotR = width - right;
    const per = 12, padB = 9;
    const bandH = (g) => Math.max(26, g.items.length * per + padB * 2);

    // x range from the core series, always including the null line.
    let lo = o.nullAt, hi = o.nullAt;
    for (const g of groups) for (const it of g.items) {
      if (o.core && !o.core.includes(it.s)) continue;
      if (it.lo !== null) lo = Math.min(lo, it.lo); else lo = Math.min(lo, it.est);
      if (it.hi !== null) hi = Math.max(hi, it.hi); else hi = Math.max(hi, it.est);
    }
    if (o.core) for (const g of groups) for (const it of g.items) { lo = Math.min(lo, it.est); hi = Math.max(hi, it.est); }
    const pad = (hi - lo) * 0.06 || 0.1;
    lo -= pad; hi += pad;
    const ticks = niceTicks(lo, hi, Math.max(3, Math.round((plotR - plotL) / 90)));
    lo = Math.min(lo, ticks[0]); hi = Math.max(hi, ticks[ticks.length - 1]);
    const X = (v) => plotL + ((v - lo) / (hi - lo)) * (plotR - plotL);
    const clampX = (v) => Math.max(plotL, Math.min(plotR, X(v)));

    let y = axisH;
    const H = axisH * 2 + 20 + groups.reduce((a, g) => a + bandH(g), 0);
    let svg = `<svg class="fp-svg" viewBox="0 0 ${width} ${H}" width="${width}" height="${H}" role="img" aria-label="${esc(o.title)}">`;
    // gridlines + ticks top and bottom
    for (const tv of ticks) {
      const x = X(tv);
      svg += `<line class="grid" x1="${x}" x2="${x}" y1="${axisH - 4}" y2="${H - axisH - 16}"/>`;
      svg += `<text class="tick" x="${x}" y="${axisH - 10}" text-anchor="middle">${fmt(tv)}</text>`;
      svg += `<text class="tick" x="${x}" y="${H - axisH - 2}" text-anchor="middle">${fmt(tv)}</text>`;
    }
    svg += `<line class="zero" x1="${X(o.nullAt)}" x2="${X(o.nullAt)}" y1="${axisH - 4}" y2="${H - axisH - 16}"/>`;

    const flat = [];
    groups.forEach((g, gi) => {
      const h = bandH(g);
      if (gi % 2 === 0) svg += `<rect class="band" x="0" y="${y}" width="${width}" height="${h}" rx="3"/>`;
      const label = g.label.length > 40 ? g.label.slice(0, 38) + "…" : g.label;
      svg += `<text class="lab" x="8" y="${y + h / 2}" dominant-baseline="middle"><title>${esc(g.label)}</title>${esc(label)}</text>`;
      const n = g.items.length;
      g.items.forEach((it, i) => {
        const cy = y + h / 2 + (i - (n - 1) / 2) * per;
        const k = flat.push(it) - 1;
        const sc = `s${it.s}`;
        let ci = "";
        if (it.lo !== null && it.hi !== null) {
          const x1 = clampX(it.lo), x2 = clampX(it.hi);
          ci += `<line class="ci ${sc}" x1="${x1}" x2="${x2}" y1="${cy}" y2="${cy}"/>`;
          if (X(it.lo) < plotL) ci += `<path class="f${it.s}" d="M${plotL},${cy}l6,-4v8z"/>`;
          if (X(it.hi) > plotR) ci += `<path class="f${it.s}" d="M${plotR},${cy}l-6,-4v8z"/>`;
        }
        const mx = clampX(it.est);
        svg += `<rect class="hit" data-k="${k}" x="${plotL - 4}" y="${cy - per / 2}" width="${plotR - plotL + 8}" height="${per}"/>`;
        svg += `<g data-g="${k}">${ci}${marker(o.series[it.s - 1].shape, mx, cy, 4.5, `mk ${sc} ${it.sig ? "f" + it.s : "hollow"}`)}</g>`;
      });
      y += h;
    });
    svg += `<text class="axt" x="${(plotL + plotR) / 2}" y="${H - 2}" text-anchor="middle">${esc(o.xLabel)}</text>`;
    svg += `</svg>`;

    const legend = (o.series.length > 1 ? o.series.map((s, i) => `<span>${legendMark(s.shape, i + 1, false)}${esc(s.label)}</span>`).join("") : "") +
      `<span>${legendMark(o.series[0].shape, 1, false)} filled: ${esc(o.sigLabel)}</span><span>${legendMark(o.series[0].shape, 1, true)} hollow: not significant</span><span>Bars: 95% CI</span>`;
    el.insertAdjacentHTML("beforeend", `<div class="fp-legend">${legend}</div><div class="fp-plot">${svg}</div>${o.caption ? `<p class="fp-caption">${o.caption}</p>` : ""}`);

    // Hover tooltip + click to open the full row.
    const tip = document.createElement("div");
    tip.className = "fp-tip"; tip.hidden = true;
    el.appendChild(tip);
    const plot = el.querySelector(".fp-plot svg");
    let on = null;
    plot.addEventListener("mousemove", (e) => {
      const h = e.target.closest(".hit");
      if (!h) { tip.hidden = true; if (on) on.classList.remove("on"); on = null; return; }
      const it = flat[+h.dataset.k];
      const g = plot.querySelector(`[data-g="${h.dataset.k}"]`);
      if (on !== g) { if (on) on.classList.remove("on"); on = g; g.classList.add("on"); }
      tip.innerHTML = it.tip;
      tip.hidden = false;
      const r = el.getBoundingClientRect();
      let left = e.clientX - r.left + 14, top = e.clientY - r.top + 12;
      if (left + tip.offsetWidth > r.width - 8) left = e.clientX - r.left - tip.offsetWidth - 14;
      tip.style.left = left + "px"; tip.style.top = top + "px";
    });
    plot.addEventListener("mouseleave", () => { tip.hidden = true; if (on) on.classList.remove("on"); on = null; });
    plot.addEventListener("click", (e) => {
      const h = e.target.closest(".hit"); if (!h) return;
      helpers.openDrawer(flat[+h.dataset.k].row);
    });
  }

  // Redraw on resize (width-dependent layout).
  let last = null, rt;
  window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(() => { if (last && document.body.contains(last[0]) && !last[0].hidden) last[5](...last.slice(0, 5)); }, 150); });

  /* ---------- LDSC genetic correlations (S13, S20) ---------- */
  function ldsc(el, t, rows, st, helpers) {
    last = [el, t, rows, st, helpers, ldsc];
    el.innerHTML = "";
    const cT = col(t, "trait"), cR = col(t, "genetic correlation", "rg"), cSE = col(t, "se"), cP = col(t, "pvalue", "p"), cF = col(t, "fdr_bh", "fdr");
    const items = rows.map((r) => {
      const est = num(r[cR]), se = num(r[cSE]);
      if (est === null) return null;
      const fdr = num(r[cF]), p = num(r[cP]);
      return {
        label: String(r[cT]),
        it: {
          est, lo: se === null ? null : est - 1.96 * se, hi: se === null ? null : est + 1.96 * se, s: 1,
          sig: fdr !== null && fdr < 0.05, row: r,
          tip: `<b>${esc(r[cT])}</b><br>rg ${fmt(est)} (95% CI ${fmt(est - 1.96 * se)} to ${fmt(est + 1.96 * se)})<br>p ${fmt(p)} · FDR ${fmt(fdr)}`,
        },
      };
    }).filter(Boolean).sort((a, b) => b.it.est - a.it.est);
    const who = t.id === "S13" ? "BAGHan" : "BAG factor";
    el.insertAdjacentHTML("beforeend", `<div class="fp-head"><span class="fp-title">Genetic correlation of ${who} with each trait</span></div>`);
    forest(el, {
      title: `Genetic correlations of ${who}`,
      groups: items.map((x) => ({ label: x.label, items: [x.it] })),
      series: [{ label: "rg", shape: "circle" }],
      nullAt: 0, xLabel: "Genetic correlation (rg)", sigLabel: "FDR < 0.05",
      caption: "Sorted by rg. 95% CI computed as rg ± 1.96 × SE. Hover a point for values; click it to see the full row.",
    }, helpers);
  }

  /* ---------- Mendelian randomisation (S22, S23, S24) ---------- */
  const METHODS = [
    { key: "inverse variance weighted", label: "Inverse variance weighted", shape: "circle" },
    { key: "mr egger", label: "MR Egger", shape: "diamond" },
    { key: "weighted median", label: "Weighted median", shape: "square" },
  ];
  function mr(el, t, rows, st, helpers) {
    last = [el, t, rows, st, helpers, mr];
    el.innerHTML = "";
    const cM = col(t, "method"), cO = col(t, "outcome"), cE = col(t, "exposure"), cB = col(t, "b"), cSE = col(t, "se"),
      cL = col(t, "lci"), cH = col(t, "hci"), cP = col(t, "pval");
    const isBag = (v) => /^bag/i.test(String(v || ""));
    const recs = [];
    for (const r of rows) {
      const fwd = isBag(r[cO]), rev = isBag(r[cE]);
      if (!fwd && !rev) continue;
      const mi = METHODS.findIndex((m) => String(r[cM]).toLowerCase() === m.key);
      if (mi < 0) continue;
      recs.push({ r, dir: fwd ? "forward" : "reverse", model: String(fwd ? r[cO] : r[cE]), trait: String(fwd ? r[cE] : r[cO]), mi });
    }
    // Options from all rows of the table, so choices don't vanish while filtering.
    const allDirs = new Set(), allModels = [];
    for (const r of t.rows) {
      const fwd = isBag(r[cO]), rev = isBag(r[cE]);
      if (!fwd && !rev) continue;
      allDirs.add(fwd ? "forward" : "reverse");
      const m = String(fwd ? r[cO] : r[cE]);
      if (!allModels.includes(m)) allModels.push(m);
    }
    allModels.sort((a, b) => (a === "BAG factor" ? -1 : b === "BAG factor" ? 1 : a.localeCompare(b)));
    st.chartOpts = st.chartOpts || {};
    const opt = st.chartOpts;
    if (!opt.dir || !allDirs.has(opt.dir)) opt.dir = allDirs.has("forward") ? "forward" : "reverse";
    if (!opt.model || !allModels.includes(opt.model)) opt.model = allModels[0];

    const dirLabel = { forward: "Trait → brain age gap (forward)", reverse: "Brain age gap → trait (reverse)" };
    el.insertAdjacentHTML("beforeend", `<div class="fp-head">
      <div class="fp-controls">
        ${allDirs.size > 1 ? `<label for="fp-dir">Direction <select id="fp-dir">${[...allDirs].map((d) => `<option value="${d}" ${d === opt.dir ? "selected" : ""}>${dirLabel[d]}</option>`).join("")}</select></label>` : `<span class="fp-title">${dirLabel[opt.dir]}</span>`}
        ${allModels.length > 1 ? `<label for="fp-model">BAG model <select id="fp-model">${allModels.map((m) => `<option ${m === opt.model ? "selected" : ""}>${esc(m)}</option>`).join("")}</select></label>` : ""}
      </div></div>`);
    const redraw = () => mr(el, t, rows, st, helpers);
    const d = el.querySelector("#fp-dir"), m = el.querySelector("#fp-model");
    if (d) d.addEventListener("change", () => { opt.dir = d.value; redraw(); });
    if (m) m.addEventListener("change", () => { opt.model = m.value; redraw(); });

    const sel = recs.filter((x) => x.dir === opt.dir && x.model === opt.model);
    const byTrait = new Map();
    for (const x of sel) {
      if (!byTrait.has(x.trait)) byTrait.set(x.trait, []);
      const b = num(x.r[cB]), se = num(x.r[cSE]);
      if (b === null) continue;
      let lo = num(x.r[cL]), hi = num(x.r[cH]);
      if (lo === null && se !== null) lo = b - 1.96 * se;
      if (hi === null && se !== null) hi = b + 1.96 * se;
      const p = num(x.r[cP]);
      byTrait.get(x.trait).push({
        est: b, lo, hi, s: x.mi + 1, sig: p !== null && p < 0.05, row: x.r, mi: x.mi,
        tip: `<b>${esc(x.trait)}</b><br>${METHODS[x.mi].label}<br>b ${fmt(b)} (95% CI ${fmt(lo)} to ${fmt(hi)})<br>p ${fmt(p)}`,
      });
    }
    const groups = [...byTrait].map(([label, items]) => ({ label, items: items.sort((a, b) => a.mi - b.mi) }));
    const ivw = (g) => { const i = g.items.find((x) => x.mi === 0); return i ? i.est : g.items[0] ? g.items[0].est : 0; };
    groups.sort((a, b) => ivw(b) - ivw(a));

    const xl = opt.dir === "forward" ? `Causal estimate of each trait on ${opt.model} (b)` : `Causal estimate of ${opt.model} on each trait (b)`;
    forest(el, {
      title: xl, groups, series: METHODS.map((x) => ({ label: x.label, shape: x.shape })), core: [1, 3],
      nullAt: 0, xLabel: xl, sigLabel: "p < 0.05",
      caption: "Sorted by the inverse-variance-weighted estimate. Estimates are on each trait's own scale, so compare direction and whether the CI crosses zero rather than sizes across traits. The x range is set by the IVW and weighted-median intervals; wider MR Egger intervals are cut off with an arrow. Hover a point for values; click it to see the full row.",
    }, helpers);
  }

  /* ---------- shared hover helper for point/bar charts ---------- */
  function attachTip(el, svg, find, helpers) {
    const tip = document.createElement("div");
    tip.className = "fp-tip"; tip.hidden = true;
    el.appendChild(tip);
    const pt = svg.createSVGPoint();
    const locate = (e) => { pt.x = e.clientX; pt.y = e.clientY; return pt.matrixTransform(svg.getScreenCTM().inverse()); };
    let hl = null;
    const clear = () => { tip.hidden = true; if (hl) hl.remove(); hl = null; svg.style.cursor = ""; };
    svg.addEventListener("mousemove", (e) => {
      const p = locate(e);
      const hit = find(p.x, p.y);
      if (!hit) { clear(); return; }
      tip.innerHTML = hit.tip; tip.hidden = false;
      svg.style.cursor = "pointer";
      if (hl) hl.remove();
      if (hit.ring) { svg.insertAdjacentHTML("beforeend", hit.ring); hl = svg.lastElementChild; }
      const r = el.getBoundingClientRect();
      let left = e.clientX - r.left + 14, top = e.clientY - r.top + 12;
      if (left + tip.offsetWidth > r.width - 8) left = e.clientX - r.left - tip.offsetWidth - 14;
      tip.style.left = left + "px"; tip.style.top = top + "px";
    });
    svg.addEventListener("mouseleave", clear);
    svg.addEventListener("click", (e) => { const p = locate(e); const hit = find(p.x, p.y); if (hit && hit.row) helpers.openDrawer(hit.row); });
  }
  const colG = (t, group, name) => t.columns.findIndex((c) => (c.group || "").toLowerCase().includes(group.toLowerCase()) && c.name.toLowerCase() === name.toLowerCase());
  const pct = (v) => (v * 100).toFixed(v * 100 < 0.1 ? 3 : 2) + "%";

  /* ---------- PGS prediction by ancestry (S29) ---------- */
  function pgs(el, t, rows, st, helpers) {
    last = [el, t, rows, st, helpers, pgs];
    el.innerHTML = "";
    const cT = col(t, "target sample"), cB = col(t, "base gwas");
    const targets = [...new Set(t.columns.filter((c) => c.group).map((c) => c.group))];
    st.chartOpts = st.chartOpts || {};
    const opt = st.chartOpts;
    if (!targets.includes(opt.target)) opt.target = targets[0];
    el.insertAdjacentHTML("beforeend", `<div class="fp-head"><div class="fp-controls">
      <label for="fp-target">Brain age gap measured with <select id="fp-target">${targets.map((g) => `<option ${g === opt.target ? "selected" : ""}>${esc(g.replace(/^Brain age gap \((.*)\)$/, "$1"))}</option>`).join("")}</select></label>
      </div></div>`);
    el.querySelector("#fp-target").addEventListener("change", (e) => { opt.target = targets[e.target.selectedIndex]; pgs(el, t, rows, st, helpers); });
    const cR2 = colG(t, opt.target, "r2"), cRho = colG(t, opt.target, "rho"), cN = colG(t, opt.target, "n"), cP = colG(t, opt.target, "p");

    const anc = [];
    for (const r of rows) {
      const a = String(r[cT]); let g = anc.find((x) => x.name === a);
      if (!g) anc.push(g = { name: a, bars: [] });
      const r2 = num(r[cR2]); if (r2 === null) continue;
      g.bars.push({ base: String(r[cB]), r2, p: num(r[cP]), n: num(r[cN]), rho: num(r[cRho]), row: r });
    }
    const panels = anc.filter((a) => a.bars.length);
    if (!panels.length) { el.insertAdjacentHTML("beforeend", `<div class="fp-empty">No rows to chart. Clear the search or filters.</div>`); return; }
    const bases = []; for (const a of panels) for (const b of a.bars) if (!bases.includes(b.base)) bases.push(b.base);
    const isFactor = (b) => /factor/i.test(b);

    const width = Math.max(320, el.clientWidth - 32);
    const labelW = 112, gap = 18;
    const perRow = width > 980 ? panels.length : width > 600 ? 2 : 1;
    const panelW = Math.floor((width - labelW - gap * (perRow - 1)) / perRow);
    const barH = 14, rowH = 22, headH = 34, axisH = 26;
    const panelH = headH + bases.length * rowH + axisH;
    const nRows = Math.ceil(panels.length / perRow);
    const H = nRows * (panelH + 14);
    let xmax = 0; for (const a of panels) for (const b of a.bars) xmax = Math.max(xmax, b.r2);
    const ticks = niceTicks(0, xmax * 1.08, Math.max(2, Math.round(panelW / 80)));
    xmax = ticks[ticks.length - 1];
    const hits = [];
    let svg = `<svg class="fp-svg" viewBox="0 0 ${width} ${H}" width="${width}" height="${H}" role="img" aria-label="Polygenic score prediction of brain age gap by ancestry">`;
    panels.forEach((a, i) => {
      const row = Math.floor(i / perRow), c = i % perRow;
      const x0 = labelW + c * (panelW + gap), y0 = row * (panelH + 14);
      const X = (v) => x0 + (v / xmax) * (panelW - 8);
      const n = a.bars.find((b) => b.n !== null);
      svg += `<text class="lab" x="${x0}" y="${y0 + 14}" style="font-weight:600">${esc(a.name.replace(/ ancestry/i, ""))}</text>`;
      if (n) svg += `<text class="tick" x="${x0}" y="${y0 + 28}">n = ${n.n.toLocaleString()}</text>`;
      for (const tv of ticks) {
        svg += `<line class="grid" x1="${X(tv)}" x2="${X(tv)}" y1="${y0 + headH - 4}" y2="${y0 + headH + bases.length * rowH}"/>`;
        svg += `<text class="tick" x="${X(tv)}" y="${y0 + headH + bases.length * rowH + 16}" text-anchor="middle">${(tv * 100).toFixed(tv * 100 < 1 && tv > 0 ? 1 : 0)}%</text>`;
      }
      bases.forEach((base, j) => {
        const cy = y0 + headH + j * rowH + rowH / 2;
        if (c === 0) svg += `<text class="lab" x="${labelW - 10}" y="${cy}" text-anchor="end" dominant-baseline="middle" style="${isFactor(base) ? "font-weight:600" : ""}">${esc(base)}</text>`;
        const b = a.bars.find((x) => x.base === base); if (!b) return;
        const s = isFactor(base) ? 2 : 1, sig = b.p !== null && b.p < 0.05;
        const w = Math.max(1.5, X(b.r2) - x0);
        svg += `<rect class="mk s${s} ${sig ? "f" + s : "hollow"}" x="${x0}" y="${cy - barH / 2}" width="${w}" height="${barH}" rx="3"/>`;
        if (isFactor(base)) svg += `<text class="tick" x="${x0 + w + 6}" y="${cy}" dominant-baseline="middle">${pct(b.r2)}</text>`;
        hits.push({ x1: x0, x2: x0 + Math.max(w, 40), y1: cy - rowH / 2, y2: cy + rowH / 2, row: b.row,
          ring: `<rect x="${x0 - 2}" y="${cy - barH / 2 - 2}" width="${w + 4}" height="${barH + 4}" rx="4" fill="none" style="stroke:var(--ink);stroke-width:1.5"/>`,
          tip: `<b>${esc(base)}</b> → ${esc(a.name)}<br>R² ${pct(b.r2)} · rho ${fmt(b.rho)}<br>p ${fmt(b.p)} · n ${b.n === null ? "–" : b.n.toLocaleString()}` });
      });
    });
    svg += `</svg>`;
    el.insertAdjacentHTML("beforeend", `<div class="fp-legend">
      <span>${legendMark("square", 2, false)}BAG factor</span><span>${legendMark("square", 1, false)}Single-model BAG GWAS</span>
      <span>${legendMark("square", 1, false)} filled: p &lt; 0.05</span><span>${legendMark("square", 1, true)} hollow: not significant</span></div>
      <div class="fp-plot">${svg}</div>
      <p class="fp-caption">Variance in brain age gap explained (partial R²) by each polygenic score, per ancestry group of the UK Biobank hold-out sample. All panels share one x-axis. Hover a bar for values; click it to see the full row.</p>`);
    const sv = el.querySelector(".fp-plot svg");
    attachTip(el, sv, (x, y) => hits.find((h) => x >= h.x1 && x <= h.x2 && y >= h.y1 && y <= h.y2), helpers);
  }

  /* ---------- PheWAS (S32) ---------- */
  function phewas(el, t, rows, st, helpers) {
    last = [el, t, rows, st, helpers, phewas];
    el.innerHTML = "";
    const models = [...new Set(t.columns.filter((c) => c.group).map((c) => c.group))];
    st.chartOpts = st.chartOpts || {};
    const opt = st.chartOpts;
    if (!models.includes(opt.model)) opt.model = models.find((m) => /factor/i.test(m)) || models[0];
    if (!opt.cap) opt.cap = 30;
    el.insertAdjacentHTML("beforeend", `<div class="fp-head"><div class="fp-controls">
      <label for="fp-model">Polygenic score <select id="fp-model">${models.map((m) => `<option ${m === opt.model ? "selected" : ""}>${esc(m)}</option>`).join("")}</select></label>
      <label for="fp-cap">Y axis up to <select id="fp-cap">${[10, 30, 60, 0].map((c) => `<option value="${c}" ${c === opt.cap ? "selected" : ""}>${c ? "−log10 p = " + c : "show all"}</option>`).join("")}</select></label>
      </div></div>`);
    el.querySelector("#fp-model").addEventListener("change", (e) => { opt.model = e.target.value; phewas(el, t, rows, st, helpers); });
    el.querySelector("#fp-cap").addEventListener("change", (e) => { opt.cap = +e.target.value; phewas(el, t, rows, st, helpers); });

    const cID = col(t, "phenotype id"), cD = col(t, "phenotype description"), cC = col(t, "category");
    const cP = colG(t, opt.model, "p"), cF = colG(t, opt.model, "fdr"), cRho = colG(t, opt.model, "rho");
    const pts = [];
    for (const r of rows) {
      const p = num(r[cP]); if (p === null || p <= 0) continue;
      pts.push({ r, p, lp: -Math.log10(p), fdr: num(r[cF]), rho: num(r[cRho]), cat: String(r[cC] || "Other"), id: String(r[cID]) });
    }
    if (!pts.length) { el.insertAdjacentHTML("beforeend", `<div class="fp-empty">No rows to chart. Clear the search or filters.</div>`); return; }
    const cats = [...new Set(pts.map((p) => p.cat))].sort((a, b) => a.localeCompare(b));
    pts.sort((a, b) => cats.indexOf(a.cat) - cats.indexOf(b.cat) || a.id.localeCompare(b.id, undefined, { numeric: true }));

    const width = Math.max(320, el.clientWidth - 32);
    const left = 100, right = 12, top = 12, plotH = 360, bottom = 150;
    const H = top + plotH + bottom;
    const maxLp = Math.max(...pts.map((p) => p.lp));
    const yTop = opt.cap ? Math.min(opt.cap, Math.max(maxLp, 2)) : maxLp;
    const yt = niceTicks(0, yTop, 6); const yMax = Math.max(yTop, yt[yt.length - 1]);
    // Each category gets half its share by count and half an equal share, so small categories stay readable.
    const plotW = width - left - right;
    const counts = cats.map((c) => pts.filter((p) => p.cat === c).length);
    const bandW = counts.map((n) => plotW * (0.5 * n / pts.length + 0.5 / cats.length));
    const bandX = []; bandW.reduce((acc, w, k) => (bandX[k] = acc, acc + w), left);
    const posInCat = new Map(); { const seen = {}; pts.forEach((p) => { const k = cats.indexOf(p.cat); seen[k] = (seen[k] || 0) + 1; posInCat.set(p, seen[k] - 1); }); }
    const X = (p) => { const k = cats.indexOf(p.cat); return bandX[k] + ((posInCat.get(p) + 0.5) / counts[k]) * bandW[k]; };
    const Y = (v) => top + plotH - (Math.min(v, yMax) / yMax) * plotH;
    const nSig = pts.filter((p) => p.fdr !== null && p.fdr < 0.05).length;
    const bonf = -Math.log10(0.05 / pts.length);

    let svg = `<svg class="fp-svg" viewBox="0 0 ${width} ${H}" width="${width}" height="${H}" role="img" aria-label="PheWAS of ${esc(opt.model)} polygenic score">`;
    // category bands
    cats.forEach((c, ci) => {
      const xa = bandX[ci], xb = bandX[ci] + bandW[ci];
      if (ci % 2 === 0) svg += `<rect class="band" x="${xa}" y="${top}" width="${xb - xa}" height="${plotH}"/>`;
      svg += `<text class="tick" transform="translate(${(xa + xb) / 2},${top + plotH + 10}) rotate(-40)" text-anchor="end" style="font-family:var(--font-body)">${esc(c)}</text>`;
    });
    for (const v of yt) {
      if (v > yMax) continue;
      svg += `<line class="grid" x1="${left}" x2="${width - right}" y1="${Y(v)}" y2="${Y(v)}"/><text class="tick" x="${left - 6}" y="${Y(v)}" text-anchor="end" dominant-baseline="middle">${v}</text>`;
    }
    svg += `<text class="axt" transform="translate(${left - 40},${top + plotH / 2}) rotate(-90)" text-anchor="middle">−log10(p)</text>`;
    if (bonf < yMax) svg += `<line class="zero" x1="${left}" x2="${width - right}" y1="${Y(bonf)}" y2="${Y(bonf)}"/><text class="tick" x="${width - right - 4}" y="${Y(bonf) - 5}" text-anchor="end">Bonferroni</text>`;
    // non-significant first (muted), then significant on top
    const ns = [], sg = [];
    pts.forEach((p) => { p.x = X(p); p.y = Y(p.lp); p.capped = p.lp > yMax; (p.fdr !== null && p.fdr < 0.05 ? sg : ns).push(p); });
    let dots = "";
    for (const p of ns) dots += `<circle cx="${p.x.toFixed(1)}" cy="${p.y.toFixed(1)}" r="2.2" style="fill:var(--muted);opacity:${cats.indexOf(p.cat) % 2 ? 0.55 : 0.32}"/>`;
    for (const p of sg) {
      const up = (p.rho || 0) >= 0, s = up ? 1 : 2, r = 5;
      const d = up ? `M${p.x},${p.y - r}L${p.x + r},${p.y + r * 0.8}L${p.x - r},${p.y + r * 0.8}Z` : `M${p.x},${p.y + r}L${p.x + r},${p.y - r * 0.8}L${p.x - r},${p.y - r * 0.8}Z`;
      dots += `<path class="mk f${s}" d="${d}" style="stroke:var(--surface);stroke-width:1"/>`;
      if (p.capped) dots += `<text class="tick" x="${p.x}" y="${p.y - 8}" text-anchor="middle" style="font-size:9px">↑</text>`;
    }
    svg += dots + `</svg>`;

    el.insertAdjacentHTML("beforeend", `<div class="fp-legend">
      <span>${legendTri(true)} FDR &lt; 0.05, positive association</span><span>${legendTri(false)} FDR &lt; 0.05, negative association</span>
      <span><svg width="10" height="10" aria-hidden="true"><circle cx="5" cy="5" r="2.5" style="fill:var(--muted);opacity:.5"/></svg> not FDR-significant</span>
      <span>${pts.length.toLocaleString()} phenotypes · ${nSig.toLocaleString()} FDR-significant</span></div>
      <div class="fp-plot">${svg}</div>
      <p class="fp-caption">Each point is one UK Biobank phenotype, grouped by category. Direction follows the sign of rho. ${opt.cap ? `Points above −log10 p = ${yMax} are drawn at the top edge with an arrow; choose “show all” to see them to scale.` : ""} Hover for the phenotype and statistics; click to see the full row. Use the table search (for example “blood pressure”) to chart a subset.</p>`);

    const sv = el.querySelector(".fp-plot svg");
    const all = sg.concat(ns);
    attachTip(el, sv, (x, y) => {
      let best = null, bd = 64;
      for (const p of all) { const d = (p.x - x) ** 2 + (p.y - y) ** 2; if (d < bd) { bd = d; best = p; } }
      if (!best) return null;
      return { row: best.r, ring: `<circle cx="${best.x}" cy="${best.y}" r="7" fill="none" style="stroke:var(--ink);stroke-width:1.5"/>`,
        tip: `<b>${esc(best.r[cD])}</b><br>${esc(best.cat)}<br>rho ${fmt(best.rho)} · p ${fmt(best.p)}${best.capped ? " (above axis)" : ""}<br>FDR ${fmt(best.fdr)}` };
    }, helpers);
  }
  function legendTri(up) {
    const d = up ? "M7,2L12,11L2,11Z" : "M7,12L12,3L2,3Z";
    return `<svg class="fp-svg" width="14" height="14" viewBox="0 0 14 14" aria-hidden="true" style="width:14px;display:inline"><path class="f${up ? 1 : 2}" d="${d}"/></svg>`;
  }

  window.SuppCharts = {
    S13: { render: ldsc },
    S20: { render: ldsc },
    S22: { render: mr },
    S23: { render: mr },
    S24: { render: mr },
    S29: { render: pgs },
    S32: { render: phewas },
  };
})();
