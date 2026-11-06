// Built-in plots for kaks batch results: histograms of dN, dS, dN/dS and a
// dN-vs-dS scatter, drawn as SVG strings (no DOM needed, so the pure parts run
// under Node for tests). No third-party code.

/** Column names for each metric source. */
export const SOURCES = {
  yn00: { label: "YN00", dN: "dN", dS: "dS", omega: "omega" },
  ng86: { label: "NG86", dN: "NG_dN", dS: "NG_dS", omega: "NG_omega" },
  ml: { label: "ML (codeml-style)", dN: "ML_dN", dS: "ML_dS", omega: "ML_omega" },
};

/** Which sources are present in a set of result rows. */
export function availableSources(rows) {
  const r = rows[0] || {};
  return Object.keys(SOURCES).filter((k) => SOURCES[k].dN in r);
}

/**
 * Apply the filters. Returns { kept, dropped, reasons } where reasons counts
 * rows removed by each rule. A row must have finite dN, dS and omega.
 *   opts.dsMax    drop rows with dS > dsMax (saturation), null to skip
 *   opts.iqrK     drop rows whose dN, dS or omega lies outside
 *                 [Q1 - k*IQR, Q3 + k*IQR], null to skip
 *   opts.omegaMax drop rows with omega > omegaMax, null to skip
 */
export function filterRows(rows, src, opts = {}) {
  const c = SOURCES[src];
  const reasons = { nonfinite: 0, dsMax: 0, omegaMax: 0, iqr: 0 };
  let kept = rows.filter((r) => {
    const ok = [c.dN, c.dS, c.omega].every((k) => Number.isFinite(r[k]));
    if (!ok) reasons.nonfinite++;
    return ok;
  });
  if (opts.dsMax != null) {
    const n = kept.length;
    kept = kept.filter((r) => r[c.dS] <= opts.dsMax);
    reasons.dsMax = n - kept.length;
  }
  if (opts.omegaMax != null) {
    const n = kept.length;
    kept = kept.filter((r) => r[c.omega] <= opts.omegaMax);
    reasons.omegaMax = n - kept.length;
  }
  if (opts.iqrK != null && kept.length >= 4) {
    const bounds = {};
    for (const k of [c.dN, c.dS, c.omega]) {
      const q = quantiles(kept.map((r) => r[k]));
      const iqr = q.q3 - q.q1;
      bounds[k] = [q.q1 - opts.iqrK * iqr, q.q3 + opts.iqrK * iqr];
    }
    const n = kept.length;
    kept = kept.filter((r) => [c.dN, c.dS, c.omega].every((k) => r[k] >= bounds[k][0] && r[k] <= bounds[k][1]));
    reasons.iqr = n - kept.length;
  }
  return { kept, dropped: rows.length - kept.length, reasons };
}

/** Sorted quartiles (linear interpolation), mean, min, max. */
export function quantiles(values) {
  const v = values.filter(Number.isFinite).sort((a, b) => a - b);
  const n = v.length;
  if (!n) return { n: 0 };
  const at = (p) => {
    const i = (n - 1) * p;
    const lo = Math.floor(i), hi = Math.ceil(i);
    return v[lo] + (v[hi] - v[lo]) * (i - lo);
  };
  return {
    n,
    min: v[0],
    max: v[n - 1],
    q1: at(0.25),
    median: at(0.5),
    q3: at(0.75),
    mean: v.reduce((s, x) => s + x, 0) / n,
  };
}

/** Equal-width bins over [lo, hi]; values outside are ignored. */
export function histogram(values, bins, lo, hi) {
  const counts = new Array(bins).fill(0);
  const w = (hi - lo) / bins;
  for (const x of values) {
    if (!Number.isFinite(x) || x < lo || x > hi) continue;
    let b = Math.floor((x - lo) / w);
    if (b === bins) b = bins - 1; // x === hi
    counts[b]++;
  }
  return { counts, lo, hi, width: w };
}

const PALETTE = ["#1d4ed8", "#dc2626", "#059669", "#d97706", "#7c3aed", "#0891b2", "#be185d", "#4d7c0f", "#b45309", "#6d28d9", "#0f766e", "#9f1239"];

/** Stable colour per category value. */
export function colorScale(categories) {
  const map = new Map();
  [...new Set(categories)].forEach((c, i) => map.set(c, PALETTE[i % PALETTE.length]));
  return map;
}

const esc = (s) => String(s).replace(/[&<>"]/g, (ch) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[ch]));
const fmt = (x) => (Math.abs(x) >= 100 ? x.toFixed(1) : Math.abs(x) >= 1 ? x.toFixed(3) : x.toPrecision(3));

/** "Nice" tick positions for an axis. */
export function ticks(lo, hi, n = 5) {
  if (!(hi > lo)) return [lo];
  const raw = (hi - lo) / n;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => (hi - lo) / s <= n + 1) || mag;
  const out = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + 1e-9; t += step) out.push(+t.toFixed(10));
  return out;
}

/**
 * Histogram SVG. `series` is [{name, values, color}]; with several series the
 * bars are stacked. opts: {bins, log, xmin, xmax, width, height, title, xlabel}
 */
export function histogramSvg(series, opts = {}) {
  const bins = opts.bins || 20;
  const log = !!opts.log;
  const W = opts.width || 360, H = opts.height || 240;
  const m = { l: 48, r: 12, t: 28, b: 42 };
  const tx = (v) => (log ? Math.log10(v) : v);
  const all = series.flatMap((s) => s.values).filter((v) => Number.isFinite(v) && (!log || v > 0));
  if (!all.length) return `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}"><text x="12" y="24" font-size="12">${esc(opts.title || "")}: no data</text></svg>`;
  let lo = opts.xmin != null ? tx(opts.xmin) : Math.min(...all.map(tx));
  let hi = opts.xmax != null ? tx(opts.xmax) : Math.max(...all.map(tx));
  if (hi <= lo) hi = lo + 1;
  const hists = series.map((s) => histogram(s.values.filter((v) => !log || v > 0).map(tx), bins, lo, hi));
  const totals = new Array(bins).fill(0);
  hists.forEach((h) => h.counts.forEach((c, i) => (totals[i] += c)));
  const ymax = Math.max(1, ...totals);
  const px = (v) => m.l + ((v - lo) / (hi - lo)) * (W - m.l - m.r);
  const py = (c) => H - m.b - (c / ymax) * (H - m.t - m.b);
  let out = `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="system-ui,sans-serif" font-size="11">`;
  out += `<rect width="${W}" height="${H}" fill="white"/>`;
  out += `<text x="${m.l}" y="18" font-size="13" font-weight="600">${esc(opts.title || "")}</text>`;
  // bars
  const base = new Array(bins).fill(0);
  hists.forEach((h, si) => {
    for (let i = 0; i < bins; i++) {
      const c = h.counts[i];
      if (!c) continue;
      const x0 = px(lo + i * h.width), x1 = px(lo + (i + 1) * h.width);
      const y1 = py(base[i]), y0 = py(base[i] + c);
      out += `<rect x="${x0.toFixed(1)}" y="${y0.toFixed(1)}" width="${Math.max(0.5, x1 - x0 - 1).toFixed(1)}" height="${(y1 - y0).toFixed(1)}" fill="${series[si].color}"><title>${esc(series[si].name)}: ${c}</title></rect>`;
      base[i] += c;
    }
  });
  // axes
  out += `<line x1="${m.l}" y1="${H - m.b}" x2="${W - m.r}" y2="${H - m.b}" stroke="#444"/>`;
  out += `<line x1="${m.l}" y1="${m.t}" x2="${m.l}" y2="${H - m.b}" stroke="#444"/>`;
  for (const t of ticks(lo, hi, 5)) {
    const x = px(t);
    const label = log ? (Number.isInteger(t) ? `1e${t}` : fmt(10 ** t)) : fmt(t);
    out += `<line x1="${x.toFixed(1)}" y1="${H - m.b}" x2="${x.toFixed(1)}" y2="${H - m.b + 4}" stroke="#444"/><text x="${x.toFixed(1)}" y="${H - m.b + 15}" text-anchor="middle">${label}</text>`;
  }
  for (const t of ticks(0, ymax, 4)) {
    const y = py(t);
    out += `<line x1="${m.l - 4}" y1="${y.toFixed(1)}" x2="${m.l}" y2="${y.toFixed(1)}" stroke="#444"/><text x="${m.l - 6}" y="${(y + 3.5).toFixed(1)}" text-anchor="end">${t}</text>`;
  }
  out += `<text x="${(m.l + W - m.r) / 2}" y="${H - 8}" text-anchor="middle">${esc(opts.xlabel || "")}${log ? " (log10)" : ""}</text>`;
  out += `<text transform="translate(12,${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle">pairs</text>`;
  out += legendSvg(series, W - m.r, m.t + 4);
  return out + "</svg>";
}

function legendSvg(series, xRight, yTop) {
  if (series.length < 2 || series.length > 12) return "";
  let out = "";
  series.forEach((s, i) => {
    const y = yTop + i * 13;
    out += `<rect x="${xRight - 110}" y="${y}" width="9" height="9" fill="${s.color}"/><text x="${xRight - 97}" y="${y + 8}">${esc(String(s.name).slice(0, 18))}</text>`;
  });
  return out;
}

/**
 * Scatter SVG of y against x. `points` is [{x, y, label, color}].
 * opts: {logx, logy, width, height, title, xlabel, ylabel, legend:[{name,color}]}
 */
export function scatterSvg(points, opts = {}) {
  const W = opts.width || 420, H = opts.height || 320;
  const m = { l: 52, r: 12, t: 28, b: 42 };
  const fx = (v) => (opts.logx ? Math.log10(v) : v);
  const fy = (v) => (opts.logy ? Math.log10(v) : v);
  const pts = points.filter((p) => Number.isFinite(p.x) && Number.isFinite(p.y) && (!opts.logx || p.x > 0) && (!opts.logy || p.y > 0));
  if (!pts.length) return `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}"><text x="12" y="24" font-size="12">${esc(opts.title || "")}: no data</text></svg>`;
  const xs = pts.map((p) => fx(p.x)), ys = pts.map((p) => fy(p.y));
  let xlo = Math.min(...xs), xhi = Math.max(...xs), ylo = Math.min(...ys), yhi = Math.max(...ys);
  if (xhi <= xlo) xhi = xlo + 1;
  if (yhi <= ylo) yhi = ylo + 1;
  const padx = (xhi - xlo) * 0.04, pady = (yhi - ylo) * 0.04;
  xlo -= padx; xhi += padx; ylo -= pady; yhi += pady;
  const px = (v) => m.l + ((v - xlo) / (xhi - xlo)) * (W - m.l - m.r);
  const py = (v) => H - m.b - ((v - ylo) / (yhi - ylo)) * (H - m.t - m.b);
  const tickLabel = (t, log) => (log ? (Number.isInteger(t) ? `1e${t}` : fmt(10 ** t)) : fmt(t));
  let out = `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="system-ui,sans-serif" font-size="11">`;
  out += `<rect width="${W}" height="${H}" fill="white"/>`;
  out += `<text x="${m.l}" y="18" font-size="13" font-weight="600">${esc(opts.title || "")}</text>`;
  for (const t of ticks(xlo, xhi, 5)) {
    const x = px(t);
    out += `<line x1="${x.toFixed(1)}" y1="${m.t}" x2="${x.toFixed(1)}" y2="${H - m.b}" stroke="#eee"/><text x="${x.toFixed(1)}" y="${H - m.b + 15}" text-anchor="middle">${tickLabel(t, opts.logx)}</text>`;
  }
  for (const t of ticks(ylo, yhi, 5)) {
    const y = py(t);
    out += `<line x1="${m.l}" y1="${y.toFixed(1)}" x2="${W - m.r}" y2="${y.toFixed(1)}" stroke="#eee"/><text x="${m.l - 6}" y="${(y + 3.5).toFixed(1)}" text-anchor="end">${tickLabel(t, opts.logy)}</text>`;
  }
  out += `<line x1="${m.l}" y1="${H - m.b}" x2="${W - m.r}" y2="${H - m.b}" stroke="#444"/><line x1="${m.l}" y1="${m.t}" x2="${m.l}" y2="${H - m.b}" stroke="#444"/>`;
  for (const p of pts) {
    out += `<circle cx="${px(fx(p.x)).toFixed(1)}" cy="${py(fy(p.y)).toFixed(1)}" r="3.5" fill="${p.color || "#1d4ed8"}" fill-opacity="0.7"><title>${esc(p.label || "")}: x=${fmt(p.x)} y=${fmt(p.y)}</title></circle>`;
  }
  out += `<text x="${(m.l + W - m.r) / 2}" y="${H - 8}" text-anchor="middle">${esc(opts.xlabel || "x")}${opts.logx ? " (log10)" : ""}</text>`;
  out += `<text transform="translate(12,${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle">${esc(opts.ylabel || "y")}${opts.logy ? " (log10)" : ""}</text>`;
  if (opts.legend) out += legendSvg(opts.legend, W - m.r, m.t + 4);
  return out + "</svg>";
}

/** One-line summary text for a metric. */
export function summaryText(values) {
  const q = quantiles(values);
  if (!q.n) return "no values";
  return `n=${q.n}  median ${fmt(q.median)}  mean ${fmt(q.mean)}  Q1 ${fmt(q.q1)}  Q3 ${fmt(q.q3)}  min ${fmt(q.min)}  max ${fmt(q.max)}`;
}

/**
 * Build everything for the Plots section from result rows.
 * state: {src, bins, log, xmin, xmax, dsMax, iqrK, omegaMax, colorBy}
 * Returns {html, kept, dropped, reasons, svgs:{dN,dS,omega,scatter}}
 */
export function buildPlots(rows, state) {
  const src = SOURCES[state.src];
  const { kept, dropped, reasons } = filterRows(rows, state.src, {
    dsMax: state.dsMax, iqrK: state.iqrK, omegaMax: state.omegaMax,
  });
  const key = state.colorBy === "file" ? "file" : state.colorBy === "code" ? "genetic_code" : null;
  const cats = key ? kept.map((r) => r[key] ?? "?") : [];
  const scale = key ? colorScale(cats) : null;
  const groups = key
    ? [...scale.keys()].map((name) => ({ name, color: scale.get(name), rows: kept.filter((r) => (r[key] ?? "?") === name) }))
    : [{ name: "all", color: PALETTE[0], rows: kept }];
  const common = { bins: state.bins, log: state.log, xmin: state.xmin, xmax: state.xmax };
  const hist = (col, title) =>
    histogramSvg(groups.map((g) => ({ name: g.name, color: g.color, values: g.rows.map((r) => r[col]) })), { ...common, title, xlabel: title });
  const svgs = {
    dN: hist(src.dN, "dN"),
    dS: hist(src.dS, "dS"),
    omega: hist(src.omega, "dN/dS"),
    scatter: scatterSvg(
      groups.flatMap((g) => g.rows.map((r) => ({ x: r[src.dS], y: r[src.dN], color: g.color, label: `${r.file ? r.file + " " : ""}${r.seq1}-${r.seq2}` }))),
      { logx: state.log, logy: state.log, title: `dN vs dS (${src.label})`, xlabel: "dS", ylabel: "dN", legend: key ? groups : null },
    ),
  };
  const summaries = {
    dN: summaryText(kept.map((r) => r[src.dN])),
    dS: summaryText(kept.map((r) => r[src.dS])),
    omega: summaryText(kept.map((r) => r[src.omega])),
  };
  return { kept, dropped, reasons, svgs, summaries };
}
