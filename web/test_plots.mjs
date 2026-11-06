import { filterRows, quantiles, histogram, ticks, histogramSvg, scatterSvg, buildPlots, availableSources } from "./plots.js";
import { describeTable, buildRequest, extractSpec } from "./ai.js";
import assert from "node:assert/strict";
const rows = [];
for (let i = 0; i < 50; i++) rows.push({ file: `f${i % 3}.fa`, seq1: "a", seq2: "b", genetic_code: 1, dN: 0.01 * (i + 1), dS: 0.1 * (i + 1), omega: 0.1, NG_dN: 0.02, NG_dS: 0.2, NG_omega: 0.1 });
rows.push({ file: "f9.fa", seq1: "a", seq2: "b", genetic_code: 1, dN: 0.5, dS: 40, omega: 0.0125 }); // saturated
rows.push({ file: "f9.fa", seq1: "a", seq2: "b", genetic_code: 1, dN: NaN, dS: 1, omega: NaN });     // undefined
assert.deepEqual(availableSources(rows), ["yn00", "ng86"]);
const f = filterRows(rows, "yn00", { dsMax: 2 });
assert.equal(f.reasons.nonfinite, 1); assert.equal(f.reasons.dsMax, 31); assert.equal(f.kept.length, 20);
const q = quantiles([1, 2, 3, 4, 5]); assert.equal(q.median, 3); assert.equal(q.q1, 2); assert.equal(q.q3, 4); assert.equal(q.mean, 3);
const h = histogram([0, 0.5, 1, 1, 2], 4, 0, 2); assert.deepEqual(h.counts, [1, 1, 2, 1]);
assert.deepEqual(ticks(0, 1, 5), [0, 0.2, 0.4, 0.6, 0.8, 1]);
const iqr = filterRows(rows, "yn00", { iqrK: 1.5 }); assert.ok(iqr.reasons.iqr >= 1);
const svg = histogramSvg([{ name: "all", color: "#000", values: rows.map((r) => r.dS) }], { bins: 10, title: "dS", log: true });
assert.ok(svg.startsWith("<svg") && svg.endsWith("</svg>") && svg.includes("<rect"));
const sc = scatterSvg(rows.map((r) => ({ x: r.dS, y: r.dN })), { title: "t" }); assert.ok((sc.match(/<circle/g) || []).length === 51);
const bp = buildPlots(rows, { src: "yn00", bins: 15, log: false, dsMax: 2, iqrK: null, omegaMax: null, colorBy: "file" });
assert.equal(bp.kept.length, 20); assert.ok(bp.svgs.dN.includes("</svg>") && bp.summaries.dS.startsWith("n=20"));
// xss-safe labels
const bad = scatterSvg([{ x: 1, y: 1, label: "<script>" }], {}); assert.ok(!bad.includes("<script>"));
// ai helpers
const d = describeTable(rows); assert.ok(d.includes("Columns: file, seq1") && d.includes("dS: n=52"));
const g = buildRequest("gemini", { endpoint: "https://generativelanguage.googleapis.com/v1beta/", key: "K", model: "gemini-2.0-flash" }, "plot", d, null);
assert.ok(g.url.includes("/models/gemini-2.0-flash:generateContent?key=K")); assert.ok(JSON.parse(g.init.body).system_instruction);
const o = buildRequest("openai", { endpoint: "https://api.openai.com/v1/chat/completions", key: "K", model: "gpt-4o-mini" }, "plot", d, "{}");
assert.equal(o.init.headers.Authorization, "Bearer K"); assert.equal(JSON.parse(o.init.body).messages[1].role, "user");
assert.equal(o.extract({ choices: [{ message: { content: "x" } }] }), "x");
assert.equal(g.extract({ candidates: [{ content: { parts: [{ text: "y" }] } }] }), "y");
const spec = extractSpec('Here you go:\n```json\n{"mark":"point","encoding":{"x":{"field":"dS"}}}\n```');
assert.equal(spec.mark, "point"); assert.deepEqual(spec.data, { name: "table" }); assert.ok(spec.$schema.includes("vega-lite"));
assert.throws(() => extractSpec("no json here"));
console.log("plots.js + ai.js: all assertions passed");
