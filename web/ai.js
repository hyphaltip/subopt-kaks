// Optional "ask a model for a plot" panel. The user's API key is held only in
// memory; the request goes from the browser straight to the provider. What is
// sent: column names, a few sample rows, summary statistics and the request
// text. Sequences are never sent. The model answers with a Vega-Lite spec,
// which the page renders with vega-embed (loaded from a CDN only when this
// panel is first used).

export const PROVIDERS = {
  openai: {
    label: "OpenAI-compatible (OpenAI, Anthropic compat, Ollama, vLLM, ...)",
    endpoint: "https://api.openai.com/v1/chat/completions",
    model: "gpt-4o-mini",
  },
  gemini: {
    label: "Google Gemini",
    endpoint: "https://generativelanguage.googleapis.com/v1beta",
    model: "gemini-2.0-flash",
  },
};

const SYSTEM = `You write Vega-Lite v5 specifications for a dataset of pairwise dN/dS results.
Reply with ONLY one JSON object, the Vega-Lite spec, no prose and no code fences.
Use "data": {"name": "table"} exactly; the page supplies the rows.
Keep specs self-contained: "width" and "height" set, axis titles set, tooltips on.
For filtering use a "transform" with "filter" expressions on the column names given.
Use "scale": {"type": "log"} on an axis when asked for log axes; never log-transform values yourself.
Prefer "mark": "point" for scatter plots and "mark": "bar" with "bin": true for histograms.
Column meanings: dN, dS, omega are the YN00 counting estimates; NG_* are Nei-Gojobori 1986;
ML_* are maximum-likelihood (codeml-style) estimates, lnL the log-likelihood; file is the input file name;
seq1 and seq2 are the sequence names; genetic_code is the NCBI translation table id.`;

/** Describe the table for the model: columns, types, 5 sample rows, summaries. */
export function describeTable(rows) {
  const cols = Object.keys(rows[0] || {});
  const numeric = cols.filter((c) => typeof rows[0][c] === "number");
  const stats = numeric.map((c) => {
    const v = rows.map((r) => r[c]).filter(Number.isFinite).sort((a, b) => a - b);
    if (!v.length) return `${c}: no finite values`;
    return `${c}: n=${v.length} min=${v[0]} median=${v[Math.floor(v.length / 2)]} max=${v[v.length - 1]}`;
  });
  return [
    `Rows: ${rows.length}. Columns: ${cols.join(", ")}.`,
    `Numeric summaries:\n${stats.join("\n")}`,
    `Sample rows (JSON):\n${JSON.stringify(rows.slice(0, 5))}`,
  ].join("\n\n");
}

/** Build the fetch request for a provider. Exported for tests; no network here. */
export function buildRequest(provider, { endpoint, key, model }, userText, tableText, previousSpec) {
  const user = `${tableText}\n\nRequest: ${userText}` + (previousSpec ? `\n\nCurrent spec to modify:\n${previousSpec}` : "");
  if (provider === "gemini") {
    const base = endpoint.replace(/\/+$/, "");
    return {
      url: `${base}/models/${encodeURIComponent(model)}:generateContent?key=${encodeURIComponent(key)}`,
      init: {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          system_instruction: { parts: [{ text: SYSTEM }] },
          contents: [{ role: "user", parts: [{ text: user }] }],
          generationConfig: { temperature: 0.2, responseMimeType: "application/json" },
        }),
      },
      extract: (j) => j?.candidates?.[0]?.content?.parts?.map((p) => p.text).join("") ?? "",
    };
  }
  return {
    url: endpoint,
    init: {
      method: "POST",
      headers: { "Content-Type": "application/json", Authorization: `Bearer ${key}` },
      body: JSON.stringify({
        model,
        temperature: 0.2,
        messages: [
          { role: "system", content: SYSTEM },
          { role: "user", content: user },
        ],
      }),
    },
    extract: (j) => j?.choices?.[0]?.message?.content ?? "",
  };
}

/** Pull the JSON object out of a model reply (tolerates code fences and prose). */
export function extractSpec(text) {
  let t = String(text).trim();
  const fence = t.match(/```(?:json)?\s*([\s\S]*?)```/i);
  if (fence) t = fence[1].trim();
  const start = t.indexOf("{"), end = t.lastIndexOf("}");
  if (start < 0 || end <= start) throw new Error("the reply contains no JSON object");
  const spec = JSON.parse(t.slice(start, end + 1));
  if (typeof spec !== "object" || spec === null) throw new Error("the reply is not a Vega-Lite spec");
  if (!spec.$schema) spec.$schema = "https://vega.github.io/schema/vega-lite/v5.json";
  spec.data = { name: "table" };
  return spec;
}

/** Ask the model; returns the parsed spec. */
export async function askForSpec(provider, creds, userText, rows, previousSpec) {
  if (!creds.key) throw new Error("enter an API key");
  if (!rows.length) throw new Error("run an analysis first");
  const req = buildRequest(provider, creds, userText, describeTable(rows), previousSpec);
  const r = await fetch(req.url, req.init);
  if (!r.ok) {
    const body = await r.text().catch(() => "");
    throw new Error(`${provider} returned HTTP ${r.status}: ${body.slice(0, 300)}`);
  }
  const text = req.extract(await r.json());
  if (!text) throw new Error("empty reply from the model");
  return extractSpec(text);
}

let vegaLoading = null;
/** Load vega-embed (and vega, vega-lite) from cdnjs once, on first use. */
export function loadVega() {
  if (globalThis.vegaEmbed) return Promise.resolve(globalThis.vegaEmbed);
  if (vegaLoading) return vegaLoading;
  const add = (src) =>
    new Promise((res, rej) => {
      const s = document.createElement("script");
      s.src = src;
      s.onload = res;
      s.onerror = () => rej(new Error(`could not load ${src}`));
      document.head.appendChild(s);
    });
  vegaLoading = add("https://cdnjs.cloudflare.com/ajax/libs/vega/5.30.0/vega.min.js")
    .then(() => add("https://cdnjs.cloudflare.com/ajax/libs/vega-lite/5.21.0/vega-lite.min.js"))
    .then(() => add("https://cdnjs.cloudflare.com/ajax/libs/vega-embed/6.26.0/vega-embed.min.js"))
    .then(() => globalThis.vegaEmbed);
  return vegaLoading;
}

/** Render a spec into `el` with the result rows as the "table" dataset. */
export async function renderSpec(el, spec, rows) {
  const embed = await loadVega();
  const s = { ...spec, data: { values: rows } };
  return embed(el, s, { actions: { export: true, source: true, editor: false, compiled: false } });
}
