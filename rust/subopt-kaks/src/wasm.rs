//! Browser entry points (`--features wasm`, target `wasm32-unknown-unknown`).
//!
//! Everything the page needs is embedded: BLOSUM62 and all genetic codes. No
//! file or environment access happens in this module.

use crate::{align, analyze, genetic_code, matrix, report, CodonFreq, Method, Options, OutFormat};
use wasm_bindgen::prelude::*;

/// Run `kaks` on sequence text and return the rendered report.
///
/// * `method`: yn00 | ng86 | both | ml | all
/// * `codon_freq`: fequal | f1x4 | f3x4 | fcodon
/// * `aln_mode`: global | local
/// * `format`: tsv | table | json
/// * `matrix_text`: a PAML/BLAST-style scoring matrix, or empty for BLOSUM62
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn kaks_run(
    input: &str,
    prealigned: bool,
    method: &str,
    codon_freq: &str,
    genetic_code: u8,
    aln_mode: &str,
    gap_open: i32,
    gap_ext: i32,
    common_f3x4: bool,
    format: &str,
    matrix_text: &str,
) -> Result<String, JsError> {
    let method = Method::parse(method).ok_or_else(|| JsError::new("unknown method"))?;
    let format = OutFormat::parse(format).ok_or_else(|| JsError::new("unknown output format"))?;
    let codon_freq = match codon_freq.to_ascii_lowercase().as_str() {
        "fequal" => CodonFreq::Fequal,
        "f1x4" => CodonFreq::F1x4,
        "f3x4" => CodonFreq::F3x4,
        "fcodon" => CodonFreq::Fcodon,
        _ => return Err(JsError::new("unknown codon-frequency model")),
    };
    let aln_mode = match aln_mode.to_ascii_lowercase().as_str() {
        "global" => align::Mode::Global,
        "local" => align::Mode::Local,
        _ => return Err(JsError::new("unknown alignment mode")),
    };
    let matrix = if matrix_text.trim().is_empty() {
        None
    } else {
        Some(matrix::ScoreMatrix::parse(matrix_text).map_err(|e| JsError::new(&e))?)
    };
    let opts = Options {
        prealigned,
        informat: None,
        method,
        common_f3x4,
        codon_freq,
        genetic_code,
        aln_mode,
        matrix,
        gap_open,
        gap_ext,
    };
    let results = analyze(input, &opts).map_err(|e| JsError::new(&e))?;
    Ok(report::to_string(&results, method, format, true))
}

/// All supported genetic codes, one `id<TAB>name` per line.
#[wasm_bindgen]
pub fn kaks_codes() -> String {
    genetic_code::list()
}

/// Crate version.
#[wasm_bindgen]
pub fn kaks_version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}
