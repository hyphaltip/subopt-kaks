//! Substitution-matrix loading (NCBI/BLAST text format) with a dynamically
//! settable matrix directory and an embedded BLOSUM62 default.
//!
//! Resolution order for a matrix given as a bare name (e.g. `BLOSUM62`):
//!   1. `--matrix-dir` argument
//!   2. `$SUBOPT_MATRIX_DIR`
//!   3. `$SUBOPTDIR/matricies/aa`, `$SUBOPT_DIR/matricies/aa`
//!   4. `./matricies/aa`, `../matricies/aa`, `../../matricies/aa`
//!   5. the compiled-in embedded BLOSUM62 (only for the default name)
//!
//! A value containing a path separator is treated as an explicit file path.

use std::path::{Path, PathBuf};

/// The repo's own BLOSUM62, embedded at compile time so the tool stays
/// self-contained and zero-config when no matrix directory is available.
const BLOSUM62_TEXT: &str = include_str!("../../../matricies/aa/BLOSUM62");

/// A residue substitution matrix indexed by amino-acid letter (ASCII).
pub struct ScoreMatrix {
    /// `table[a][b]` for ASCII residue codes; unknown pairs fall back to `default`.
    table: Box<[[i32; 128]; 128]>,
    default: i32,
}

impl ScoreMatrix {
    /// Score for a residue pair (case-insensitive); unknown residues use the
    /// matrix's own gap/`*` score.
    #[inline]
    pub fn score(&self, a: u8, b: u8) -> i32 {
        let a = a.to_ascii_uppercase() as usize;
        let b = b.to_ascii_uppercase() as usize;
        if a < 128 && b < 128 {
            self.table[a][b]
        } else {
            self.default
        }
    }

    /// Parse a matrix in NCBI/BLAST text format (`#` comments, a header row of
    /// residue letters, then one score row per residue).
    pub fn parse(text: &str) -> Result<ScoreMatrix, String> {
        let mut header: Vec<u8> = Vec::new();
        let mut table = Box::new([[i32::MIN; 128]; 128]);

        for line in text.lines() {
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') {
                continue;
            }
            let toks: Vec<&str> = t.split_whitespace().collect();
            if header.is_empty() {
                header = toks.iter().map(|s| s.as_bytes()[0]).collect();
                continue;
            }
            let row = toks[0].as_bytes()[0].to_ascii_uppercase() as usize;
            for (k, tok) in toks[1..].iter().enumerate() {
                if k >= header.len() {
                    break;
                }
                let col = header[k].to_ascii_uppercase() as usize;
                let v: i32 = tok
                    .parse()
                    .map_err(|_| format!("bad matrix value '{tok}'"))?;
                table[row][col] = v;
            }
        }
        if header.is_empty() {
            return Err("no matrix header row found".into());
        }

        // Default (unknown/stop) score: the `*`×`*` cell if present, else -4.
        let star = table[b'*' as usize][b'*' as usize];
        let default = if star == i32::MIN { -4 } else { star };
        for i in 0..128 {
            for j in 0..128 {
                if table[i][j] == i32::MIN {
                    table[i][j] = default;
                }
            }
        }
        Ok(ScoreMatrix { table, default })
    }

    /// The embedded BLOSUM62 default.
    pub fn embedded_blosum62() -> ScoreMatrix {
        ScoreMatrix::parse(BLOSUM62_TEXT).expect("embedded BLOSUM62 must parse")
    }
}

/// Directories to search for a named matrix, in priority order.
fn search_dirs(matrix_dir: Option<&str>) -> Vec<PathBuf> {
    let mut dirs: Vec<PathBuf> = Vec::new();
    if let Some(d) = matrix_dir {
        dirs.push(PathBuf::from(d));
    }
    if let Ok(d) = std::env::var("SUBOPT_MATRIX_DIR") {
        if !d.is_empty() {
            dirs.push(PathBuf::from(d));
        }
    }
    for var in ["SUBOPTDIR", "SUBOPT_DIR"] {
        if let Ok(d) = std::env::var(var) {
            if !d.is_empty() {
                dirs.push(Path::new(&d).join("matricies").join("aa"));
            }
        }
    }
    for rel in ["matricies/aa", "../matricies/aa", "../../matricies/aa"] {
        dirs.push(PathBuf::from(rel));
    }
    dirs
}

/// Load the scoring matrix named/located by `spec`, searching `matrix_dir` and
/// the environment. `None` uses the default (BLOSUM62), which falls back to the
/// embedded copy if no matrix directory is present.
pub fn load(spec: Option<&str>, matrix_dir: Option<&str>) -> Result<ScoreMatrix, String> {
    let name = spec.unwrap_or("BLOSUM62");

    // Explicit path: read verbatim.
    if name.contains('/') || name.contains('\\') {
        let txt = std::fs::read_to_string(name)
            .map_err(|e| format!("cannot read matrix file '{name}': {e}"))?;
        return ScoreMatrix::parse(&txt);
    }

    // Bare name: search the matrix directories.
    for dir in search_dirs(matrix_dir) {
        let p = dir.join(name);
        if p.is_file() {
            let txt = std::fs::read_to_string(&p)
                .map_err(|e| format!("cannot read matrix file '{}': {e}", p.display()))?;
            return ScoreMatrix::parse(&txt);
        }
    }

    // Fall back to the embedded default only for the default matrix name.
    if name.eq_ignore_ascii_case("BLOSUM62") {
        return Ok(ScoreMatrix::embedded_blosum62());
    }

    Err(format!(
        "cannot locate substitution matrix '{name}' in any matrix directory\n\
         (use --matrix <file>, --matrix-dir <dir>, or set $SUBOPT_MATRIX_DIR)"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_matches_known_scores() {
        let m = ScoreMatrix::embedded_blosum62();
        assert_eq!(m.score(b'A', b'A'), 4);
        assert_eq!(m.score(b'W', b'W'), 11);
        assert_eq!(m.score(b'A', b'R'), -1);
        // case-insensitive
        assert_eq!(m.score(b'a', b'r'), -1);
    }
}
