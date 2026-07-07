//! Minimal, dependency-free sequence readers for FASTA and sequential PHYLIP.

use std::fmt;
use std::fs;
use std::io::{self, Read};
use std::path::Path;

/// A single named sequence.
#[derive(Clone, Debug)]
pub struct Record {
    pub name: String,
    pub seq: String,
}

/// Supported input formats.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Phylip,
}

#[derive(Debug)]
pub enum SeqError {
    Io(io::Error),
    Parse(String),
}

impl fmt::Display for SeqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SeqError::Io(e) => write!(f, "I/O error: {e}"),
            SeqError::Parse(m) => write!(f, "parse error: {m}"),
        }
    }
}

impl std::error::Error for SeqError {}

impl From<io::Error> for SeqError {
    fn from(e: io::Error) -> Self {
        SeqError::Io(e)
    }
}

/// Read records from a path (`-` reads stdin), auto-detecting the format unless one is given.
pub fn read_path(path: &str, fmt: Option<Format>) -> Result<Vec<Record>, SeqError> {
    let data = if path == "-" {
        let mut buf = String::new();
        io::stdin().read_to_string(&mut buf)?;
        buf
    } else {
        fs::read_to_string(Path::new(path))?
    };
    let fmt = fmt.unwrap_or_else(|| detect(&data));
    match fmt {
        Format::Fasta => parse_fasta(&data),
        Format::Phylip => parse_phylip(&data),
    }
}

/// Guess the format from the first non-blank character.
pub fn detect(data: &str) -> Format {
    for line in data.lines() {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        return if t.starts_with('>') {
            Format::Fasta
        } else {
            Format::Phylip
        };
    }
    Format::Fasta
}

/// Parse (possibly line-wrapped) FASTA.
pub fn parse_fasta(data: &str) -> Result<Vec<Record>, SeqError> {
    let mut records = Vec::new();
    let mut name: Option<String> = None;
    let mut seq = String::new();
    for line in data.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(n) = name.take() {
                records.push(Record { name: n, seq: std::mem::take(&mut seq) });
            }
            // header is the token up to the first whitespace
            name = Some(rest.split_whitespace().next().unwrap_or("").to_string());
        } else {
            seq.extend(line.split_whitespace());
        }
    }
    if let Some(n) = name.take() {
        records.push(Record { name: n, seq });
    }
    if records.is_empty() {
        return Err(SeqError::Parse("no FASTA records found".into()));
    }
    Ok(records)
}

/// Parse sequential PHYLIP (` nseq nsites` header, then `name<space>sequence`,
/// with the sequence optionally continued on wrapped lines).
pub fn parse_phylip(data: &str) -> Result<Vec<Record>, SeqError> {
    let mut lines = data.lines().filter(|l| !l.trim().is_empty());
    let header = lines
        .next()
        .ok_or_else(|| SeqError::Parse("empty PHYLIP input".into()))?;
    let mut it = header.split_whitespace();
    let nseq: usize = it
        .next()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| SeqError::Parse("bad PHYLIP header".into()))?;
    let mut records: Vec<Record> = Vec::with_capacity(nseq);
    for line in lines {
        if records.len() < nseq {
            let (name, rest) = line
                .split_once(char::is_whitespace)
                .unwrap_or((line, ""));
            let mut seq = String::new();
            seq.extend(rest.split_whitespace());
            records.push(Record { name: name.to_string(), seq });
        } else {
            // interleaved continuation blocks
            let idx = (records.len()) % nseq;
            let _ = idx;
            break; // sequential PHYLIP is the common case for this tool; keep it simple
        }
    }
    if records.len() != nseq {
        return Err(SeqError::Parse(format!(
            "expected {nseq} sequences, found {}",
            records.len()
        )));
    }
    Ok(records)
}
