//! Streaming FASTA reader/writer. Contigs are read one at a time so the
//! whole query genome is never held in memory.

use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};

#[derive(Clone)]
pub struct FastaRecord {
    pub name: String,
    pub seq: Vec<u8>,
}

pub struct FastaReader<R: Read> {
    reader: BufReader<R>,
    pending_header: Option<String>,
}

impl<R: Read> FastaReader<R> {
    pub fn new(inner: R) -> Self {
        FastaReader {
            reader: BufReader::with_capacity(1 << 20, inner),
            pending_header: None,
        }
    }

    pub fn next_record(&mut self) -> io::Result<Option<FastaRecord>> {
        let header = match self.pending_header.take() {
            Some(h) => h,
            None => loop {
                let mut line = String::new();
                if self.reader.read_line(&mut line)? == 0 {
                    return Ok(None);
                }
                let t = line.trim_end();
                if t.starts_with('>') {
                    break t[1..].to_string();
                }
                // tolerate junk/blank lines before the first header
            },
        };
        let name = header.split_whitespace().next().unwrap_or("").to_string();
        let mut seq: Vec<u8> = Vec::with_capacity(1 << 16);
        loop {
            let mut line = String::new();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 {
                break;
            }
            let t = line.trim_end();
            if t.starts_with('>') {
                self.pending_header = Some(t[1..].to_string());
                break;
            }
            if t.is_empty() {
                continue;
            }
            seq.extend_from_slice(t.as_bytes());
        }
        seq.make_ascii_uppercase();
        Ok(Some(FastaRecord { name, seq }))
    }
}

impl<R: Read> Iterator for FastaReader<R> {
    type Item = io::Result<FastaRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_record().transpose()
    }
}

pub struct FastaWriter<W: Write> {
    out: BufWriter<W>,
}

impl<W: Write> FastaWriter<W> {
    pub fn new(inner: W) -> Self {
        FastaWriter {
            out: BufWriter::new(inner),
        }
    }

    pub fn write_record(&mut self, name: &str, seq: &[u8]) -> io::Result<()> {
        writeln!(self.out, ">{name}")?;
        for chunk in seq.chunks(60) {
            self.out.write_all(chunk)?;
            self.out.write_all(b"\n")?;
        }
        Ok(())
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.out.flush()
    }
}
