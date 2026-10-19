//! Minimal VCF 4.2 writer (sites-only, no sample column) used for simulator
//! truth output; the real variant caller will reuse and extend this.

use std::io::{self, BufWriter, Write};

pub struct VcfWriter<W: Write> {
    out: BufWriter<W>,
}

impl<W: Write> VcfWriter<W> {
    pub fn create(inner: W) -> io::Result<Self> {
        let mut out = BufWriter::new(inner);
        writeln!(out, "##fileformat=VCFv4.2")?;
        writeln!(out, "##source=genodiff-sim")?;
        Ok(VcfWriter { out })
    }

    pub fn add_contig(&mut self, name: &str, length: u64) -> io::Result<()> {
        writeln!(self.out, "##contig=<ID={name},length={length}>")
    }

    pub fn start_body(&mut self) -> io::Result<()> {
        writeln!(self.out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    }

    /// `pos` is 1-based, as written into the file.
    pub fn record(&mut self, chrom: &str, pos: u64, ref_al: &[u8], alt_al: &[u8], info: &str) -> io::Result<()> {
        write!(self.out, "{chrom}\t{pos}\t.\t")?;
        self.out.write_all(ref_al)?;
        self.out.write_all(b"\t")?;
        self.out.write_all(alt_al)?;
        writeln!(self.out, "\t.\t.\t{info}")
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.out.flush()
    }
}
