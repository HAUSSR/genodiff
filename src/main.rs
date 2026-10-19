use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;
use std::sync::mpsc::channel;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use genodiff::call::{call_contig, CallParams, CallStats, Var, VarRec};
use genodiff::chain::{chain_contig, ChainParams};
use genodiff::fasta::{FastaReader, FastaRecord, FastaWriter};
use genodiff::index::RefIndex;
use genodiff::rng::Rng;
use genodiff::seeds::SeedStats;
use genodiff::sim::{self, SimParams};
use genodiff::vcf::VcfWriter;

const USAGE: &str = "genodiff - reference-guided genome comparison (M4: variant + SV calling)

USAGE:
  genodiff index   <ref.fa> <out.gdx> [--k 19] [--w 19]
      Build a minimizer index of the reference genome (reusable for many queries).

  genodiff map     <ref.gdx> <query.fa> [--threads N] [--max-occ 100]
                   [--drift 150] [--max-gap 100000] [--min-anchors 25]
                   [--min-span 500] [--patience 1000] [--keep-secondary]
      Stream query contigs, chain seed anchors into collinear alignment blocks,
      write them to stdout in PAF format (one line per block, not per seed).
      Short chains (< min_anchors or < min_span) are dropped; chains mostly
      contained in a higher-scoring chain are treated as repeat-copy noise and
      also dropped unless --keep-secondary.

  genodiff call    <ref.fa> <query.fa> [-o out.vcf] [--threads N]
                   [--max-occ 100] [--drift 150] [--patience 1000]
                   [--sample NAME] [--sv-seq-cap 100000]
      Full pipeline: index the reference in memory, chain each query contig,
      base-align chain gaps, and write SNP/InDel/SV records to a VCF 4.2 file
      (default: <query stem>.vcf next to the query; -o - writes stdout).
      --sample overrides the sample column name (default: query file stem).
      SV representation: DEL/INV/DUP are reference-derived and written as
      symbolic alleles (<DEL> etc.) with SVTYPE/END/SVLEN; only INS (foreign
      sequence) carries its bases, up to --sv-seq-cap, then symbolic <INS>.

  genodiff sim     <ref.fa> <outdir> [--seed 42] [--snp 500] [--indel 50]
                   [--svdel 5] [--svins 5] [--inv 2] [--dup 2]
      Simulate a mutated genome from the reference and write <outdir>/query.fa
      plus <outdir>/truth.vcf (ground truth for validation).

  genodiff randref <out.fa> <length> [--contigs 1] [--seed 1]
      Generate a random ACGT reference genome (for benchmarks and tests).
";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(args) {
        Ok(code) => code,
        Err(e) => {
            eprintln!("genodiff: error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: Vec<String>) -> Result<ExitCode, String> {
    if args.is_empty() {
        print!("{USAGE}");
        return Ok(ExitCode::SUCCESS);
    }
    match args[0].as_str() {
        "index" => cmd_index(&args[1..]),
        "map" => cmd_map(&args[1..]),
        "call" => cmd_call(&args[1..]),
        "sim" => cmd_sim(&args[1..]),
        "randref" => cmd_randref(&args[1..]),
        "-h" | "--help" | "help" => {
            print!("{USAGE}");
            Ok(ExitCode::SUCCESS)
        }
        other => Err(format!("unknown command '{other}' (try --help)")),
    }
}

// ---------- shared arg helpers ----------

/// flags that are present/absent rather than taking a value
const BOOL_FLAGS: &[&str] = &["--keep-secondary"];

fn positional<'a>(args: &'a [String]) -> Vec<&'a String> {
    // Value flags ("--flag value", "--flag=value", "-o value") consume their
    // argument; boolean flags stand alone. Flag values must not become
    // positionals.
    let mut out = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if a.starts_with("--") {
            if BOOL_FLAGS.contains(&a.as_str()) || a.contains('=') {
                i += 1;
            } else if i + 1 < args.len() {
                i += 2;
            } else {
                i += 1;
            }
        } else if a.starts_with('-') && a.len() > 1 && !a[1..].chars().all(|c| c.is_ascii_digit()) {
            // short value flags like -o
            i += 2;
        } else {
            out.push(a);
            i += 1;
        }
    }
    out
}

fn flag_value(args: &[String], name: &str) -> Option<String> {
    let mut it = args.iter();
    while let Some(a) = it.next() {
        if a == name {
            return it.next().cloned();
        }
        if let Some(v) = a.strip_prefix(&format!("{name}=")) {
            return Some(v.to_string());
        }
    }
    None
}

fn parse_num<T: std::str::FromStr>(name: &str, v: &str) -> Result<T, String> {
    v.parse::<T>()
        .map_err(|_| format!("invalid value for {name}: '{v}'"))
}

fn flag_num<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> Result<T, String> {
    match flag_value(args, name) {
        Some(v) => parse_num(name, &v),
        None => Ok(default),
    }
}

// ---------- genodiff index ----------

fn cmd_index(args: &[String]) -> Result<ExitCode, String> {
    let pos = positional(args);
    if pos.len() != 2 {
        return Err("usage: genodiff index <ref.fa> <out.gdx> [--k 19] [--w 19]".into());
    }
    let k = flag_num::<u32>(args, "--k", 19)?;
    let w = flag_num::<u32>(args, "--w", 19)?;
    if !(1..=32).contains(&k) {
        return Err("--k must be in 1..=32".into());
    }
    if w < 1 {
        return Err("--w must be >= 1".into());
    }
    let t0 = Instant::now();
    let file = File::open(Path::new(pos[0])).map_err(|e| format!("open {}: {e}", pos[0]))?;
    let idx =
        RefIndex::build(FastaReader::new(file), k, w).map_err(|e| format!("build index: {e}"))?;
    let total_bp: u64 = idx.contigs.iter().map(|c| c.length).sum();
    idx.save(Path::new(pos[1]))
        .map_err(|e| format!("save {}: {e}", pos[1]))?;
    eprintln!(
        "[genodiff] index: {} contigs, {} bp, {} minimizers (k={k}, w={w}) -> {} [{:.2}s]",
        idx.contigs.len(),
        total_bp,
        idx.hits.len(),
        pos[1],
        t0.elapsed().as_secs_f64()
    );
    Ok(ExitCode::SUCCESS)
}

// ---------- genodiff map (contig-parallel, streamed, chained) ----------

struct Job {
    seq_idx: usize,
    name: Arc<str>,
    qlen: u64,
    seq: Arc<Vec<u8>>,
}

struct JobOut {
    seq_idx: usize,
    lines: String,
    stats: SeedStats,
    chains: u64,
    anchors_in_chains: u64,
    short_dropped: u64,
    secondary_dropped: u64,
}

fn process_job(idx: &RefIndex, job: &Job, max_occ: usize, params: ChainParams) -> JobOut {
    let (chains, cstats, stats) = chain_contig(idx, &job.seq, max_occ, params, false);
    let mut lines = String::with_capacity(chains.len() * 112 + 64);
    for c in &chains {
        let rc = &idx.contigs[c.ref_contig as usize];
        let strand = if c.strand == 0 { '+' } else { '-' };
        let _ = std::fmt::Write::write_fmt(
            &mut lines,
            format_args!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tna:{}\n",
                job.name,
                job.qlen,
                c.q_start,
                c.q_end,
                strand,
                rc.name,
                rc.length,
                c.r_start,
                c.r_end,
                c.n_anchors,
                c.n_anchors
            ),
        );
    }
    JobOut {
        seq_idx: job.seq_idx,
        lines,
        stats,
        chains: cstats.chains,
        anchors_in_chains: cstats.anchors_in_chains,
        short_dropped: cstats.short_dropped,
        secondary_dropped: cstats.secondary_dropped,
    }
}

fn cmd_map(args: &[String]) -> Result<ExitCode, String> {
    let pos = positional(args);
    if pos.len() != 2 {
        return Err(
            "usage: genodiff map <ref.gdx> <query.fa> [--threads N] [--max-occ 100] \
             [--drift 150] [--max-gap 100000] [--min-anchors 25] [--min-span 500] \
             [--patience 1000] [--keep-secondary]"
                .into(),
        );
    }
    let max_occ: usize = flag_num::<usize>(args, "--max-occ", 100)?;
    let params = ChainParams {
        max_drift: flag_num::<i64>(args, "--drift", 150)?,
        tol_cap: 2000,
        max_chain_gap: flag_num::<u64>(args, "--max-gap", 100_000)?,
        min_anchors: flag_num::<u32>(args, "--min-anchors", 25)?,
        min_span: flag_num::<u64>(args, "--min-span", 500)?,
        suppress_secondary: !args.iter().any(|a| a == "--keep-secondary"),
        patience: flag_num::<u64>(args, "--patience", 1_000)?,
    };
    let threads: usize = flag_num::<usize>(
        args,
        "--threads",
        thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1),
    )?
    .max(1);
    // cap how many contigs are held in flight (each fully in memory while it
    // is being chained); 4 GiB of sequence in flight is plenty of work for
    // any thread count while bounding peak RSS
    const MAX_INFLIGHT_BYTES: u64 = 4 << 30;
    let max_inflight_contigs = threads.max(1).next_power_of_two();

    let t0 = Instant::now();
    let idx = Arc::new(
        RefIndex::load(Path::new(pos[0])).map_err(|e| format!("load index {}: {e}", pos[0]))?,
    );
    let qfile = File::open(Path::new(pos[1])).map_err(|e| format!("open {}: {e}", pos[1]))?;
    let mut reader = FastaReader::new(qfile);

    let (job_tx, job_rx) = channel::<Job>();
    let job_rx = Arc::new(Mutex::new(job_rx));
    let (res_tx, res_rx) = channel::<JobOut>();
    let mut handles = Vec::with_capacity(threads);
    for _ in 0..threads {
        let rx = Arc::clone(&job_rx);
        let tx = res_tx.clone();
        let idx = Arc::clone(&idx);
        handles.push(thread::spawn(move || loop {
            let job = {
                let guard = rx.lock().unwrap();
                guard.recv()
            };
            let job = match job {
                Ok(j) => j,
                Err(_) => break,
            };
            if tx.send(process_job(&idx, &job, max_occ, params)).is_err() {
                break;
            }
        }));
    }
    drop(res_tx);

    let stdout = io::stdout();
    let mut out = BufWriter::with_capacity(1 << 20, stdout.lock());
    let mut done_reading = false;
    let mut in_flight = 0usize;
    let mut inflight_bytes: u64 = 0;
    let mut submitted = 0usize;
    let mut next_write = 0usize;
    let mut buffer: std::collections::HashMap<usize, JobOut> = std::collections::HashMap::new();
    let mut total = SeedStats::default();
    let (mut total_chains, mut total_in_chains) = (0u64, 0u64);
    let (mut total_short, mut total_sec) = (0u64, 0u64);

    loop {
        // submit contigs while capacity allows
        while !done_reading
            && in_flight < max_inflight_contigs
            && inflight_bytes < MAX_INFLIGHT_BYTES
        {
            match reader.next_record() {
                Ok(Some(r)) => {
                    let job = Job {
                        seq_idx: submitted,
                        name: Arc::from(r.name.as_str()),
                        qlen: r.seq.len() as u64,
                        seq: Arc::new(r.seq),
                    };
                    inflight_bytes += job.qlen;
                    if job_tx.send(job).is_err() {
                        return Err("worker thread panicked".into());
                    }
                    submitted += 1;
                    in_flight += 1;
                }
                Ok(None) => {
                    done_reading = true;
                    break;
                }
                Err(e) => return Err(format!("read query: {e}")),
            }
        }
        if in_flight == 0 {
            break;
        }
        let res = res_rx.recv().map_err(|_| "worker threads died")?;
        in_flight -= 1;
        buffer.insert(res.seq_idx, res);
        while let Some(r) = buffer.remove(&next_write) {
            out.write_all(r.lines.as_bytes())
                .map_err(|e| format!("write output: {e}"))?;
            total.add(&r.stats);
            total_chains += r.chains;
            total_in_chains += r.anchors_in_chains;
            total_short += r.short_dropped;
            total_sec += r.secondary_dropped;
            next_write += 1;
        }
    }
    drop(job_tx);
    for h in handles {
        let _ = h.join();
    }
    out.flush().map_err(|e| format!("write output: {e}"))?;

    eprintln!(
        "[genodiff] map: {} contig(s) on {} thread(s): minimizers={}, matched={} ({:.1}%), \
         repeat-skipped={}, chains={} (anchors in chains: {}, short dropped: {}, \
         secondary dropped: {}) [{:.2}s]",
        submitted,
        threads,
        total.query_minimizers,
        total.matched_minimizers,
        100.0 * total.matched_fraction(),
        0,
        total_chains,
        total_in_chains,
        total_short,
        total_sec,
        t0.elapsed().as_secs_f64()
    );
    Ok(ExitCode::SUCCESS)
}

// ---------- genodiff call (variant calling -> VCF) ----------

struct CallJob {
    seq_idx: usize,
    name: Arc<str>,
    seq: Arc<Vec<u8>>,
}

struct CallOut {
    seq_idx: usize,
    name: Arc<str>,
    vars: Vec<VarRec>,
    stats: CallStats,
    #[allow(dead_code)]
    bases: u64,
}

fn write_vcf_header<W: Write>(
    out: &mut W,
    idx: &RefIndex,
    ref_seqs: &[Vec<u8>],
    ref_path: &str,
    query_path: &str,
    sample: &str,
) -> io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(out, "##source=genodiff")?;
    // the genome every CHROM/POS/END refers to, and the genome the sample
    // was called from — both needed to extract SV sequences later
    writeln!(out, "##reference=file:{ref_path}")?;
    writeln!(out, "##genodiff_query=file:{query_path}")?;
    for (c, seq) in idx.contigs.iter().zip(ref_seqs) {
        writeln!(out, "##contig=<ID={},length={}>", c.name, seq.len())?;
    }
    writeln!(out, "##ALT=<ID=DEL,Description=\"Deletion\">")?;
    writeln!(out, "##ALT=<ID=INS,Description=\"Insertion\">")?;
    writeln!(out, "##ALT=<ID=INV,Description=\"Inversion\">")?;
    writeln!(out, "##ALT=<ID=DUP,Description=\"Duplication\">")?;
    writeln!(out, "##ALT=<ID=BND,Description=\"Breakend\">")?;
    writeln!(
        out,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    )?;
    writeln!(out, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant on the reference\">")?;
    writeln!(
        out,
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the structural variant\">"
    )?;
    writeln!(out, "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome of second breakend (reference)\">")?;
    writeln!(out, "##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Position of second breakend (reference)\">")?;
    writeln!(out, "##INFO=<ID=QSEQ,Number=1,Type=String,Description=\"Query contig of this variant's query-side span (see ##genodiff_query)\">")?;
    writeln!(out, "##INFO=<ID=QSTART,Number=1,Type=Integer,Description=\"1-based start of the query-side span; extract the sequence with e.g. samtools faidx query.fa QSEQ:QSTART-QEND\">")?;
    writeln!(
        out,
        "##INFO=<ID=QEND,Number=1,Type=Integer,Description=\"1-based end of the query-side span\">"
    )?;
    writeln!(
        out,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"
    )
}

/// append QSEQ/QSTART/QEND to the INFO column when the variant carries a
/// query-side span (the exact place the involved sequence lives in the
/// query assembly)
fn write_qspan<W: Write>(out: &mut W, qname: &str, qspan: Option<(u64, u64)>) -> io::Result<()> {
    match qspan {
        Some((s, e)) => write!(out, ";QSEQ={qname};QSTART={};QEND={e}", s + 1),
        None => Ok(()),
    }
}

fn write_vcf_record<W: Write>(
    out: &mut W,
    idx: &RefIndex,
    ref_seqs: &[Vec<u8>],
    v: &VarRec,
    sv_seq_cap: usize,
    qname: &str,
) -> io::Result<()> {
    let rseq = &ref_seqs[v.contig as usize];
    let chrom = idx.contig_name(v.contig);
    match &v.var {
        Var::Snp { ref_base, alt_base } => {
            write!(out, "{chrom}\t{}\t.\t", v.pos)?;
            out.write_all(&[*ref_base])?;
            out.write_all(b"\t")?;
            out.write_all(&[*alt_base])?;
            writeln!(out, "\t.\t.\t.\tGT\t1")
        }
        Var::Del { len } => {
            let p = v.pos as usize; // anchor base is rseq[p-1]
            write!(out, "{chrom}\t{p}\t.\t")?;
            out.write_all(&rseq[p - 1..p + *len as usize])?;
            write!(out, "\t")?;
            out.write_all(&rseq[p - 1..p])?;
            writeln!(out, "\t.\t.\t.\tGT\t1")
        }
        Var::Ins { seq } => {
            let p = v.pos as usize;
            write!(out, "{chrom}\t{p}\t.\t")?;
            out.write_all(&rseq[p - 1..p])?;
            write!(out, "\t")?;
            out.write_all(&rseq[p - 1..p])?;
            out.write_all(seq)?;
            writeln!(out, "\t.\t.\t.\tGT\t1")
        }
        Var::SvDel { seq } => {
            // the deleted bases are known reference sequence: describe the
            // span instead of copying it (v.pos is the anchor base, 1-based;
            // the deleted span is the seq.len() bases after it)
            let p = (v.pos.max(1)) as usize;
            writeln!(
                out,
                "{chrom}\t{p}\t.\t{}\t<DEL>\t.\t.\tSVTYPE=DEL;END={};SVLEN={}\tGT\t1",
                rseq[p - 1] as char,
                p + seq.len(),
                seq.len()
            )
        }
        Var::SvIns { seq } => {
            // foreign sequence absent from the reference: the only SV class
            // whose bases must be written out (REF/ALT = anchor + inserted);
            // QSEQ/QSTART/QEND say where the same bases sit in the query
            let p = (v.pos.max(1)) as usize;
            write!(out, "{chrom}\t{p}\t.\t")?;
            out.write_all(&rseq[p - 1..p])?;
            write!(out, "\t")?;
            if seq.len() <= sv_seq_cap {
                out.write_all(&rseq[p - 1..p])?;
                out.write_all(seq)?;
            } else {
                write!(out, "<INS>")?;
            }
            write!(out, "\t.\t.\tSVTYPE=INS;SVLEN={}", seq.len())?;
            write_qspan(out, qname, v.qspan)?;
            writeln!(out, "\tGT\t1")
        }
        Var::SvInv { seq } => {
            // the inverted bases are reference sequence (reversed): describe
            // the span; v.pos is the 0-based first inverted base, so VCF POS
            // is s+1 and END the 1-based last inverted base
            let s = (v.pos.max(1)) as usize;
            write!(
                out,
                "{chrom}\t{}\t.\t{}\t<INV>\t.\t.\tSVTYPE=INV;END={};SVLEN={}",
                s + 1,
                rseq[s] as char,
                s + seq.len(),
                seq.len()
            )?;
            write_qspan(out, qname, v.qspan)?;
            writeln!(out, "\tGT\t1")
        }
        Var::SvDup { seq } => {
            // duplicated span is reference sequence: describe it (anchor
            // form: the copy is inserted after the anchor base; the span is
            // the seq.len() ref bases from 0-based v.pos, 1-based END
            // = v.pos + seq.len())
            let p = (v.pos.max(1)) as usize;
            write!(
                out,
                "{chrom}\t{p}\t.\t{}\t<DUP>\t.\t.\tSVTYPE=DUP;END={};SVLEN={}",
                rseq[p - 1] as char,
                p + seq.len(),
                seq.len()
            )?;
            write_qspan(out, qname, v.qspan)?;
            writeln!(out, "\tGT\t1")
        }
        Var::SvTra { chr2, pos2 } => {
            let p = (v.pos.max(1)) as usize;
            writeln!(
                out,
                "{chrom}\t{p}\t.\t{}\t<BND>\t.\t.\tSVTYPE=BND;CHR2={};POS2={pos2}\tGT\t1",
                rseq[p - 1] as char,
                idx.contig_name(*chr2)
            )
        }
    }
}

fn cmd_call(args: &[String]) -> Result<ExitCode, String> {
    let pos = positional(args);
    if pos.len() != 2 {
        return Err(
            "usage: genodiff call <ref.fa> <query.fa> [-o out.vcf] [--threads N] \
             [--max-occ 100] [--drift 150] [--patience 1000] [--sample NAME]"
                .into(),
        );
    }
    let (ref_path, query_path) = (pos[0], pos[1]);
    let out_path = flag_value(args, "-o");
    let sample_override = flag_value(args, "--sample");
    let sv_seq_cap: usize = flag_num::<usize>(args, "--sv-seq-cap", 1000)?;
    let threads: usize = flag_num::<usize>(
        args,
        "--threads",
        thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1),
    )?
    .max(1);
    let params = CallParams {
        chain: ChainParams {
            max_drift: flag_num::<i64>(args, "--drift", 150)?,
            max_chain_gap: flag_num::<u64>(args, "--max-gap", 100_000)?,
            min_anchors: flag_num::<u32>(args, "--min-anchors", 25)?,
            min_span: flag_num::<u64>(args, "--min-span", 500)?,
            suppress_secondary: !args.iter().any(|a| a == "--keep-secondary"),
            patience: flag_num::<u64>(args, "--patience", 1_000)?,
            ..ChainParams::default()
        },
        max_occ: flag_num::<usize>(args, "--max-occ", 100)?,
        ..CallParams::default()
    };

    let t0 = Instant::now();
    // 1) reference: sequences + in-memory index
    let ref_file = File::open(Path::new(ref_path)).map_err(|e| format!("open {ref_path}: {e}"))?;
    let refs: Vec<FastaRecord> = FastaReader::new(ref_file)
        .collect::<io::Result<Vec<_>>>()
        .map_err(|e| format!("read reference: {e}"))?;
    if refs.is_empty() {
        return Err("reference contains no sequences".into());
    }
    let ref_seqs: Vec<Vec<u8>> = refs.iter().map(|r| r.seq.clone()).collect();
    eprintln!(
        "[genodiff] call: indexing reference ({} contigs, {} bp) ...",
        refs.len(),
        ref_seqs.iter().map(|s| s.len()).sum::<usize>()
    );
    let idx = Arc::new(
        RefIndex::build(refs.iter().map(|r| Ok(r.clone())), 19, 19)
            .map_err(|e| format!("build index: {e}"))?,
    );

    // 2) parallel per-contig calling, results buffered in query order
    let (job_tx, job_rx) = channel::<CallJob>();
    let job_rx = Arc::new(Mutex::new(job_rx));
    let (res_tx, res_rx) = channel::<CallOut>();
    let ref_seqs_arc = Arc::new(ref_seqs.clone());
    let mut handles = Vec::with_capacity(threads);
    for _ in 0..threads {
        let rx = Arc::clone(&job_rx);
        let tx = res_tx.clone();
        let idx = Arc::clone(&idx);
        let ref_seqs = Arc::clone(&ref_seqs_arc);
        let params = params.clone();
        handles.push(thread::spawn(move || loop {
            let job = {
                let guard = rx.lock().unwrap();
                guard.recv()
            };
            let job = match job {
                Ok(j) => j,
                Err(_) => break,
            };
            let (vars, stats) = call_contig(&idx, &ref_seqs, &job.seq, &params);
            if tx
                .send(CallOut {
                    seq_idx: job.seq_idx,
                    name: Arc::clone(&job.name),
                    vars,
                    stats,
                    bases: job.seq.len() as u64,
                })
                .is_err()
            {
                break;
            }
        }));
    }
    drop(res_tx);

    // 3) stream query contigs in, collect outputs in order
    let qfile = File::open(Path::new(query_path)).map_err(|e| format!("open {query_path}: {e}"))?;
    let mut reader = FastaReader::new(qfile);
    let sample = sample_override.unwrap_or_else(|| {
        Path::new(query_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("query")
            .to_string()
    });
    let default_out = format!("{}.vcf", Path::new(query_path).with_extension("").display());
    let vcf_path = out_path.unwrap_or(default_out);

    let mut submitted = 0usize;
    let mut next_write = 0usize;
    let mut buffer: std::collections::HashMap<usize, CallOut> = std::collections::HashMap::new();
    let mut total = CallStats::default();
    let mut done_reading = false;
    let mut in_flight = 0usize;

    let stdout = io::stdout();
    let mut vcf: Box<dyn Write> = if vcf_path == "-" {
        Box::new(BufWriter::new(stdout.lock()))
    } else {
        Box::new(BufWriter::with_capacity(
            1 << 20,
            File::create(Path::new(&vcf_path)).map_err(|e| format!("create {vcf_path}: {e}"))?,
        ))
    };
    write_vcf_header(&mut vcf, &idx, &ref_seqs, ref_path, query_path, &sample)
        .map_err(|e| format!("write VCF: {e}"))?;

    loop {
        while !done_reading && in_flight < threads.max(2) * 2 {
            match reader.next_record() {
                Ok(Some(r)) => {
                    let job = CallJob {
                        seq_idx: submitted,
                        name: Arc::from(r.name.as_str()),
                        seq: Arc::new(r.seq),
                    };
                    if job_tx.send(job).is_err() {
                        return Err("worker thread panicked".into());
                    }
                    submitted += 1;
                    in_flight += 1;
                }
                Ok(None) => {
                    done_reading = true;
                    break;
                }
                Err(e) => return Err(format!("read query: {e}")),
            }
        }
        if in_flight == 0 {
            break;
        }
        let res = res_rx.recv().map_err(|_| "worker threads died")?;
        in_flight -= 1;
        buffer.insert(res.seq_idx, res);
        while let Some(r) = buffer.remove(&next_write) {
            // variants within a contig come out unordered: sort by position
            let mut vars = r.vars;
            vars.sort_by_key(|v| (v.contig, v.pos));
            for v in &vars {
                write_vcf_record(&mut vcf, &idx, &ref_seqs, v, sv_seq_cap, &r.name)
                    .map_err(|e| format!("write VCF: {e}"))?;
            }
            total.add(&r.stats);
            next_write += 1;
        }
    }
    drop(job_tx);
    for h in handles {
        let _ = h.join();
    }
    vcf.flush().map_err(|e| format!("write VCF: {e}"))?;

    eprintln!(
        "[genodiff] call: {} contig(s) on {} thread(s): SNP={}, InDel={}+{}, \
         SV: DEL={} INS={} INV={} DUP={} BND={} \
         [chains used {} (junk {}, repeat-flagged {}), gaps: aligned={} pure={} slip={} sv={} skipped={} weird={}] \
         -> {} [{:.2}s]",
        submitted,
        threads,
        total.snps,
        total.dels,
        total.ins,
        total.sv_dels,
        total.sv_ins,
        total.sv_inv,
        total.sv_dup,
        total.sv_tra,
        total.chains_used,
        total.chains_junk,
        total.chains_repeat,
        total.gaps_aligned,
        total.gaps_pure,
        total.gaps_slip,
        total.gaps_sv,
        total.gaps_skipped,
        total.gaps_weird,
        if vcf_path == "-" { "(stdout)".to_string() } else { vcf_path.clone() },
        t0.elapsed().as_secs_f64()
    );
    Ok(ExitCode::SUCCESS)
}

// ---------- genodiff sim ----------

fn cmd_sim(args: &[String]) -> Result<ExitCode, String> {
    let pos = positional(args);
    if pos.len() != 2 {
        return Err(
            "usage: genodiff sim <ref.fa> <outdir> [--seed 42] [--snp 500] [--indel 50] \
             [--svdel 5] [--svins 5] [--inv 2] [--dup 2]"
                .into(),
        );
    }
    let params = SimParams {
        snp: flag_num::<u64>(args, "--snp", 500)?,
        indel: flag_num::<u64>(args, "--indel", 50)?,
        del: flag_num::<u64>(args, "--svdel", 5)?,
        ins: flag_num::<u64>(args, "--svins", 5)?,
        inv: flag_num::<u64>(args, "--inv", 2)?,
        dup: flag_num::<u64>(args, "--dup", 2)?,
        seed: flag_num::<u64>(args, "--seed", 42)?,
    };
    let file = File::open(Path::new(pos[0])).map_err(|e| format!("open {}: {e}", pos[0]))?;
    let refs: Vec<FastaRecord> = FastaReader::new(file)
        .collect::<io::Result<Vec<_>>>()
        .map_err(|e| format!("read reference: {e}"))?;
    if refs.is_empty() {
        return Err("reference contains no sequences".into());
    }

    let t0 = Instant::now();
    let res = sim::simulate(&refs, &params);

    let outdir = Path::new(pos[1]);
    fs::create_dir_all(outdir).map_err(|e| format!("create {}: {e}", pos[1]))?;

    {
        let mut fw = FastaWriter::new(
            File::create(outdir.join("query.fa")).map_err(|e| format!("write query.fa: {e}"))?,
        );
        for q in &res.queries {
            fw.write_record(&q.name, &q.seq)
                .map_err(|e| format!("write query.fa: {e}"))?;
        }
        fw.flush().map_err(|e| format!("write query.fa: {e}"))?;
    }
    {
        let vf =
            File::create(outdir.join("truth.vcf")).map_err(|e| format!("write truth.vcf: {e}"))?;
        let mut vw = VcfWriter::create(vf).map_err(|e| format!("write truth.vcf: {e}"))?;
        sim::write_truth_vcf(&mut vw, &refs, &res.variants)
            .map_err(|e| format!("write truth.vcf: {e}"))?;
    }

    use genodiff::sim::VarKind;
    let mut counts = std::collections::BTreeMap::new();
    for (_, v) in &res.variants {
        let label = match v {
            VarKind::Snp { .. } => "SNP",
            VarKind::Del { len, .. } if *len < 50 => "small-DEL",
            VarKind::Del { .. } => "DEL",
            VarKind::Ins { seq, .. } if seq.len() < 50 => "small-INS",
            VarKind::Ins { .. } => "INS",
            VarKind::Inv { .. } => "INV",
            VarKind::Dup { .. } => "DUP",
        };
        *counts.entry(label).or_insert(0usize) += 1;
    }
    let summary: Vec<String> = counts.iter().map(|(k, n)| format!("{k}={n}")).collect();
    eprintln!(
        "[genodiff] simulated {} variant(s) [{}] -> {} [{:.2}s]",
        res.variants.len(),
        summary.join(", "),
        outdir.display(),
        t0.elapsed().as_secs_f64()
    );
    Ok(ExitCode::SUCCESS)
}

// ---------- genodiff randref ----------

fn cmd_randref(args: &[String]) -> Result<ExitCode, String> {
    let pos = positional(args);
    if pos.len() != 2 {
        return Err("usage: genodiff randref <out.fa> <length> [--contigs 1] [--seed 1]".into());
    }
    let total: u64 = parse_num("<length>", pos[1])?;
    let ncontigs: u64 = flag_num::<u64>(args, "--contigs", 1)?.max(1);
    let seed: u64 = flag_num::<u64>(args, "--seed", 1)?;
    let mut rng = Rng::new(seed);

    let mut fw = FastaWriter::new(
        File::create(Path::new(pos[0])).map_err(|e| format!("create {}: {e}", pos[0]))?,
    );
    let base = total / ncontigs;
    let rem = total % ncontigs;
    let mut written = 0u64;
    for i in 0..ncontigs {
        let len = base + u64::from(i < rem);
        let seq: Vec<u8> = (0..len).map(|_| rng.base()).collect();
        fw.write_record(&format!("ctg{}", i + 1), &seq)
            .map_err(|e| format!("write: {e}"))?;
        written += len;
    }
    fw.flush().map_err(|e| format!("write: {e}"))?;
    eprintln!(
        "[genodiff] random reference: {} bp in {} contig(s) -> {}",
        written, ncontigs, pos[0]
    );
    Ok(ExitCode::SUCCESS)
}
