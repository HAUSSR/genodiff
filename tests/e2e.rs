//! End-to-end check: random reference -> simulated divergent genome ->
//! in-memory index -> seed mapping. Exercises the same code paths as the CLI.

use genodiff::fasta::FastaRecord;
use genodiff::index::RefIndex;
use genodiff::rng::Rng;
use genodiff::seeds::{collect_anchors, SeedStats};
use genodiff::sim::{self, SimParams};

fn rand_seq(rng: &mut Rng, n: usize) -> Vec<u8> {
    (0..n).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
}

#[test]
fn simulated_genome_maps_back_to_reference() {
    let mut rng = Rng::new(123);
    let refs = vec![
        FastaRecord {
            name: "chr1".into(),
            seq: rand_seq(&mut rng, 1_500_000),
        },
        FastaRecord {
            name: "chr2".into(),
            seq: rand_seq(&mut rng, 500_000),
        },
    ];
    // ~0.1% divergence plus SVs: harder than a typical "near-relative" pair
    let p = SimParams {
        snp: 1500,
        indel: 150,
        del: 10,
        ins: 10,
        inv: 3,
        dup: 3,
        seed: 7,
    };
    let simres = sim::simulate(&refs, &p);

    let idx = RefIndex::build(refs.iter().map(|r| Ok(r.clone())), 15, 10).unwrap();
    let mut total = SeedStats::default();
    for q in &simres.queries {
        let mut st = SeedStats::default();
        let anchors = collect_anchors(&idx, &q.seq, 0, 1000, &mut st);
        assert!(!anchors.is_empty(), "no anchors for {}", q.name);
        total.add(&st);
    }

    let frac = total.matched_fraction();
    assert!(frac >= 0.85, "matched fraction too low: {frac:.3}");
    // sanity on the absolute scale: 2 Mb at w=19 must give a healthy seed count
    assert!(
        total.query_minimizers > 100_000,
        "unexpectedly few minimizers"
    );
}
