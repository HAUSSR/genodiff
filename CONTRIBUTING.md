# Contributing to genodiff

Thank you for considering contributing!

## Development setup

```bash
git clone https://github.com/<your-username>/genodiff.git
cd genodiff
cargo build
cargo test           # all 43 tests must pass
```

## Code style

- No external dependencies — the entire project uses only the Rust standard library
- All public functions and structs have doc comments (`///`)
- Every algorithm change must be validated by `cargo test --release`
- Variant-calling accuracy is validated against the built-in simulator:
  run `genodiff sim` + `genodiff call` + `scripts/compare_vcf_norm.py`
  and confirm ≥95% precision/recall before submitting

## Pull request process

1. Fork the repository and create a feature branch
2. Make your changes with tests
3. Run `cargo test --release` and `cargo clippy`
4. Submit a PR with a clear description

## Reporting issues

- Use the issue templates
- For bugs: include the exact command, input file characteristics
  (genome size, number of contigs), and full stderr output
- For feature requests: describe the use case first

## Code of conduct

Be respectful and constructive.
