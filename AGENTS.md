# Repository Guidelines

## Project Structure & Module Organization
- `src/lib.rs` is the library entry point and re-exports core types.
- `src/algebra/` defines abstract traits like `Ring`, `Field`, and `Group`.
- `src/structures/` contains concrete implementations such as `Fp`, `Poly`, and NTT helpers.
- `tests/` holds property-based tests (currently `tests/field_properties.rs`).
- `examples/` includes runnable demos (e.g., `cargo run --example f17`).

## Build, Test, and Development Commands
- `cargo build`: compile the library.
- `cargo test`: run the full test suite.
- `cargo test <test_name>`: run a specific test by name.
- `cargo test --features serde|rand|subtle`: exercise optional feature sets.
- `cargo test --all-features`: validate all optional features together.
- `cargo clippy`: run lints for style and correctness checks.
- `cargo run --example f17`: run the prime field demo.

## Coding Style & Naming Conventions
- Rust 2021 edition with standard `rustfmt` defaults (4-space indent).
- Use Rust idioms: `snake_case` for functions/modules, `UpperCamelCase` for types,
  and `UPPER_SNAKE_CASE` for constants (e.g., `F17::ONE`).
- Prefer explicit, small helper functions for arithmetic operations; keep
  trait implementations close to their types.

## Testing Guidelines
- Tests use `proptest` for property-based coverage (see `tests/field_properties.rs`).
- Name tests in `snake_case` and keep each property focused on one algebraic rule.
- Run targeted checks with `cargo test addition_commutative` or feature-specific
  suites with `cargo test --features subtle`.

## Commit & Pull Request Guidelines
- Commit messages use short, imperative sentences without prefixes
  (e.g., "Add polynomial factorization over finite fields").
- PRs should include a concise summary, tests run, and any feature flags used.
- Add example output or notes when behavior changes (e.g., new arithmetic ops).

## Configuration Notes
- `Fp<P>` expects `P` to be prime; division by zero will panic.
- Optional features are `serde`, `rand`, and `subtle` (constant-time operations).
