# Initial Backlog Seed

This backlog is curated from verified code/test/docs evidence plus clearly-marked roadmap items.
Issue count: 20

---

## ISSUE: Feature-gate `fp.rs` and `ext.rs` tests/doctests so `--no-default-features` builds pass
Labels: bug, P0-correctness/blocker, area:no-std, area:field-fp, area:ext-fields
Milestone: v0.1-hardening
Context:
Several test blocks and doctests call APIs that are gated behind `alloc` or `std`, but the
tests themselves lack corresponding `#[cfg(...)]` guards. This breaks `cargo test` and
`cargo test --doc` under `--no-default-features`.
Evidence:
- `src/structures/fp.rs:1376-1560`: unguarded tests calling `discrete_log*`, `multiplicative_order`, `primitive_root` — methods gated at `fp.rs:630` (`std`), `fp.rs:659` (`std`), `fp.rs:730` (`alloc`), `fp.rs:850` (`alloc`).
- `src/structures/ext.rs:1100-1145`: test module uses `Poly` and `vec!` without `alloc` gate.
- `src/structures/ext.rs:25-44`, `ext.rs:699-720`: doctests import `Poly`/`mul_mod` which require `alloc`.
- Reproducer: `cargo test -q --no-default-features` fails; `cargo test -q --doc --no-default-features` fails.
Impact:
- Breaks advertised `no_std` workflow. Blocks CI matrix expansion.
Definition of Done:
- [ ] Add `#[cfg(feature = "alloc")]` / `#[cfg(feature = "std")]` guards to affected test blocks in `fp.rs`.
- [ ] Add `#[cfg(feature = "alloc")]` guard to `ext.rs` test module.
- [ ] Fence alloc-dependent doctests in `ext.rs` with `#[cfg(feature = "alloc")]`.
- [ ] `cargo test -q --no-default-features --lib` passes.
- [ ] `cargo test -q --doc --no-default-features` passes.

## ISSUE: Feature-gate examples to avoid `--no-default-features` compile breakage
Labels: bug, P1-important, area:docs-examples, area:no-std
Milestone: v0.1-hardening
Context:
Examples import alloc-gated modules/types unconditionally, so `cargo test --no-default-features`
tries to compile them and fails.
Evidence:
- `examples/field_basics.rs:114-140` uses APIs gated by `alloc`/`std`.
- `examples/extension_fields.rs:13-15`, `poly_arithmetic.rs:12`, `ntt_mul.rs:10`, `irreducible_polys.rs:11` require alloc-exported types.
- `examples/factorization.rs:11-16` already demonstrates the correct pattern (runtime feature check with message).
Impact:
- Breaks `cargo test --no-default-features`. Hinders contributor onboarding.
Definition of Done:
- [ ] Add compile-time `#[cfg]` gates or `required-features` in `Cargo.toml` for each affected example.
- [ ] `cargo test -q --no-default-features` no longer fails on examples.

## ISSUE: Fix tautological `clippy` assertion in `poly.rs` test
Labels: tech-debt, P1-important, area:polynomials, area:infra-ci
Milestone: v0.1-hardening
Context:
A tautological assertion causes `cargo clippy --all-targets --all-features` to fail.
Evidence:
- `src/structures/poly.rs:2556`: `assert!(is_prim || !is_prim)` — always true.
Impact:
- Blocks clean clippy CI gate.
Definition of Done:
- [ ] Replace with a meaningful assertion or remove.
- [ ] `cargo clippy -q --all-targets --all-features` passes.

## ISSUE: Eliminate unused-import warnings under non-default feature combinations
Labels: tech-debt, P1-important, area:field-fp, area:ext-fields, area:ntt-fft
Milestone: v0.1-hardening
Context:
Feature-conditional imports emit warnings when compiled without their gate features.
Evidence:
- Warnings observed for imports in `src/structures/gf.rs:12`, `ntt.rs:19`, `ext.rs:8` under various feature combinations.
Impact:
- Noisy build output. Undermines warning-free CI.
Definition of Done:
- [ ] Add correct `#[cfg(...)]` to feature-dependent imports.
- [ ] `cargo build -q --no-default-features` and `cargo build -q` both produce zero warnings.

## ISSUE: Add CI workflow matrix for default/all-features/no-default/doc/clippy checks
Labels: tech-debt, P1-important, area:infra-ci
Milestone: v0.1-hardening
Context:
No CI workflows exist in the repository. Multiple regressions (feature-gating, clippy) are
only caught by manual local runs.
Evidence:
- `.github/workflows/` directory is missing.
Definition of Done:
- [ ] Add GitHub Actions workflow with matrix: `cargo test`, `--all-features`, `--no-default-features`, doctest, clippy.
- [ ] Workflow fails on warnings and lint violations.

## ISSUE: Update README capability matrix and feature prerequisites
Labels: docs, P1-important, area:docs-examples
Milestone: v0.1-hardening
Context:
README still describes the crate as a minimal prime-field library despite significant expansion
(extension fields, NTT, factorization, discrete log, tower fields).
Evidence:
- `README.md:3-15` frames crate as minimal.
Impact:
- Misleading onboarding for new users.
Definition of Done:
- [ ] Add module capability table (Fp / Poly / Ext / GF / NTT) with required features.
- [ ] Document which examples require which features.
- [ ] Note `no_std`/`alloc` limitations.

## ISSUE: Add checked-division APIs for `GF` and `TowerGF`
Labels: feature, P1-important, area:ext-fields, area:traits-algebra
Milestone: v0.1-hardening
Context:
Division operators on `GF` and `TowerGF` unconditionally panic on zero divisors (or unsupported
tower degree). Callers have no fallible path.
Evidence:
- `src/structures/fp.rs:965` panics on zero divisor.
- `src/structures/gf.rs:427`, `gf.rs:440` panic on zero.
- `src/structures/gf.rs:843`, `gf.rs:866` panic on zero or `D2 != 2`.
Impact:
- No error-recovery path for caller-controlled inputs.
Definition of Done:
- [ ] Add `try_div` / checked-division methods for `GF` and `TowerGF`.
- [ ] Document `/` operator as panicking.
- [ ] Add tests for both success and error paths.

## ISSUE: Add fallible NTT APIs to replace panic pathways
Labels: feature, P1-important, area:ntt-fft
Milestone: v0.1-hardening
Context:
NTT functions panic when the prime is incompatible or polynomials are too large.
Evidence:
- `src/structures/ntt.rs:527`: `expect("Prime does not support NTT")`.
- `src/structures/ntt.rs:530-535`: asserts on size support.
Impact:
- No error-recovery path for unsupported prime/size combinations.
Definition of Done:
- [ ] Introduce `try_mul_ntt` and fallible plan/transform APIs returning `Result`.
- [ ] Define typed NTT error enum.
- [ ] Keep panic APIs as convenience wrappers.

## ISSUE: Investigate `TowerField::reduce_tower` correctness for `D2 > 2`
Labels: bug, P0-correctness/blocker, area:ext-fields, status:needs-investigation
Milestone: v0.2-math-quality
Context:
The reduction algorithm in `reduce_tower` drops high-degree terms with a comment suggesting
it may need additional passes for `D2 > 2`. No tests exercise this path.
Evidence:
- `src/structures/ext.rs:855-864`: reduction drops terms with `target >= D2`; comment says "another pass would be needed."
- All current tests use `D2 = 2`.
Impact:
- Potential silent arithmetic corruption for tower fields with `D2 > 2`.
Definition of Done:
- [ ] Add tests for `D2 = 3` and `D2 = 4` multiplication/reduction.
- [ ] Confirm algorithm correctness or produce a reproducer showing incorrect results.
- [ ] Fix the algorithm or restrict the API to `D2 = 2` with a compile-time or runtime check.

## ISSUE: Implement generic tower inversion for `D2 > 2` or restrict the API
Labels: feature, P1-important, area:ext-fields
Milestone: v0.2-math-quality
Context:
Tower field inversion only works for `D2 = 2`; other values silently return `None`,
and the division wrappers then panic.
Evidence:
- `src/structures/ext.rs:888-893`: returns `None` for `D2 != 2`.
- `src/structures/gf.rs:843`, `gf.rs:866`: division panics in this case.
Impact:
- Misleading API surface — type permits `D2 > 2` but arithmetic is broken.
Definition of Done:
- [ ] Implement general inversion (extended-GCD over extension polynomials) OR add compile-time/runtime restriction to `D2 = 2`.
- [ ] Ensure non-panicking error path for unsupported cases.
- [ ] Add tests for both supported and unsupported `D2` values.

## ISSUE: Resolve `ExtFieldWithModulus` doc/implementation gap
Labels: docs, P1-important, area:ext-fields
Milestone: v0.2-math-quality
Context:
The `ExtFieldWithModulus` wrapper's documentation claims it enables Ring/Field-style
operations, but only constructor and accessor methods are implemented.
Evidence:
- `src/structures/ext.rs:653-657`: doc claims Ring/Field support.
- `src/structures/ext.rs:668-682`: only `new()`, `elem()`, `modulus()` exist.
Impact:
- Documentation misleads users into expecting trait-based arithmetic.
Definition of Done:
- [ ] Either implement `Ring`/`Field` traits for the wrapper OR update docs to reflect current scope.
- [ ] Add usage examples for the chosen path.

## ISSUE: Define and document constant-time guarantees per API surface
Labels: security, P1-important, area:constant-time, area:field-fp
Milestone: v0.2-math-quality
Context:
The crate exposes `subtle`-gated constant-time helpers but core operations still include
branches and variable-time algorithms. No documentation clarifies which APIs are safe for
secret-dependent inputs.
Evidence:
- `src/structures/fp.rs:280-333`: `pow_ct`, `inverse_ct` (constant-time).
- `src/structures/fp.rs:913-930`: branch-based arithmetic.
- `src/structures/fp.rs:982-1002`: variable-time inversion (extended GCD).
Impact:
- Misuse hazard: users may assume all `Fp` operations are constant-time when `subtle` is enabled.
Definition of Done:
- [ ] Add crate-level and per-function documentation of constant-time vs. variable-time status.
- [ ] Mark variable-time APIs prominently (e.g., `// WARNING: variable-time` in docs).
- [ ] Consider a side-channel review checklist in `SECURITY.md` or crate docs.

## ISSUE: Add strict irreducibility-checked deserialization for `GF` modulus
Labels: security, P1-important, area:serialization, area:ext-fields
Milestone: v0.2-math-quality
Context:
The `to_gf` deserialization path uses `Modulus::new_unchecked`, skipping irreducibility
validation. A reducible modulus silently produces non-field arithmetic.
Evidence:
- `src/structures/gf.rs:1269-1271`: comment notes irreducibility is NOT checked.
- `src/structures/gf.rs:1288`: uses `Modulus::new_unchecked`.
Impact:
- Untrusted input can create a `GF` with broken field invariants.
Definition of Done:
- [ ] Add strict `to_gf_checked` (or similar) that validates irreducibility via `Modulus::new`.
- [ ] Keep existing unchecked path with explicit `_unchecked` naming.
- [ ] Add test: reducible modulus is rejected by strict path.

## ISSUE: Verify and correct binary-field polynomial reference table
Labels: bug, P1-important, area:ext-fields
Milestone: v0.2-math-quality
Context:
The binary-field irreducible polynomial lookup table contains entries marked as approximate
or placeholder, presented alongside correct entries without distinction.
Evidence:
- `src/structures/gf.rs:1108-1114`: comments include "actually ...", "approx", and "Placeholder" for degrees 14, 16, 32, 64, 128.
Impact:
- Incorrect polynomial data for downstream binary-field construction.
Definition of Done:
- [ ] Verify all entries against standard references (e.g., NIST, Lidl & Niederreiter tables).
- [ ] Correct or remove inaccurate entries.
- [ ] Add test for at least AES polynomial (`x^8 + x^4 + x^3 + x + 1`).

## ISSUE: Expand property-based tests to Poly, ExtField, and NTT
Labels: tech-debt, P1-important, area:polynomials, area:ext-fields, area:ntt-fft
Milestone: v0.2-math-quality
Context:
Property-based testing coverage is concentrated on `Fp` in `tests/field_properties.rs`.
Polynomial, extension-field, and NTT modules rely on deterministic unit tests only.
Evidence:
- `tests/field_properties.rs`: strong `Fp` proptests.
- `src/structures/poly.rs`, `ext.rs`, `ntt.rs`: deterministic tests only.
Impact:
- Lower confidence in algebraic invariants outside `Fp`.
Definition of Done:
- [ ] Add proptests for polynomial ring laws (associativity, distributivity, division invariants).
- [ ] Add proptests for extension-field arithmetic identities.
- [ ] Add NTT round-trip and convolution consistency proptests.

## ISSUE: Add criterion benchmarks for factorization and discrete-log
Labels: perf, P1-important, area:benchmarks, area:factorization, area:discrete-log
Milestone: v0.2-math-quality
Context:
Existing benchmarks in `benches/benchmarks.rs` cover basic field and polynomial operations
but lack coverage for factorization and discrete-log, which are the most expensive algorithms.
Evidence:
- `benches/benchmarks.rs:10-181`: no factorization or discrete-log benchmark groups.
Impact:
- No regression detection for the most performance-sensitive code paths.
Definition of Done:
- [ ] Add benchmark groups for `Poly::factor` / `Poly::roots` at representative sizes.
- [ ] Add benchmark groups for `Fp::discrete_log*` with varying subgroup sizes.
- [ ] Verify baselines are stable and suitable for PR comparison.

## ISSUE: Establish NTT vs. naive multiplication crossover benchmarks
Labels: perf, P1-important, area:benchmarks, area:ntt-fft
Milestone: v0.2-math-quality
Context:
The existing NTT benchmark is too narrow to guide callers on when to prefer `mul_ntt` over
naive multiplication.
Evidence:
- `benches/benchmarks.rs:102-114`: limited comparison, no crossover analysis.
Impact:
- Users cannot make informed performance decisions.
Definition of Done:
- [ ] Benchmark across wider degree ranges and multiple primes.
- [ ] Record crossover points.
- [ ] Document a selection heuristic in `ntt` module docs.

## ISSUE: Roadmap: multi-limb field backend for 256-bit primes
Labels: feature, P2-later, area:field-fp, area:traits-algebra
Milestone: v0.3-crypto-roadmap
Context:
The current field core is `Fp<const P: u64>` with `u64`/`u128` arithmetic. Modern
cryptographic applications require 256-bit (or larger) prime fields.
Evidence:
- `src/structures/fp.rs:90-99`: Montgomery form using `u64`.
- `src/structures/fp.rs:940-942`: multiplication via `u128`.
Definition of Done:
- [ ] Define limb-based field element type and trait integration strategy.
- [ ] Provide basic arithmetic and reduction for at least one 256-bit prime.
- [ ] Add conformance tests against known test vectors.

## ISSUE: Roadmap: dedicated `QuadExt<P>` helper for degree-2 extensions
Labels: feature, P2-later, area:ext-fields
Milestone: v0.3-crypto-roadmap
Context:
Degree-2 extensions are the most common in cryptographic constructions, but currently
require the general `ExtField<P, 2>` machinery with explicit modulus passing.
Evidence:
- No specialized `QuadExt` type in the codebase.
Definition of Done:
- [ ] Design `QuadExt` API optimized for common cryptographic patterns.
- [ ] Provide conversions to/from `ExtField<P, 2>`.
- [ ] Add examples and comparative benchmarks.

## ISSUE: Roadmap: scaffold ECC module with curve traits and scalar multiplication
Labels: feature, P2-later, area:traits-algebra
Milestone: v0.3-crypto-roadmap
Context:
No elliptic-curve module exists. This is the primary capability gap for general
cryptography use.
Evidence:
- No ECC module in `src/` layout.
Definition of Done:
- [ ] Define curve point and scalar traits.
- [ ] Add one reference curve implementation skeleton with scalar multiplication.
- [ ] Include constant-time and misuse-hazard notes.
