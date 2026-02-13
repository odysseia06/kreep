# Initial Backlog Seed

This backlog is curated from verified code/test/docs evidence plus clearly-marked roadmap items.

## ISSUE: Gate `fp.rs` tests that require `alloc/std` so `--no-default-features` test build succeeds
Labels: bug, P0-correctness/blocker, area:no-std, area:field-fp, area:infra-ci
Milestone: v0.1-hardening
Evidence:
- `src/structures/fp.rs:1376-1560` has unguarded tests calling `discrete_log*`, `multiplicative_order`, and `primitive_root`.
- These methods are feature-gated (`std`/`alloc`) in `src/structures/fp.rs:630`, `src/structures/fp.rs:659`, `src/structures/fp.rs:730`, `src/structures/fp.rs:850`.
- `cargo test -q --no-default-features` currently fails.
Impact:
- Portability and correctness of advertised no_std workflow.
Definition of Done:
- [ ] Add appropriate `#[cfg(...)]` guards for feature-dependent test blocks in `fp.rs`.
- [ ] `cargo test -q --no-default-features --lib` passes.
- [ ] `cargo test -q --no-default-features` no longer fails due to `fp.rs` feature mismatches.

## ISSUE: Gate `ext.rs` tests/doctests that depend on `alloc` APIs
Labels: bug, P0-correctness/blocker, area:no-std, area:ext-fields, area:docs-examples
Milestone: v0.1-hardening
Evidence:
- `src/structures/ext.rs:1100-1145` test module uses `Poly` and `vec!` while no `alloc` gating is present.
- `src/structures/ext.rs:25-44` and `src/structures/ext.rs:699-720` doctest examples import/use `Poly` and `mul_mod` APIs gated by `alloc`.
- `cargo test -q --doc --no-default-features` fails on those doctests.
Impact:
- Portability and docs correctness for no_std users.
Definition of Done:
- [ ] Add `cfg` guards for alloc-dependent test/doctest content.
- [ ] Use `#[cfg(feature = "alloc")]` doctest fences where needed.
- [ ] `cargo test -q --doc --no-default-features` passes.

## ISSUE: Feature-gate examples to avoid no-default-features compile breakage
Labels: bug, P1-important, area:docs-examples, area:no-std, area:infra-ci
Milestone: v0.1-hardening
Evidence:
- `examples/field_basics.rs:114-140` uses APIs gated by `alloc/std`.
- `examples/extension_fields.rs:13-15`, `examples/poly_arithmetic.rs:12`, `examples/ntt_mul.rs:10`, `examples/irreducible_polys.rs:11` require alloc-exported modules/types.
- `cargo test -q --no-default-features` currently tries to compile incompatible examples and fails.
Impact:
- Portability and contributor workflow reliability.
Definition of Done:
- [ ] Add compile-time `cfg` gates or split examples by feature set.
- [ ] For gated-out examples, print clear run instructions (similar to `examples/factorization.rs:11-16`).
- [ ] `cargo test -q --no-default-features` no longer fails on examples.

## ISSUE: Add CI matrix for default/all-features/no-default/doc/clippy checks
Labels: tech-debt, P1-important, area:infra-ci, area:no-std, area:docs-examples
Milestone: v0.1-hardening
Evidence:
- No CI workflows currently in repository (`.github/workflows/*` missing).
- Multiple regressions are only visible when running non-default matrices locally.
Impact:
- Correctness and release confidence.
Definition of Done:
- [ ] Add a GitHub Actions workflow covering `cargo test`, `--all-features`, and no-default builds/tests.
- [ ] Add doctest and clippy jobs.
- [ ] Fail fast on warnings/lints that currently regress quality gates.

## ISSUE: Fix `clippy --all-targets --all-features` failure in `poly.rs` tests
Labels: tech-debt, P1-important, area:polynomials, area:infra-ci
Milestone: v0.1-hardening
Evidence:
- `cargo clippy -q --all-targets --all-features` fails.
- Tautological assertion at `src/structures/poly.rs:2556` (`assert!(is_prim || !is_prim)`).
Impact:
- Tooling quality gate reliability.
Definition of Done:
- [ ] Replace tautological assertion with meaningful deterministic assertion.
- [ ] Re-run clippy until clean for this error path.
- [ ] Add regression test intent comment or stronger assertion logic.

## ISSUE: Provide checked-division APIs alongside panic-based `/` operators
Labels: bug, P1-important, area:traits-algebra, area:field-fp, area:ext-fields
Milestone: v0.1-hardening
Evidence:
- `src/structures/fp.rs:965` panics on division by zero.
- `src/structures/gf.rs:427`, `src/structures/gf.rs:440` panic on zero divisors.
- `src/structures/gf.rs:843`, `src/structures/gf.rs:866` panic on zero divisors and unsupported tower degree.
Impact:
- Correctness and API ergonomics in caller-controlled error paths.
Definition of Done:
- [ ] Add `try_div`/checked division style APIs for `GF` and `TowerGF`.
- [ ] Update docs/examples to use checked APIs where appropriate.
- [ ] Keep `/` semantics explicit and documented as panic-prone.

## ISSUE: Add fallible NTT APIs to replace panic pathways
Labels: feature, P1-important, area:ntt-fft, area:polynomials
Milestone: v0.1-hardening
Evidence:
- `src/structures/ntt.rs:527` uses `expect("Prime does not support NTT")`.
- `src/structures/ntt.rs:530-535` asserts on size support.
- `NttPlan` methods also assert on length/size assumptions.
Impact:
- Correctness and API ergonomics for production callers.
Definition of Done:
- [ ] Introduce `Result`-returning variants (e.g., `try_mul_ntt`, fallible plan/transform APIs).
- [ ] Define and expose typed NTT errors.
- [ ] Keep panic APIs as thin wrappers or deprecate with migration guidance.

## ISSUE: Investigate `TowerField::reduce_tower` correctness for `D2 > 2`
Labels: bug, P0-correctness/blocker, area:ext-fields, area:polynomials, status:needs-investigation
Milestone: v0.2-math-quality
Evidence:
- In `src/structures/ext.rs:855-864`, reduction drops terms with `target >= D2` and comments that another pass would be needed.
- Current tests focus on `Tower2x2` (`D2 = 2`) and do not cover higher `D2` paths.
Impact:
- Potential arithmetic correctness bug for generic tower configurations.
Definition of Done:
- [ ] Add focused tests for `D2 > 2` multiplication/reduction invariants.
- [ ] Confirm whether current algorithm is correct or produce reproducer.
- [ ] Implement fix or explicitly constrain API if unsupported.

## ISSUE: Implement generic tower inversion (or explicit type restriction) for `D2 != 2`
Labels: feature, P1-important, area:ext-fields, area:traits-algebra
Milestone: v0.2-math-quality
Evidence:
- `src/structures/ext.rs:888-893` returns `None` for `D2 != 2`.
- Division wrappers in `src/structures/gf.rs:843` and `src/structures/gf.rs:866` panic in this case.
Impact:
- Correctness expectations and API completeness.
Definition of Done:
- [ ] Implement inversion algorithm for general `D2` OR restrict type-level API to `D2 = 2` where required.
- [ ] Ensure division errors are non-panicking in checked APIs.
- [ ] Add tests for supported and unsupported paths with explicit behavior.

## ISSUE: Resolve `ExtFieldWithModulus` doc/implementation mismatch
Labels: docs, P1-important, area:ext-fields, area:traits-algebra, area:docs-examples
Milestone: v0.2-math-quality
Evidence:
- `src/structures/ext.rs:653-657` claims the wrapper allows Ring/Field-style operations.
- Only constructor/accessor impls exist at `src/structures/ext.rs:668-682`; no `Ring`/`Field` impls for wrapper type.
Impact:
- API ergonomics and documentation trust.
Definition of Done:
- [ ] Either implement missing trait/operator support or update docs to reflect current scope.
- [ ] Add examples for the chosen path.
- [ ] Add tests validating documented behavior.

## ISSUE: Define and document constant-time guarantees per API surface
Labels: security, P1-important, area:constant-time, area:field-fp, area:docs-examples
Milestone: v0.2-math-quality
Evidence:
- `Cargo.toml` exposes `subtle` as constant-time feature.
- `src/structures/fp.rs:280-333` has ct helpers, but core ops include branches (`src/structures/fp.rs:913-930`) and variable-time inversion (`src/structures/fp.rs:982-1002`).
Impact:
- Constant-time safety and misuse risk in cryptographic contexts.
Definition of Done:
- [ ] Publish explicit constant-time contract by API/function class.
- [ ] Mark variable-time APIs prominently in docs.
- [ ] Add side-channel-focused review checklist/tests where feasible.

## ISSUE: Add strict irreducibility-checked deserialization for `GFWithModulus`
Labels: security, P1-important, area:serialization, area:ext-fields
Milestone: v0.2-math-quality
Evidence:
- `src/structures/gf.rs:1269-1271` notes `to_gf` skips irreducibility checks.
- `src/structures/gf.rs:1288` uses `Modulus::new_unchecked`.
Impact:
- Misuse hazard: reducible modulus can deserialize into non-field arithmetic context.
Definition of Done:
- [ ] Add strict API path using `Modulus::new` validation.
- [ ] Keep existing unchecked path explicitly named/documented for trusted inputs.
- [ ] Add tests for reducible modulus rejection in strict mode.

## ISSUE: Validate and correct binary-field polynomial presets and placeholders
Labels: bug, P1-important, area:ext-fields, area:docs-examples
Milestone: v0.2-math-quality
Evidence:
- `src/structures/gf.rs:1108-1114` contains “actually/approx/placeholder” comments for several entries.
- The table is presented as standard irreducible polynomial reference.
Impact:
- Correctness and documentation reliability for downstream binary-field work.
Definition of Done:
- [ ] Verify each listed polynomial against standard references.
- [ ] Correct inaccurate entries and remove placeholders/approximations.
- [ ] Add tests for known canonical entries (e.g., AES polynomial).

## ISSUE: Update README capability matrix and feature prerequisites
Labels: docs, P1-important, area:docs-examples, area:infra-ci
Milestone: v0.1-hardening
Evidence:
- README intro (`README.md:3-15`) still frames crate as minimal prime-field library.
- Current crate includes extension fields, NTT, factorization, discrete log, and multiple feature-gated APIs.
Impact:
- API ergonomics and onboarding accuracy.
Definition of Done:
- [ ] Add module capability matrix (Fp/Poly/Ext/GF/NTT) with feature requirements.
- [ ] Clarify which examples require which features.
- [ ] Document no_std and alloc limitations clearly.

## ISSUE: Expand property-based tests beyond `Fp`
Labels: tech-debt, P1-important, area:polynomials, area:ext-fields, area:ntt-fft
Milestone: v0.2-math-quality
Evidence:
- Property-based coverage is concentrated in `tests/field_properties.rs`.
- Poly/ext/ntt rely mostly on deterministic unit tests in module files.
Impact:
- Correctness confidence for algebraic invariants.
Definition of Done:
- [ ] Add proptests for polynomial ring laws and division invariants.
- [ ] Add proptests for extension-field arithmetic identities.
- [ ] Add NTT round-trip and convolution consistency proptests.

## ISSUE: Add criterion benchmarks for factorization and discrete-log paths
Labels: perf, P1-important, area:benchmarks, area:factorization, area:discrete-log
Milestone: v0.2-math-quality
Evidence:
- `benches/benchmarks.rs:10-181` lacks dedicated factorization and discrete-log benchmarks.
Impact:
- Performance visibility and regression detection.
Definition of Done:
- [ ] Add benchmark groups for `Poly::factor`/`roots` under representative sizes.
- [ ] Add benchmark groups for `Fp::discrete_log*` with varying subgroup sizes.
- [ ] Emit baseline reports suitable for PR comparison.

## ISSUE: Establish NTT vs naive crossover benchmarks and heuristics
Labels: perf, P1-important, area:benchmarks, area:ntt-fft, area:polynomials
Milestone: v0.2-math-quality
Evidence:
- Existing comparison in `benches/benchmarks.rs:102-114` is limited and does not produce actionable crossover guidance.
Impact:
- Performance and API ergonomics.
Definition of Done:
- [ ] Benchmark wider degree ranges and multiple primes.
- [ ] Record crossover points for `mul_ntt` vs naive multiplication.
- [ ] Document a selection heuristic for callers.

## ISSUE: Add deterministic RNG hooks for factorization/root APIs
Labels: feature, P2-later, area:randomness, area:factorization
Milestone: v0.2-math-quality
Evidence:
- Randomized factorization/root APIs are tied to `rand::Rng` (e.g., `src/structures/poly.rs:1095`, `src/structures/poly.rs:1213`, `src/structures/poly.rs:1276`).
Impact:
- Test reproducibility and portability for alternative RNG backends.
Definition of Done:
- [ ] Support deterministic seeded workflows in docs/examples.
- [ ] Add helpers or trait abstractions to reduce direct dependency friction.
- [ ] Add regression tests with fixed seeds.

## ISSUE: Clean warning set under non-default feature combinations
Labels: tech-debt, P1-important, area:field-fp, area:ext-fields, area:ntt-fft
Milestone: v0.1-hardening
Evidence:
- Warnings observed for unused imports in `src/structures/gf.rs:12`, `src/structures/ntt.rs:19`, and `src/structures/ext.rs:8` under different feature combinations.
Impact:
- Maintainability and CI signal quality.
Definition of Done:
- [ ] Eliminate feature-dependent unused-import warnings.
- [ ] Add targeted cfg-aware imports where needed.
- [ ] Keep clean `cargo test/build` output for supported feature sets.

## ISSUE: Roadmap: add multi-limb field backend for 256-bit primes
Labels: feature, P2-later, area:field-fp, area:traits-algebra
Milestone: v0.3-crypto-roadmap
Evidence:
- Current field core is `Fp<const P: u64>` with `u64/u128` arithmetic (`src/structures/fp.rs:90-99`, `src/structures/fp.rs:940-942`).
- No multi-limb representation exists in repo.
Impact:
- Capability expansion for modern cryptographic prime fields.
Definition of Done:
- [ ] Define limb-based field element type and trait integration strategy.
- [ ] Provide basic arithmetic and reduction for at least one 256-bit prime.
- [ ] Add conformance tests against known vectors.

## ISSUE: Roadmap: add dedicated quadratic extension helper `QuadExt<P, N>`
Labels: feature, P2-later, area:ext-fields, area:docs-examples
Milestone: v0.3-crypto-roadmap
Evidence:
- Extension arithmetic uses generic `ExtField`/`TowerField`; no specialized `QuadExt` type exists.
Impact:
- API ergonomics and potential performance wins for common degree-2 paths.
Definition of Done:
- [ ] Design `QuadExt` API around common cryptographic patterns.
- [ ] Provide conversions/interoperability with existing ext field types.
- [ ] Add docs/examples and comparative benchmarks.

## ISSUE: Roadmap: add additional discrete-log algorithms (Pollard rho/Pohlig-Hellman)
Labels: feature, P2-later, area:discrete-log, area:field-fp
Milestone: v0.3-crypto-roadmap
Evidence:
- `Fp` discrete log implementation documents BSGS (`src/structures/fp.rs:613-614`); no alternative algorithms present.
Impact:
- Capability and performance across larger subgroup scenarios.
Definition of Done:
- [ ] Add at least one additional algorithm with clear API boundaries.
- [ ] Provide correctness tests against known small/medium cases.
- [ ] Add benchmark comparisons vs BSGS.

## ISSUE: Roadmap: scaffold ECC module (curve traits + scalar multiplication)
Labels: feature, P2-later, area:traits-algebra, area:docs-examples
Milestone: v0.3-crypto-roadmap
Evidence:
- No ECC module currently exists in `src/` module layout.
Impact:
- Capability expansion toward general cryptography goals.
Definition of Done:
- [ ] Define foundational curve/scalar traits and module boundaries.
- [ ] Add one reference curve implementation skeleton with scalar multiplication.
- [ ] Add safety notes (constant-time/misuse constraints) and examples.
