# Roadmap and API Design

This roadmap proposes concrete, incremental features with API sketches.
Ordering balances impact, correctness, and developer effort.

## Phase 1: Safer Modulus Handling (High Priority)
**Goal:** Prevent mixed-modulus errors and make extension fields easier to use.

**API design:**
```rust
pub enum ModulusError {
    WrongDegree { expected: usize, got: usize },
    NotMonic,
    NotIrreducible,
}

pub struct Modulus<const P: u64, const D: usize> {
    poly: Poly<P>,
    // Potentially cache validation results if needed.
}

impl<const P: u64, const D: usize> Modulus<P, D> {
    // Requires irreducibility check; may be feature-gated.
    pub fn new(poly: Poly<P>) -> Result<Self, ModulusError>;

    // Skips irreducibility check, still validates degree and monic.
    pub fn new_unchecked(poly: Poly<P>) -> Result<Self, ModulusError>;

    pub fn poly(&self) -> &Poly<P>;
}

pub struct GF<const P: u64, const D: usize> {
    elem: ExtField<P, D>,
    modulus: Rc<Modulus<P, D>>,
}

impl<const P: u64, const D: usize> GF<P, D> {
    pub fn new(coeffs: [Fp<P>; D], modulus: Rc<Modulus<P, D>>) -> Self;
    pub fn modulus(&self) -> &Modulus<P, D>;
}
```

**Validation rules:**
- Degree must be exactly `D`.
- Polynomial must be monic.
- Polynomial must be irreducible (feature-gated or optional).

**Work items:**
- Add `ModulusError` and `Modulus` type.
- Decide whether irreducibility requires `rand` or is optional.
- Refactor `GF` to accept `Modulus`.
- Update examples/tests to use validated modulus.

---

## Phase 2: Serialization for Core Types
**Goal:** Provide consistent `serde` support across types.

**API design:**
```rust
// Feature: serde
impl<const P: u64> Serialize for Poly<P> { /* coeffs */ }
impl<const P: u64, const D: usize> Serialize for ExtField<P, D> { /* coeff array */ }
impl<const P: u64, const D: usize> Serialize for GF<P, D> { /* coeffs only */ }

#[derive(Serialize, Deserialize)]
pub struct GFWithModulus<const P: u64, const D: usize> {
    coeffs: [u64; D],
    modulus: [u64; D + 1],
}
```

**Notes:**
- Default `GF` serialization is coefficients-only; caller provides modulus.
- Optional `GFWithModulus` for self-contained data.
- Document that deserializing `GF` requires the same modulus.

---

## Phase 3: NTT Precomputation API
**Goal:** Improve repeated polynomial multiplication throughput.

**API design:**
```rust
pub struct NttPlan<const P: u64> {
    n: usize,
    omega: Fp<P>,
    omega_inv: Fp<P>,
    n_inv: Fp<P>,
    twiddles: Vec<Fp<P>>,
    inv_twiddles: Vec<Fp<P>>,
}

impl<const P: u64> NttPlan<P> {
    pub fn new(n: usize) -> Result<Self, &'static str>;
    pub fn ntt(&self, a: &mut [Fp<P>]);
    pub fn intt(&self, a: &mut [Fp<P>]);
    pub fn mul_ntt_with_plan(&self, a: &Poly<P>, b: &Poly<P>) -> Poly<P>;
    // Accumulates acc += a * b where a/b are already in NTT domain.
    pub fn mul_accumulate(&self, acc: &mut [Fp<P>], a: &[Fp<P>], b: &[Fp<P>]);
}
```

**Work items:**
- Create `NttPlan` with cached twiddles.
- Provide plan-based multiplication helpers.

---

## Phase 4: Polynomial Utilities
**Goal:** Add targeted helpers that reduce boilerplate.

**API design ideas:**
```rust
impl<const P: u64> Poly<P> {
    pub fn div_by_linear(&self, r: Fp<P>) -> (Self, Fp<P>); // quotient + remainder
    pub fn eval_many(&self, xs: &[Fp<P>]) -> Vec<Fp<P>>;
    pub fn interpolate(points: &[(Fp<P>, Fp<P>)]) -> Option<Self>;
    pub fn from_roots(roots: &[Fp<P>]) -> Self;
}
```

**Notes:**
- Prefer `div_by_linear` or `synthetic_div` naming to avoid confusion.
- `interpolate` name is sufficient; "many" is implied by the slice.

---

## Phase 5: `no_std` + `alloc` Support
**Goal:** Make the core usable in embedded/crypto contexts.

**Approach:**
- `Fp` and fixed-size types (`ExtField`) are `no_std` without `alloc`.
- `Poly` requires `alloc` for `Vec`.
- Add `no_std` and `alloc` feature gates; document constraints.

---

## Phase 6: Extension Field Extras
**Goal:** Deeper field theory utilities for extensions.

**API design ideas:**
```rust
impl<const P: u64, const D: usize> ExtField<P, D> {
    pub fn minimal_polynomial(&self, modulus: &Poly<P>) -> Poly<P>;
    pub fn conjugates(&self, modulus: &Poly<P>) -> Vec<Self>;
    pub fn is_primitive(&self, modulus: &Poly<P>) -> bool;
    pub fn order(&self, modulus: &Poly<P>) -> u64;
}
```

---

## Cross-Cutting Items (Suggested Mini-Phases)
1. **Error handling strategy**: align on `Result` vs panics for public APIs.
2. **Documentation**: add `#[doc]` examples for public types and methods.
3. **Benchmarks**: measure NTT and polynomial factorization improvements.
4. **TowerField parity**: align `TowerField` ergonomics with `GF`.

---

## Deliverable Checklist (Suggested Order)
1. Modulus type + `GF` refactor
2. `serde` support for `Poly`/`ExtField`/`GF`
3. `NttPlan` and plan-based multiplication
4. Polynomial helpers (`div_by_linear`, `eval_many`)
5. `no_std`/`alloc` split
6. Extension-field extras (min poly, conjugates, order)
7. Cross-cutting items (errors, docs, benchmarks, TowerField)
