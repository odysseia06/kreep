# Example Ideas and Implementation Details

This document outlines example programs that demonstrate the library's features
and provides implementation notes for each. These are intended for
`examples/*.rs` files.

## 1) Finite Field Basics
**Goal:** Show core `Fp<P>` operations and number theory helpers.

**Features covered:** `add`, `mul`, `pow`, `inverse`, `legendre`, `sqrt`, `discrete_log`.

**Implementation outline:**
- Define `type F17 = Fp<17>`.
- Create a few elements, show arithmetic and inverses.
- Compute Legendre symbol and square roots.
- Demonstrate discrete log in a small field.

**Example sketch:**
- Pick `base = F17::new(3)` (primitive root for 17).
- Compute `target = base.pow(5)`.
- Print `target.discrete_log(base)`.

## 2) Polynomial Arithmetic
**Goal:** Demonstrate `Poly<P>` operations and Euclidean algorithms.

**Features covered:** `add`, `mul`, `div_rem`, `gcd`, `extended_gcd`, `derivative`.

**Implementation outline:**
- Build two polynomials with `Poly::new`.
- Show multiplication and division with remainder.
- Compute gcd and extended gcd, verify `s*a + t*b = g`.
- Compute derivative and evaluate at a point.

**Example sketch:**
- `p(x) = 1 + 2x + x^2`, `q(x) = 1 + x`.
- Print `(p * q)` and `p.div_rem(&q)`.

## 3) Irreducible and Primitive Polynomials
**Goal:** Show irreducibility/primitive checks and generators.

**Features covered:** `is_irreducible`, `is_primitive`, `find_irreducible`.

**Implementation outline:**
- Use a known reducible polynomial and an irreducible one.
- Call `is_irreducible` and `is_primitive`.
- Call `find_irreducible::<P>(degree)` and show its degree.

**Example sketch:**
- Over `Fp<17>`, compare `x^2 - 3` and `x^2 - 4`.
- Use `find_irreducible::<5>(4)`.

## 4) Extension Fields with `ExtField`
**Goal:** Show manual extension field arithmetic with explicit modulus.

**Features covered:** `mul_mod`, `inverse_mod`, `pow_mod`, `norm`, `trace`.

**Implementation outline:**
- Build modulus `m(x)` as an irreducible polynomial.
- Create `ExtField::<P, D>` elements.
- Multiply and invert using `mul_mod` and `inverse_mod`.
- Compute `norm` and `trace`.

**Example sketch:**
- `F17[x]/(x^2 - 3)` with `ExtField<17, 2>`.
- `a = 2 + 3x`, `b = 1 + x`.

## 5) `GF` Wrapper for Natural Arithmetic
**Goal:** Show extension field arithmetic without passing modulus explicitly.

**Features covered:** `GF` operations, `frobenius`, `pow`, `norm`, `trace`.

**Implementation outline:**
- Create modulus with `gf::irreducible_poly_deg2::<17>()`.
- Build `GF::<17, 2>` values sharing the modulus.
- Use `+`, `*`, `/`, `pow`, `frobenius`.

**Example sketch:**
- Construct `a` and `b`, compute `(a * b) / a`.
- Show `a.frobenius()` and `a.trace()`.

## 6) NTT-Based Polynomial Multiplication
**Goal:** Compare `mul_ntt` with naive multiplication.

**Features covered:** `ntt_info`, `mul_ntt`.

**Implementation outline:**
- Use NTT-friendly prime `998244353`.
- Multiply two polynomials with `mul_ntt` and `Poly` multiplication.
- Assert they match.

**Example sketch:**
- `a = 1 + 2x + 3x^2`, `b = 4 + 5x`.
- Print both results and equality.

## 7) Factorization and Roots (rand feature)
**Goal:** Showcase polynomial factorization and root finding.

**Features covered:** `factor`, `roots` (feature `rand`).

**Implementation outline:**
- Enable `rand` feature for the example.
- Build a polynomial with known linear factors and multiplicities.
- Use `factor` and `roots` to recover them.

**Example sketch:**
- `f(x) = (x - 1)^2 (x - 2)` over `Fp<17>`.
- Print factor list and root multiplicities.

## 8) Square Roots and Quadratic Residues
**Goal:** Show `sqrt` and residue checks in a small prime field.

**Features covered:** `legendre`, `is_quadratic_residue`, `sqrt`.

**Implementation outline:**
- List all quadratic residues mod 17.
- For each, compute and verify a square root.
- Demonstrate a non-residue returning `None`.

**Example sketch:**
- Scan `0..17`, collect residues and square roots.

## 9) Discrete Log and Multiplicative Order
**Goal:** Demonstrate BSGS discrete log and order computation.

**Features covered:** `discrete_log`, `discrete_log_with_order`, `multiplicative_order`.

**Implementation outline:**
- Pick a primitive root `g` and compute order.
- Compute discrete logs for a few elements.
- Show behavior in a proper subgroup via `discrete_log_with_order`.

**Example sketch:**
- In `Fp<17>`, `g = 3` (order 16), show logs for `g^k`.

## 10) Property Checks (lightweight)
**Goal:** Provide a small runnable demo that checks key algebraic identities.

**Features covered:** `Ring`/`Field` identities, basic invariants.

**Implementation outline:**
- Generate a few elements in `Fp<17>`.
- Check distributivity and inverse properties in a loop.

**Example sketch:**
- Loop over `0..17`, assert `(a / b) * b == a` for nonzero `b`.

## Suggested File Map
- `examples/field_basics.rs`
- `examples/poly_arithmetic.rs`
- `examples/irreducible_polys.rs`
- `examples/extfield.rs`
- `examples/gf_wrapper.rs`
- `examples/ntt_mul.rs`
- `examples/factorization.rs` (requires `--features rand`)
- `examples/quadratic_residues.rs`
- `examples/discrete_log.rs`
- `examples/field_identities.rs`
