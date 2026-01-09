//! Finite field arithmetic library.
//!
//! This crate provides types and operations for finite field arithmetic,
//! including prime fields, extension fields, and polynomials.
//!
//! # Features
//!
//! - `std` (default): Enables standard library support. Implies `alloc`.
//! - `alloc`: Enables types that require heap allocation (`Poly`, `GF`, `NTT`).
//! - `serde`: Enables serialization support.
//! - `rand`: Enables random number generation support.
//! - `subtle`: Enables constant-time operations.
//!
//! # `no_std` Support
//!
//! This crate supports `no_std` environments. Without the `alloc` feature,
//! only `Fp` and `ExtField` (with fixed-size coefficients) are available.
//! With `alloc` but without `std`, all types are available but some
//! functionality may be limited.

#![no_std]

#[cfg(feature = "std")]
extern crate std;

#[cfg(feature = "alloc")]
extern crate alloc;

pub mod algebra;
pub mod structures;
pub mod utils;

pub use algebra::field::Field;
pub use algebra::group::Group;
pub use algebra::ring::Ring;

pub use structures::ext::ExtField;
pub use structures::ext::TowerField;
pub use structures::fp::Fp;
#[cfg(feature = "alloc")]
pub use structures::gf;
#[cfg(all(feature = "serde", feature = "alloc"))]
pub use structures::gf::GFWithModulus;
#[cfg(feature = "alloc")]
pub use structures::gf::{Modulus, ModulusError};
#[cfg(feature = "alloc")]
pub use structures::ntt;
#[cfg(feature = "alloc")]
pub use structures::poly::Poly;
pub use utils::{ceil_sqrt_u64, gcd, is_prime};
