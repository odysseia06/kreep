//! Serde serialization/deserialization tests
//!
//! Run with: cargo test --features serde --test serde_tests

#![cfg(feature = "serde")]

use kreep::gf::{irreducible_poly_deg2, GFWithModulus, Modulus, GF};
use kreep::{ExtField, Fp, Poly};
use std::rc::Rc;

type F17 = Fp<17>;

#[test]
fn fp_roundtrip() {
    let a = F17::new(7);
    let json = serde_json::to_string(&a).unwrap();
    assert_eq!(json, "7");
    let b: F17 = serde_json::from_str(&json).unwrap();
    assert_eq!(a, b);
}

#[test]
fn poly_roundtrip() {
    // 3 + 2x + x^2
    let p = Poly::new(vec![F17::new(3), F17::new(2), F17::new(1)]);
    let json = serde_json::to_string(&p).unwrap();
    assert_eq!(json, "[3,2,1]");
    let q: Poly<17> = serde_json::from_str(&json).unwrap();
    assert_eq!(p, q);
}

#[test]
fn poly_zero_roundtrip() {
    let p = Poly::<17>::zero();
    let json = serde_json::to_string(&p).unwrap();
    assert_eq!(json, "[]");
    let q: Poly<17> = serde_json::from_str(&json).unwrap();
    assert_eq!(p, q);
}

#[test]
fn extfield_roundtrip() {
    let a = ExtField::<17, 2>::new([F17::new(2), F17::new(3)]);
    let json = serde_json::to_string(&a).unwrap();
    assert_eq!(json, "[2,3]");
    let b: ExtField<17, 2> = serde_json::from_str(&json).unwrap();
    assert_eq!(a, b);
}

#[test]
fn extfield_wrong_length_fails() {
    let json = "[1,2,3]"; // 3 elements, but ExtField<17, 2> expects 2
    let result: Result<ExtField<17, 2>, _> = serde_json::from_str(json);
    assert!(result.is_err());
}

#[test]
fn gf_serialize() {
    let modulus = Rc::new(Modulus::<17, 2>::new_unchecked(irreducible_poly_deg2::<17>()).unwrap());
    let a = GF::<17, 2>::new([F17::new(5), F17::new(7)], modulus);
    let json = serde_json::to_string(&a).unwrap();
    assert_eq!(json, "[5,7]");
}

#[test]
fn gf_with_modulus_roundtrip() {
    let modulus = Rc::new(Modulus::<17, 2>::new_unchecked(irreducible_poly_deg2::<17>()).unwrap());
    let a = GF::<17, 2>::new([F17::new(5), F17::new(7)], modulus);

    // Convert to self-contained form
    let with_mod = GFWithModulus::from_gf(&a);
    let json = serde_json::to_string(&with_mod).unwrap();

    // Deserialize back
    let with_mod2: GFWithModulus<17, 2> = serde_json::from_str(&json).unwrap();
    let b = with_mod2.to_gf().unwrap();

    assert_eq!(a, b);
}

#[test]
fn gf_with_modulus_preserves_modulus() {
    let modulus = Rc::new(Modulus::<17, 2>::new_unchecked(irreducible_poly_deg2::<17>()).unwrap());
    let a = GF::<17, 2>::new([F17::new(5), F17::new(7)], modulus);

    let with_mod = GFWithModulus::from_gf(&a);
    let json = serde_json::to_string(&with_mod).unwrap();

    // Check that modulus is included
    assert!(json.contains("modulus"));
    assert!(json.contains("coeffs"));

    let with_mod2: GFWithModulus<17, 2> = serde_json::from_str(&json).unwrap();
    let b = with_mod2.to_gf().unwrap();

    // Verify the modulus was correctly restored
    assert_eq!(a.modulus().poly(), b.modulus().poly());
}
