use proptest::prelude::*;

use kreep::{Field, Fp, Ring};

type F17 = Fp<17>;

fn arb_f17() -> impl Strategy<Value = F17> {
    (0u64..17).prop_map(F17::new)
}

fn arb_f17_nonzero() -> impl Strategy<Value = F17> {
    (1u64..17).prop_map(F17::new)
}

// ===== Addition properties =====

proptest! {
    #[test]
    fn addition_commutative(a in arb_f17(), b in arb_f17()) {
        prop_assert_eq!(a + b, b + a);
    }
}

proptest! {
    #[test]
    fn addition_associative(a in arb_f17(), b in arb_f17(), c in arb_f17()) {
        prop_assert_eq!((a + b) + c, a + (b + c));
    }
}

proptest! {
    #[test]
    fn additive_identity(a in arb_f17()) {
        prop_assert_eq!(a + F17::ZERO, a);
        prop_assert_eq!(F17::ZERO + a, a);
    }
}

proptest! {
    #[test]
    fn additive_inverse(a in arb_f17()) {
        prop_assert_eq!(a + (-a), F17::ZERO);
        prop_assert_eq!((-a) + a, F17::ZERO);
    }
}

proptest! {
    #[test]
    fn double_negation(a in arb_f17()) {
        prop_assert_eq!(-(-a), a);
    }
}

// ===== Subtraction properties =====

proptest! {
    #[test]
    fn subtraction_definition(a in arb_f17(), b in arb_f17()) {
        prop_assert_eq!(a - b, a + (-b));
    }
}

proptest! {
    #[test]
    fn subtraction_self_is_zero(a in arb_f17()) {
        prop_assert_eq!(a - a, F17::ZERO);
    }
}

// ===== Multiplication properties =====

proptest! {
    #[test]
    fn multiplication_commutative(a in arb_f17(), b in arb_f17()) {
        prop_assert_eq!(a * b, b * a);
    }
}

proptest! {
    #[test]
    fn multiplication_associative(a in arb_f17(), b in arb_f17(), c in arb_f17()) {
        prop_assert_eq!((a * b) * c, a * (b * c));
    }
}

proptest! {
    #[test]
    fn multiplicative_identity(a in arb_f17()) {
        prop_assert_eq!(a * F17::ONE, a);
        prop_assert_eq!(F17::ONE * a, a);
    }
}

proptest! {
    #[test]
    fn multiplicative_zero(a in arb_f17()) {
        prop_assert_eq!(a * F17::ZERO, F17::ZERO);
        prop_assert_eq!(F17::ZERO * a, F17::ZERO);
    }
}

// ===== Distributivity =====

proptest! {
    #[test]
    fn left_distributive(a in arb_f17(), b in arb_f17(), c in arb_f17()) {
        prop_assert_eq!(a * (b + c), a * b + a * c);
    }
}

proptest! {
    #[test]
    fn right_distributive(a in arb_f17(), b in arb_f17(), c in arb_f17()) {
        prop_assert_eq!((a + b) * c, a * c + b * c);
    }
}

// ===== Field properties (inverse) =====

proptest! {
    #[test]
    fn multiplicative_inverse(a in arb_f17_nonzero()) {
        let inv = a.inverse().unwrap();
        prop_assert_eq!(a * inv, F17::ONE);
        prop_assert_eq!(inv * a, F17::ONE);
    }
}

proptest! {
    #[test]
    fn double_inverse(a in arb_f17_nonzero()) {
        let inv = a.inverse().unwrap();
        let inv_inv = inv.inverse().unwrap();
        prop_assert_eq!(inv_inv, a);
    }
}

#[test]
fn inverse_of_one() {
    let one_inv = F17::ONE.inverse().unwrap();
    assert_eq!(one_inv, F17::ONE);
}

// ===== Division properties =====

proptest! {
    #[test]
    fn division_by_self(a in arb_f17_nonzero()) {
        prop_assert_eq!(a / a, F17::ONE);
    }
}

proptest! {
    #[test]
    fn division_consistency(a in arb_f17(), b in arb_f17_nonzero()) {
        prop_assert_eq!((a / b) * b, a);
    }
}

proptest! {
    #[test]
    fn division_definition(a in arb_f17(), b in arb_f17_nonzero()) {
        let b_inv = b.inverse().unwrap();
        prop_assert_eq!(a / b, a * b_inv);
    }
}

// ===== Value representation =====

proptest! {
    #[test]
    fn value_in_range(a in arb_f17()) {
        prop_assert!(a.value() < F17::modulus());
    }
}

proptest! {
    #[test]
    fn new_reduces_mod_p(v in 0u64..1000u64) {
        let a = F17::new(v);
        prop_assert_eq!(a.value(), v % 17);
    }
}

// Test zero has no inverse
#[test]
fn zero_has_no_inverse() {
    assert!(F17::ZERO.inverse().is_none());
}

// ===== Tests with a larger prime =====

mod larger_prime {
    use super::*;

    type F101 = Fp<101>;

    fn arb_f101() -> impl Strategy<Value = F101> {
        (0u64..101).prop_map(F101::new)
    }

    fn arb_f101_nonzero() -> impl Strategy<Value = F101> {
        (1u64..101).prop_map(F101::new)
    }

    proptest! {
        #[test]
        fn addition_commutative(a in arb_f101(), b in arb_f101()) {
            prop_assert_eq!(a + b, b + a);
        }
    }

    proptest! {
        #[test]
        fn distributive(a in arb_f101(), b in arb_f101(), c in arb_f101()) {
            prop_assert_eq!(a * (b + c), a * b + a * c);
        }
    }

    proptest! {
        #[test]
        fn multiplicative_inverse(a in arb_f101_nonzero()) {
            let inv = a.inverse().unwrap();
            prop_assert_eq!(a * inv, F101::ONE);
        }
    }

    proptest! {
        #[test]
        fn division_consistency(a in arb_f101(), b in arb_f101_nonzero()) {
            prop_assert_eq!((a / b) * b, a);
        }
    }
}
