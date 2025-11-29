use core::ops::{Add, Mul, Neg, Sub};

/// A (not-necessarily commutative) ring.
///
/// This trait assumes:
/// - (R, +) is an abelian group with identity ZERO
/// - (R, ·) is a monoid with identity ONE
/// - multiplication distributes over addition.
pub trait Ring:
    Sized
    + Copy
    + Eq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
{
    const ZERO: Self;
    const ONE: Self;

    #[inline]
    fn is_zero(&self) -> bool {
        *self == Self::ZERO
    }
}
