use super::ring::Ring;

/// A (commutative) field.
///
/// Extends `Ring` with multiplicative inverses for all non-zero elements.
pub trait Field: Ring {
    /// Multiplicative inverse `a⁻¹`, if it exists.
    ///
    /// For a true field:
    /// - `self == ZERO`  ⇒  `None`
    /// - otherwise       ⇒  `Some(a⁻¹)`
    fn inverse(self) -> Option<Self>;

    /// Safe division: returns `None` on division by zero.
    #[inline]
    fn try_div(self, rhs: Self) -> Option<Self> {
        rhs.inverse().map(|inv| self * inv)
    }
}
