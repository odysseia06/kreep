/// Abstract group with a single binary operation, written multiplicatively.
///
/// Laws (you should test these for concrete types):
/// - associativity: (ab)c = a(bc)
/// - identity: e * a = a * e = a
/// - inverse: a * a⁻¹ = a⁻¹ * a = e
pub trait Group: Sized + Clone + Eq {
    /// Identity element `e`.
    fn identity() -> Self;

    /// Group operation `*` (often written as `·` or `+` in math).
    fn op(&self, rhs: &Self) -> Self;

    /// Inverse element a⁻¹.
    fn inverse(&self) -> Self;

    /// Exponentiation by a non-negative integer using square-and-multiply.
    fn pow(&self, exp: u64) -> Self {
        let mut base = self.clone();
        let mut result = Self::identity();

        let mut e = exp;
        while e > 0 {
            if e & 1 == 1 {
                result = result.op(&base);
            }
            base = base.op(&base);
            e >>= 1;
        }
        result
    }
}

/// Marker trait for abelian (commutative) groups.
///
/// For types implementing this, you *should* have:
///   a · b = b · a  for all a, b.
pub trait AbelianGroup: Group {}
