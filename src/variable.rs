//! Defines the struct [Variable], the variable for the Pairing Product Equation.
//! i.e. the (`X`, `r`) and (`Y`, `s`) notated in the paper.

use ark_ec::CurveGroup;
use ark_std::rand::Rng;

use crate::Randomness;

/// Variable `X` or `Y` for the Pairing Product Equation, which represents the values that the prover wants to hide.
/// It carries the randomness `r` or `s` for being used in commitment scheme for proof construction.
#[derive(Copy, Clone, Debug)]
pub struct Variable<G: CurveGroup> {
    pub value: G::Affine,
    pub(crate) rand: Randomness<G>,
}

impl<G: CurveGroup> Variable<G> {
    /// Constructs a new variable `X` or `Y` with the given `value` and internal randomness `r` or `s`.
    pub fn new<R: Rng>(rng: &mut R, value: G::Affine) -> Self {
        Self::with_randomness(value, Randomness::rand(rng))
    }

    /// Constructs a new variable `X` or `Y` with the given `value` and internal randomness `r` or `s`,
    /// where the randomness is set to zero.
    pub fn with_zero_randomness(value: G::Affine) -> Self {
        Self::with_randomness(value, Randomness::zero())
    }

    /// Constructs a new variable `X` or `Y` with the given `value` and internal randomness `r` or `s`,
    /// where the randomness is set to the given randomness.
    pub fn with_randomness(value: G::Affine, rand: Randomness<G>) -> Self {
        Self { value, rand }
    }
}
