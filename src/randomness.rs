//! Defines the struct [Randomness], the randomness commonly used in the entire scheme. i.e. the `r` and `s`
//! notated in the paper.

use std::ops::Add;

use ark_ec::Group;
use ark_ff::Zero;
use ark_std::{rand::Rng, UniformRand};

/// Randomness used in the entire scheme. i.e. the `r` and `s`.
#[derive(Copy, Clone, Debug)]
pub struct Randomness<G: Group>(pub G::ScalarField, pub G::ScalarField);

impl<G: Group> Randomness<G> {
    /// Generates a random `Randomness` using the given `rng`.
    pub fn rand<R: Rng>(rng: &mut R) -> Self {
        Self(G::ScalarField::rand(rng), G::ScalarField::rand(rng))
    }

    /// Returns a `Randomness` with both fields set to zero.
    pub fn zero() -> Self {
        Self(G::ScalarField::zero(), G::ScalarField::zero())
    }
}

impl<G: Group> Add for Randomness<G> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0, self.1 + other.1)
    }
}
