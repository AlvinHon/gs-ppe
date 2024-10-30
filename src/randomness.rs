use std::ops::Add;

use ark_ec::Group;
use ark_std::{rand::Rng, UniformRand};

#[derive(Clone, Debug)]
pub struct Randomness<G: Group>(pub G::ScalarField, pub G::ScalarField);

impl<G: Group> Randomness<G> {
    pub fn rand<R: Rng>(rng: &mut R) -> Self {
        Self(G::ScalarField::rand(rng), G::ScalarField::rand(rng))
    }
}

impl<G: Group> Add for Randomness<G> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0, self.1 + other.1)
    }
}
