use ark_ec::{CurveGroup, Group};
use ark_std::{rand::Rng, UniformRand};

#[derive(Clone, Debug)]
pub struct Variable<G: CurveGroup> {
    pub value: G::Affine,
    pub(crate) rand: VarRandomness<G>,
}

impl<G: CurveGroup> Variable<G> {
    pub fn new<R: Rng>(rng: &mut R, value: G::Affine) -> Self {
        let rand = VarRandomness(G::ScalarField::rand(rng), G::ScalarField::rand(rng));
        Self { value, rand }
    }
}

#[derive(Clone, Debug)]
pub struct VarRandomness<G: Group>(pub G::ScalarField, pub G::ScalarField);
