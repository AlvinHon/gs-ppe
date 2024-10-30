use ark_ec::CurveGroup;
use ark_std::rand::Rng;

use crate::Randomness;

#[derive(Clone, Debug)]
pub struct Variable<G: CurveGroup> {
    pub value: G::Affine,
    pub(crate) rand: Randomness<G>,
}

impl<G: CurveGroup> Variable<G> {
    pub fn new<R: Rng>(rng: &mut R, value: G::Affine) -> Self {
        let rand = Randomness::rand(rng);
        Self { value, rand }
    }
}
