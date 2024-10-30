use ark_ec::CurveGroup;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::Rng;
use std::ops::Mul;

use crate::{commit::CommitmentKey, randomness::Randomness};

#[derive(Copy, Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Com<G: CurveGroup>(pub G::Affine, pub G::Affine);

impl<G: CurveGroup> Com<G> {
    pub fn randomize<R: Rng>(&mut self, rng: &mut R, ck: &CommitmentKey<G>) -> ComRandomness<G> {
        let Randomness(r1, r2) = Randomness::<G>::rand(rng);
        let original = *self;

        let a = ck.0 .0.mul(r1) + ck.1 .0.mul(r2);
        let b = ck.0 .1.mul(r1) + ck.1 .1.mul(r2);
        self.0 = (self.0 + a).into();
        self.1 = (self.1 + b).into();

        (original, Randomness(r1, r2))
    }
}

pub type ComRandomness<G> = (Com<G>, Randomness<G>);
