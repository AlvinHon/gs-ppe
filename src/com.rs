use ark_ec::CurveGroup;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use std::ops::Mul;

use crate::{
    crs::CommitmentKey,
    variable::{VarRandomness, Variable},
};

#[derive(Copy, Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct Com<G: CurveGroup>(pub G::Affine, pub G::Affine);

impl<G: CurveGroup> Com<G> {
    pub fn new(ck: &CommitmentKey<G>, x: Variable<G>) -> Self {
        let VarRandomness(r1, r2) = x.rand;
        let x = x.value;

        let a = ck.0 .0.mul(r1) + ck.1 .0.mul(r2);
        let b = ck.0 .1.mul(r1) + ck.1 .1.mul(r2);
        Com(a.into(), (x + b).into())
    }

    pub fn randomize(&mut self, ck: &CommitmentKey<G>, r: VarRandomness<G>) {
        let VarRandomness(r1, r2) = r;

        let a = ck.0 .0.mul(r1) + ck.1 .0.mul(r2);
        let b = ck.0 .1.mul(r1) + ck.1 .1.mul(r2);
        self.0 = (self.0 + a).into();
        self.1 = (self.1 + b).into();
    }
}
