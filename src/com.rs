use ark_ec::{CurveGroup, Group};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use std::ops::Mul;

use crate::crs::CK;

#[derive(Copy, Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct Com<G: CurveGroup>(pub G::Affine, pub G::Affine);

impl<G: CurveGroup> Com<G> {
    pub fn new(
        ck: &CK<G>,
        x: G::Affine,
        r1: <G as Group>::ScalarField,
        r2: <G as Group>::ScalarField,
    ) -> Self {
        let a = ck.0 .0.mul(r1) + ck.1 .0.mul(r2);
        let b = ck.0 .1.mul(r1) + ck.1 .1.mul(r2);
        Com(a.into(), (x + b).into())
    }
    pub fn randomize(
        &mut self,
        ck: &CK<G>,
        r1: <G as Group>::ScalarField,
        r2: <G as Group>::ScalarField,
    ) {
        let a = ck.0 .0.mul(r1) + ck.1 .0.mul(r2);
        let b = ck.0 .1.mul(r1) + ck.1 .1.mul(r2);
        self.0 = (self.0 + a).into();
        self.1 = (self.1 + b).into();
    }
}
