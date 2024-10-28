use ark_ec::pairing::Pairing;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use std::ops::{Mul, Neg};

use crate::com::Com;

#[derive(Copy, Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct ExtractionKey<E: Pairing>(pub E::ScalarField, pub E::ScalarField);

impl<E: Pairing> ExtractionKey<E> {
    pub fn extract_1(&self, c: &Com<<E as Pairing>::G1>) -> E::G1Affine {
        (c.0.mul(&self.0.neg()) + c.1).into()
    }

    pub fn extract_2(&self, c: &Com<<E as Pairing>::G2>) -> E::G2Affine {
        (c.0.mul(&self.1.neg()) + c.1).into()
    }
}
