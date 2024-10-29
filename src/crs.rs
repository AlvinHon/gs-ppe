use ark_ec::{pairing::Pairing, CurveGroup};
use ark_std::{rand::Rng, One, UniformRand};
use std::ops::{Mul, Sub};

use crate::extract::ExtractionKey;

pub struct CRS<E: Pairing> {
    pub g1: <E as Pairing>::G1Affine,
    pub g2: <E as Pairing>::G2Affine,
    pub u: CommitmentKey<E::G1>,
    pub v: CommitmentKey<E::G2>,
}

pub struct CommitmentKey<G: CurveGroup>(pub (G::Affine, G::Affine), pub (G::Affine, G::Affine));

impl<E: Pairing> CRS<E> {
    fn new(
        g1: <E as Pairing>::G1Affine,
        g2: <E as Pairing>::G2Affine,
        a1: E::ScalarField,
        a2: E::ScalarField,
        t1: E::ScalarField,
        t2: E::ScalarField,
    ) -> Self {
        // u1 = (g1, g1^a1)
        let u1 = (g1, g1.mul(a1).into());
        // u2 = (g1^t1, g1^(a1*t1))
        let u2 = (g1.mul(t1).into(), u1.1.mul(t1).into());
        // v1 = (g2, g2^a2)
        let v1 = (g2, g2.mul(a2).into());
        // v2 = (g2^t2, g2^(a2*t2))
        let v2 = (g2.mul(t2).into(), v1.1.mul(t2).into());

        Self {
            g1,
            g2,
            u: CommitmentKey(u1, u2),
            v: CommitmentKey(v1, v2),
        }
    }

    fn new_wi(
        g1: <E as Pairing>::G1Affine,
        g2: <E as Pairing>::G2Affine,
        a1: E::ScalarField,
        a2: E::ScalarField,
        t1: E::ScalarField,
        t2: E::ScalarField,
    ) -> Self {
        // u1 = (g1, g1^a1)
        let u1 = (g1, g1.mul(a1).into());
        // u2 = (g1^t1, g1^(a1*t1 - 1))
        let u2 = (
            g1.mul(t1).into(),
            g1.mul(a1.mul(t1).sub(&E::ScalarField::one())).into(),
        );
        // v1 = (g2, g2^a2)
        let v1 = (g2, g2.mul(a2).into());
        // v2 = (g2^t2, g2^(a2*t2 - 1))
        let v2 = (
            g2.mul(t2).into(),
            g2.mul(a2.mul(t2).sub(&E::ScalarField::one())).into(),
        );
        Self {
            g1,
            g2,
            u: CommitmentKey(u1, u2),
            v: CommitmentKey(v1, v2),
        }
    }
}

pub fn setup<E: Pairing, R: Rng>(
    rng: &mut R,
    g1: <E as Pairing>::G1Affine,
    g2: <E as Pairing>::G2Affine,
) -> CRS<E> {
    let a1 = E::ScalarField::rand(rng);
    let a2 = E::ScalarField::rand(rng);
    let t1 = E::ScalarField::rand(rng);
    let t2 = E::ScalarField::rand(rng);
    CRS::new(g1, g2, a1, a2, t1, t2)
}

pub fn setup_ex<E: Pairing, R: Rng>(
    rng: &mut R,
    g1: <E as Pairing>::G1Affine,
    g2: <E as Pairing>::G2Affine,
) -> (CRS<E>, ExtractionKey<E>) {
    let a1 = E::ScalarField::rand(rng);
    let a2 = E::ScalarField::rand(rng);
    let t1 = E::ScalarField::rand(rng);
    let t2 = E::ScalarField::rand(rng);
    (CRS::new(g1, g2, a1, a2, t1, t2), ExtractionKey(a1, a2))
}

pub fn setup_wi<E: Pairing, R: Rng>(
    rng: &mut R,
    g1: <E as Pairing>::G1Affine,
    g2: <E as Pairing>::G2Affine,
) -> CRS<E> {
    let a1 = E::ScalarField::rand(rng);
    let a2 = E::ScalarField::rand(rng);
    let t1 = E::ScalarField::rand(rng);
    let t2 = E::ScalarField::rand(rng);
    CRS::new_wi(g1, g2, a1, a2, t1, t2)
}
