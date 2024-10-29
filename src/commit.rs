use ark_ec::{pairing::Pairing, CurveGroup};
use ark_ff::One;
use ark_std::{rand::Rng, UniformRand};
use std::ops::{Mul, Sub};

use crate::{
    com::Com,
    variable::{VarRandomness, Variable},
    ExtractionKey,
};

pub struct CommitmentKeys<E: Pairing> {
    pub u: CommitmentKey<E::G1>,
    pub v: CommitmentKey<E::G2>,
}

impl<E: Pairing> CommitmentKeys<E> {
    pub fn rand<R: Rng>(rng: &mut R) -> CommitmentKeys<E> {
        let g1 = E::G1Affine::rand(rng);
        let g2 = E::G2Affine::rand(rng);
        Self::setup(rng, g1, g2)
    }

    pub fn setup<R: Rng>(
        rng: &mut R,
        g1: <E as Pairing>::G1Affine,
        g2: <E as Pairing>::G2Affine,
    ) -> CommitmentKeys<E> {
        let a1 = E::ScalarField::rand(rng);
        let a2 = E::ScalarField::rand(rng);
        let t1 = E::ScalarField::rand(rng);
        let t2 = E::ScalarField::rand(rng);
        CommitmentKeys::new(g1, g2, a1, a2, t1, t2)
    }

    pub fn rand_ex<R: Rng>(rng: &mut R) -> (CommitmentKeys<E>, ExtractionKey<E>) {
        let g1 = E::G1Affine::rand(rng);
        let g2 = E::G2Affine::rand(rng);
        Self::setup_ex(rng, g1, g2)
    }

    pub fn setup_ex<R: Rng>(
        rng: &mut R,
        g1: <E as Pairing>::G1Affine,
        g2: <E as Pairing>::G2Affine,
    ) -> (CommitmentKeys<E>, ExtractionKey<E>) {
        let a1 = E::ScalarField::rand(rng);
        let a2 = E::ScalarField::rand(rng);
        let t1 = E::ScalarField::rand(rng);
        let t2 = E::ScalarField::rand(rng);
        (
            CommitmentKeys::new(g1, g2, a1, a2, t1, t2),
            ExtractionKey(a1, a2),
        )
    }

    pub fn rand_wi<R: Rng>(rng: &mut R) -> CommitmentKeys<E> {
        let g1 = E::G1Affine::rand(rng);
        let g2 = E::G2Affine::rand(rng);
        Self::setup_wi(rng, g1, g2)
    }

    pub fn setup_wi<R: Rng>(
        rng: &mut R,
        g1: <E as Pairing>::G1Affine,
        g2: <E as Pairing>::G2Affine,
    ) -> CommitmentKeys<E> {
        let a1 = E::ScalarField::rand(rng);
        let a2 = E::ScalarField::rand(rng);
        let t1 = E::ScalarField::rand(rng);
        let t2 = E::ScalarField::rand(rng);
        CommitmentKeys::new_wi(g1, g2, a1, a2, t1, t2)
    }

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
            u: CommitmentKey(u1, u2),
            v: CommitmentKey(v1, v2),
        }
    }
}

pub struct CommitmentKey<G: CurveGroup>(pub (G::Affine, G::Affine), pub (G::Affine, G::Affine));

impl<G: CurveGroup> CommitmentKey<G> {
    pub fn commit(&self, x: Variable<G>) -> Com<G> {
        let VarRandomness(r1, r2) = x.rand;
        let x = x.value;

        let a = self.0 .0.mul(r1) + self.1 .0.mul(r2);
        let b = self.0 .1.mul(r1) + self.1 .1.mul(r2);
        Com(a.into(), (x + b).into())
    }

    pub fn randomize(&self, com: &mut Com<G>, r: VarRandomness<G>) {
        let VarRandomness(r1, r2) = r;

        let a = self.0 .0.mul(r1) + self.1 .0.mul(r2);
        let b = self.0 .1.mul(r1) + self.1 .1.mul(r2);
        com.0 = (com.0 + a).into();
        com.1 = (com.1 + b).into();
    }
}
