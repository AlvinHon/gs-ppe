//! Defines the struct [CommitmentKeys], the commitment key `ck` for `SXDH Commitments`` defined in section 6.2 in the paper [Fuc10](https://eprint.iacr.org/2010/233.pdf).

use ark_ec::{pairing::Pairing, CurveGroup};
use ark_std::{rand::Rng, One, UniformRand};
use std::ops::{Mul, Sub};

use crate::{com::Com, randomness::Randomness, variable::Variable, ExtractKey};

/// Contains commitment keys `u` and `v` for the `SXDH Commitments`, where
/// `u` and `v` belong to Group G1 and G2 respectively.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct CommitmentKeys<E: Pairing> {
    pub u: CommitmentKey<E::G1>,
    pub v: CommitmentKey<E::G2>,
}

impl<E: Pairing> CommitmentKeys<E> {
    /// Generates random commitment keys for standard setup of Commitment Scheme.
    pub fn rand<R: Rng>(rng: &mut R) -> CommitmentKeys<E> {
        let g1 = E::G1Affine::rand(rng);
        let g2 = E::G2Affine::rand(rng);
        Self::setup(rng, g1, g2)
    }

    /// Generates random commitment keys for standard setup of Commitment Scheme,
    /// given the generators `g1` and `g2`.
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

    /// Generates random commitment keys for standard setup of Commitment Scheme,
    /// and in addition, returns the extract key.
    pub fn rand_ex<R: Rng>(rng: &mut R) -> (CommitmentKeys<E>, ExtractKey<E>) {
        let g1 = E::G1Affine::rand(rng);
        let g2 = E::G2Affine::rand(rng);
        Self::setup_ex(rng, g1, g2)
    }

    /// Construct commitment keys for standard setup and in addition returns the extract key,
    /// given the generators `g1` and `g2`.
    pub fn setup_ex<R: Rng>(
        rng: &mut R,
        g1: <E as Pairing>::G1Affine,
        g2: <E as Pairing>::G2Affine,
    ) -> (CommitmentKeys<E>, ExtractKey<E>) {
        let a1 = E::ScalarField::rand(rng);
        let a2 = E::ScalarField::rand(rng);
        let t1 = E::ScalarField::rand(rng);
        let t2 = E::ScalarField::rand(rng);
        (
            CommitmentKeys::new(g1, g2, a1, a2, t1, t2),
            ExtractKey(a1, a2),
        )
    }

    /// Generates random commitment keys for perfectly hiding setup (also indistinguishable
    /// by SDXH) of Commitment Scheme.
    pub fn rand_wi<R: Rng>(rng: &mut R) -> CommitmentKeys<E> {
        let g1 = E::G1Affine::rand(rng);
        let g2 = E::G2Affine::rand(rng);
        Self::setup_wi(rng, g1, g2)
    }

    /// Generates random commitment keys for perfectly hiding setup (also indistinguishable
    /// by SDXH) of Commitment Scheme, given the generators `g1` and `g2`.
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

    /// Implements the `Setup` function in section 6.2 of the paper.
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

    /// Implements the `WISetup` function in section 6.2 of the paper.
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

/// The component in commitment keys, either `u` or `v` in [CommitmentKeys].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct CommitmentKey<G: CurveGroup>(pub (G::Affine, G::Affine), pub (G::Affine, G::Affine));

impl<G: CurveGroup> CommitmentKey<G> {
    /// The commitment function `Com`. Returns the commitment of the variable `x` or `y` according to
    /// which group G the commitment key belongs to.
    pub fn commit(&self, x: &Variable<G>) -> Com<G> {
        let Randomness(r1, r2) = x.rand;
        let x = x.value;

        // Com(ck, X, r) = (u11^r1 + u21^r2, x + u12^r1 + u22^r2)
        let a = self.0 .0.mul(r1) + self.1 .0.mul(r2);
        let b = self.0 .1.mul(r1) + self.1 .1.mul(r2);
        Com(a.into(), (x + b).into())
    }
}
