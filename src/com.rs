//! Defines the struct [Com], the `SXDH Commitments`` defined in section 6.2 in the paper [Fuc10](https://eprint.iacr.org/2010/233.pdf).

use ark_ec::CurveGroup;
use ark_std::rand::Rng;
use std::ops::Mul;

use crate::{commit::CommitmentKey, randomness::Randomness};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Com<G: CurveGroup>(pub G::Affine, pub G::Affine);

impl<G: CurveGroup> Com<G> {
    /// The commitment randomization function `RdCom(ck, c, r)`. Randomizes the commitment
    /// and return [ComRandomness] for Proof Adaption.
    ///
    /// The method modifies this commitment, while the original commitment is returned with randomness
    /// in the [ComRandomness].
    pub fn randomize<R: Rng>(&mut self, rng: &mut R, ck: &CommitmentKey<G>) -> ComRandomness<G> {
        let Randomness(r1, r2) = Randomness::<G>::rand(rng);
        let original = *self;

        // RdCom(ck, c, r)
        // = c * Com(ck, 0, r)
        // = (c1 + u11^r1 + u21^r2, c2 + u12^r1 + u22^r2)
        let a = ck.0 .0.mul(r1) + ck.1 .0.mul(r2);
        let b = ck.0 .1.mul(r1) + ck.1 .1.mul(r2);
        self.0 = (self.0 + a).into();
        self.1 = (self.1 + b).into();

        (original, Randomness(r1, r2))
    }
}

// TODO implement homonorphic properties of the commitment
// "Remark3. Comcommitments are homomorphic: Com(ck, X, r) + Com(ck, X', r') = Com(ck, X + X', r+r');"

/// A tuple of a commitment and its randomness. It is used in Proof Adaption as a the input
/// `(c, r)` or `(d, s)` in the proof adaption function `RdProof`.
pub type ComRandomness<G> = (Com<G>, Randomness<G>);
