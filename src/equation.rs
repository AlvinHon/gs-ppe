use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ff::Zero;
use std::ops::Mul;

use crate::{arith::matrix::Matrix, com::Com, commit::CommitmentKey, prove::Proof};

pub struct Equation<E: Pairing> {
    pub a: Vec<<E as Pairing>::G1Affine>, // size = n
    pub b: Vec<<E as Pairing>::G2Affine>, // size = m
    pub gamma: Matrix<E::ScalarField>,    // dim = (m, n)
    pub target: PairingOutput<E>,
}

impl<E: Pairing> Equation<E> {
    pub fn verify(
        &self,
        u: &CommitmentKey<<E as Pairing>::G1>,
        v: &CommitmentKey<<E as Pairing>::G2>,
        c: &[Com<<E as Pairing>::G1>],
        d: &[Com<<E as Pairing>::G2>],
        proof: &Proof<E>,
    ) -> bool {
        let (m, n) = self.gamma.dim();
        if c.len() != n || d.len() != m || proof.phi.dim() != (2, 2) || proof.theta.dim() != (2, 2)
        {
            return false;
        }

        // Check Equation 1
        let lhs = c
            .iter()
            .enumerate()
            .fold(PairingOutput::zero(), |acc, (i, c_i)| {
                let d_product = d
                    .iter()
                    .enumerate()
                    .fold(<E as Pairing>::G2::zero(), |acc, (j, d_j)| {
                        acc + d_j.0.mul(self.gamma[(i, j)])
                    })
                    .into();

                acc + E::pairing(c_i.0, d_product)
            });
        let rhs = E::pairing(u.0 .0, proof.phi[(0, 0)])
            + E::pairing(u.1 .0, proof.phi[(1, 0)])
            + E::pairing(proof.theta[(0, 0)], v.0 .0)
            + E::pairing(proof.theta[(1, 0)], v.1 .0);

        if lhs != rhs {
            return false;
        }

        // Check Equation 2
        let lhs = c
            .iter()
            .enumerate()
            .fold(PairingOutput::zero(), |acc, (i, c_i)| {
                let d_product = d
                    .iter()
                    .enumerate()
                    .fold(<E as Pairing>::G2::zero(), |acc, (j, d_j)| {
                        acc + d_j.1.mul(self.gamma[(i, j)])
                    })
                    .into();

                acc + E::pairing(c_i.0, self.b[i] + d_product)
            });
        let rhs = E::pairing(u.0 .0, proof.phi[(0, 1)])
            + E::pairing(u.1 .0, proof.phi[(1, 1)])
            + E::pairing(proof.theta[(0, 0)], v.0 .1)
            + E::pairing(proof.theta[(1, 0)], v.1 .1);
        if lhs != rhs {
            return false;
        }

        // Check Equation 3
        let lhs = d
            .iter()
            .enumerate()
            .fold(PairingOutput::zero(), |acc, (j, d_j)| {
                let c_product = c
                    .iter()
                    .enumerate()
                    .fold(<E as Pairing>::G1::zero(), |acc, (i, c_i)| {
                        acc + c_i.1.mul(self.gamma[(i, j)])
                    })
                    .into();

                acc + E::pairing(self.a[j] + c_product, d_j.0)
            });
        let rhs = E::pairing(u.0 .1, proof.phi[(0, 0)])
            + E::pairing(u.1 .1, proof.phi[(1, 0)])
            + E::pairing(proof.theta[(0, 1)], v.0 .0)
            + E::pairing(proof.theta[(1, 1)], v.1 .0);
        if lhs != rhs {
            return false;
        }

        // Check Equation 4
        let lhs = {
            let a_d = self
                .a
                .iter()
                .zip(d.iter())
                .fold(PairingOutput::zero(), |acc, (a, d)| {
                    acc + E::pairing(a, d.1)
                });
            // TODO: optimize this - the 'bd' part was computed before
            let c_bd = c
                .iter()
                .enumerate()
                .fold(PairingOutput::zero(), |acc, (i, c_i)| {
                    let d_product = d
                        .iter()
                        .enumerate()
                        .fold(<E as Pairing>::G2::zero(), |acc, (j, d_j)| {
                            acc + d_j.1.mul(self.gamma[(i, j)])
                        })
                        .into();

                    acc + E::pairing(c_i.1, self.b[i] + d_product)
                });
            a_d + c_bd
        };
        let rhs = self.target
            + E::pairing(u.0 .1, proof.phi[(0, 1)])
            + E::pairing(u.1 .1, proof.phi[(1, 1)])
            + E::pairing(proof.theta[(0, 1)], v.0 .1)
            + E::pairing(proof.theta[(1, 1)], v.1 .1);

        lhs != rhs
    }
}
