#![doc = include_str!("../README.md")]

pub mod com;
pub use com::Com;

pub mod commit;
pub use commit::CommitmentKeys;

pub mod equation;
pub use equation::Equation;

pub mod extract;
pub use extract::ExtractKey;

pub mod matrix;
pub use matrix::Matrix;

pub mod prove;
pub use prove::Proof;

pub mod randomness;
pub use randomness::Randomness;

pub mod variable;
pub use variable::Variable;

use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ff::Zero;
use ark_std::rand::Rng;
use std::ops::Mul;

pub fn setup<E: Pairing, R: Rng>(
    rng: &mut R,
    cks: &CommitmentKeys<E>,
    ay: &[(<E as Pairing>::G1Affine, Variable<<E as Pairing>::G2>)],
    xb: &[(Variable<<E as Pairing>::G1>, <E as Pairing>::G2Affine)],
    gamma: &Matrix<E::ScalarField>,
) -> ProofSystem<E> {
    assert_eq!(gamma.dim(), (xb.len(), ay.len()));

    let ay_product = ay.iter().fold(PairingOutput::zero(), |acc, (a, y)| {
        acc + E::pairing(a, y.value)
    });
    let xb_product = xb.iter().fold(PairingOutput::zero(), |acc, (x, b)| {
        acc + E::pairing(x.value, b)
    });

    let x: Vec<Variable<_>> = xb.iter().map(|(x, _)| x.clone()).collect();
    let y: Vec<Variable<_>> = ay.iter().map(|(_, y)| y.clone()).collect();

    let mut xy_product = PairingOutput::zero();
    for (j, y_j) in y.iter().enumerate() {
        for (i, x_i) in x.iter().enumerate() {
            xy_product += E::pairing(x_i.value, y_j.value).mul(gamma[(i, j)]);
        }
    }
    let target = ay_product + xb_product + xy_product;

    let a = ay.iter().map(|(a, _)| *a).collect();
    let b = xb.iter().map(|(_, b)| *b).collect();

    let equation = Equation::<E> {
        a,
        b,
        gamma: gamma.clone(),
        target,
    };
    let c = x.iter().map(|x_i| cks.u.commit(x_i)).collect();
    let d = y.iter().map(|y_i| cks.v.commit(y_i)).collect();
    let proof = Proof::new(rng, cks, &equation, &x, &y);
    ProofSystem {
        equation,
        c,
        d,
        proof,
    }
}

pub struct ProofSystem<E: Pairing> {
    pub equation: Equation<E>,
    pub c: Vec<Com<<E as Pairing>::G1>>,
    pub d: Vec<Com<<E as Pairing>::G2>>,
    pub proof: Proof<E>,
}
