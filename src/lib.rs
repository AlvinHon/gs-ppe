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
use ark_std::{rand::Rng, Zero};
use std::ops::{Add, Mul};

/// Setup the proof system over the Pairing Product Equation:
///
/// ∏e(a, y) ∏e(x, b) ∏∏e(x, y)^gamma = T
///
/// where `a` and `b` are the public constants, `x` and `y` are the private variables,
/// and `gamma` is the matrix of the exponents.
///
/// It returns the Proof System that contains the equation, commitments, and the proof.
///
/// ## Panics
/// Panics if dimension of gamma does not match the length of `xb` and `ay`. i.e. gamma.dim() != (xb.len(), ay.len())
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

    let x: Vec<Variable<_>> = xb.iter().map(|(x, _)| *x).collect();
    let y: Vec<Variable<_>> = ay.iter().map(|(_, y)| *y).collect();

    let mut xy_product = PairingOutput::zero();
    for (j, y_j) in y.iter().enumerate() {
        for (i, x_i) in x.iter().enumerate() {
            xy_product += E::pairing(x_i.value, y_j.value).mul(gamma[(i, j)]);
        }
    }
    let target = ay_product + xb_product + xy_product;

    let a = ay.iter().map(|(a, _)| (*a).into()).collect();
    let b = xb.iter().map(|(_, b)| (*b).into()).collect();

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

/// The Proof System over the Pairing Product Equation. It consists of
/// - The specified pairing product `equation`.
/// - The commitments `c` and `d` which commit to the variables `x` and `y` respectively.
/// - The `proof` for proving `c` and `d` are committing to the variables `x` and `y` satisfying the `equation`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProofSystem<E: Pairing> {
    pub equation: Equation<E>,
    pub c: Vec<Com<<E as Pairing>::G1>>,
    pub d: Vec<Com<<E as Pairing>::G2>>,
    pub proof: Proof<E>,
}

impl<E: Pairing> ProofSystem<E> {
    /// Randomize the commitments `c` and `d` and the proof by applying the functions `RdCom` and `RdProof`
    /// define in the paper [Fuc10](https://eprint.iacr.org/2010/233.pdf).
    pub fn randomize<R: Rng>(mut self, rng: &mut R, cks: &CommitmentKeys<E>) -> Self {
        let cr = self
            .c
            .iter_mut()
            .map(|c_i| c_i.randomize(rng, &cks.u))
            .collect::<Vec<_>>();
        let ds = self
            .d
            .iter_mut()
            .map(|d_j| d_j.randomize(rng, &cks.v))
            .collect::<Vec<_>>();

        self.proof.randomize(rng, cks, &self.equation, &cr, &ds);
        self
    }
}

/// Homomorphic addition of two Proof Systems, defined in section 7.2 of the paper.
impl<E: Pairing> Add for ProofSystem<E> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let equation = self.equation + other.equation;
        let c = self.c.into_iter().chain(other.c).collect();
        let d = self.d.into_iter().chain(other.d).collect();
        let proof = self.proof + other.proof;
        ProofSystem {
            equation,
            c,
            d,
            proof,
        }
    }
}
