use ark_ec::pairing::{Pairing, PairingOutput};

use crate::arith::matrix::Matrix;

pub struct Equation<E: Pairing> {
    pub a: Vec<<E as Pairing>::G1Affine>, // size = n
    pub b: Vec<<E as Pairing>::G2Affine>, // size = m
    pub gamma: Matrix<E::ScalarField>,    // dim = (m, n)
    pub target: PairingOutput<E>,
}
