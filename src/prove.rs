use ark_ec::pairing::Pairing;
use ark_ff::Zero;
use ark_std::rand::Rng;
use std::ops::{Mul, Neg};

use crate::{
    arith::matrix::Matrix, commit::CommitmentKeys, equation::Equation, variable::Variable,
};

pub struct Proof<E: Pairing> {
    pub(crate) phi: Matrix<<E as Pairing>::G2Affine>,
    pub(crate) theta: Matrix<<E as Pairing>::G1Affine>,
}

pub fn prove<E: Pairing, R: Rng>(
    rng: &mut R,
    crs: &CommitmentKeys<E>,
    equ: &Equation<E>,
    x: &[Variable<<E as Pairing>::G1>],
    y: &[Variable<<E as Pairing>::G2>],
) -> Proof<E> {
    assert_eq!(equ.a.len(), x.len());
    assert_eq!(equ.b.len(), y.len());

    let z = Matrix::<E::ScalarField>::rand(rng, 2, 2);
    let z_u = Matrix::new(&[
        [
            crs.u.0 .0.mul(z[(0, 0)]) + crs.u.1 .0.mul(z[(0, 1)]),
            crs.u.0 .1.mul(z[(0, 0)]) + crs.u.1 .1.mul(z[(0, 1)]),
        ],
        [
            crs.u.0 .0.mul(z[(1, 0)]) + crs.u.1 .0.mul(z[(1, 1)]),
            crs.u.0 .1.mul(z[(1, 0)]) + crs.u.1 .1.mul(z[(1, 1)]),
        ],
    ]);
    let z_v = Matrix::new(&[
        [
            crs.v.0 .0.mul(z[(0, 0)].neg()) + crs.v.1 .0.mul(z[(0, 1)].neg()),
            crs.v.0 .1.mul(z[(0, 0)].neg()) + crs.v.1 .1.mul(z[(0, 1)].neg()),
        ],
        [
            crs.v.0 .0.mul(z[(1, 0)].neg()) + crs.v.1 .0.mul(z[(1, 1)].neg()),
            crs.v.0 .1.mul(z[(1, 0)].neg()) + crs.v.1 .1.mul(z[(1, 1)].neg()),
        ],
    ]);
    let (m, n) = equ.gamma.dim();

    let mut t11 = E::ScalarField::zero();
    let mut t12 = E::ScalarField::zero();
    let mut t21 = E::ScalarField::zero();
    let mut t22 = E::ScalarField::zero();

    for (i, x_i) in x.iter().enumerate() {
        for (j, y_j) in y.iter().enumerate() {
            t11 += equ.gamma[(i, j)].mul(x_i.rand.0).mul(y_j.rand.0);
            t12 += equ.gamma[(i, j)].mul(x_i.rand.0).mul(y_j.rand.1);
            t21 += equ.gamma[(i, j)].mul(x_i.rand.1).mul(y_j.rand.0);
            t22 += equ.gamma[(i, j)].mul(x_i.rand.1).mul(y_j.rand.1);
        }
    }

    let phi11 = crs.v.0 .0.mul(t11) + crs.v.1 .0.mul(t12);

    let phi12 = {
        let b_product = equ
            .b
            .iter()
            .zip(y.iter())
            .fold(<E as Pairing>::G2::zero(), |acc, (b_i, y_i)| {
                acc + b_i.mul(y_i.rand.0)
            });
        let y_product = y
            .iter()
            .enumerate()
            .fold(<E as Pairing>::G2::zero(), |acc, (j, y_j)| {
                let mut exp = E::ScalarField::zero();
                for i in 0..m {
                    exp += equ.gamma[(i, j)].mul(y_j.rand.0);
                }
                acc + y_j.value.mul(exp)
            });
        crs.v.0 .1.mul(t11) + crs.v.1 .1.mul(t12) + b_product + y_product
    };

    let phi21 = crs.v.0 .0.mul(t21) + crs.v.1 .0.mul(t22);

    let phi22 = {
        let b_product = equ
            .b
            .iter()
            .zip(y.iter())
            .fold(<E as Pairing>::G2::zero(), |acc, (b_i, y_i)| {
                acc + b_i.mul(y_i.rand.1)
            });
        let y_product = y
            .iter()
            .enumerate()
            .fold(<E as Pairing>::G2::zero(), |acc, (j, y_j)| {
                let mut exp = E::ScalarField::zero();
                for i in 0..m {
                    exp += equ.gamma[(i, j)].mul(y_j.rand.1);
                }
                acc + y_j.value.mul(exp)
            });
        crs.v.0 .1.mul(t21) + crs.v.1 .1.mul(t22) + b_product + y_product
    };

    let phi = Matrix::new(&[[phi11, phi12], [phi21, phi22]]) + z_v;

    let theta11 = <E as Pairing>::G1::zero();

    let theta12 = {
        let a_product = equ
            .a
            .iter()
            .zip(x.iter())
            .fold(<E as Pairing>::G1::zero(), |acc, (a_i, x_i)| {
                acc + a_i.mul(x_i.rand.0)
            });
        let x_product = x
            .iter()
            .enumerate()
            .fold(<E as Pairing>::G1::zero(), |acc, (i, x_i)| {
                let mut exp = E::ScalarField::zero();
                for j in 0..n {
                    exp += equ.gamma[(i, j)].mul(x_i.rand.0);
                }
                acc + x_i.value.mul(exp)
            });
        a_product + x_product
    };

    let theta21 = <E as Pairing>::G1::zero();

    let theta22 = {
        let a_product = equ
            .a
            .iter()
            .zip(x.iter())
            .fold(<E as Pairing>::G1::zero(), |acc, (a_i, x_i)| {
                acc + a_i.mul(x_i.rand.1)
            });
        let x_product = x
            .iter()
            .enumerate()
            .fold(<E as Pairing>::G1::zero(), |acc, (i, x_i)| {
                let mut exp = E::ScalarField::zero();
                for j in 0..n {
                    exp += equ.gamma[(i, j)].mul(x_i.rand.1);
                }
                acc + x_i.value.mul(exp)
            });
        a_product + x_product
    };

    let theta = Matrix::new(&[[theta11, theta12], [theta21, theta22]]) + z_u;

    Proof {
        phi: phi.into(),
        theta: theta.into(),
    }
}
