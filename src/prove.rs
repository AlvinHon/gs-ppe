use ark_ec::pairing::Pairing;
use ark_ff::Zero;
use ark_std::rand::Rng;
use std::ops::{Mul, Neg};

use crate::{
    com::ComRandomness, commit::CommitmentKeys, equation::Equation, matrix::Matrix,
    variable::Variable, Randomness,
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Proof<E: Pairing> {
    pub(crate) phi: Matrix<<E as Pairing>::G2Affine>,
    pub(crate) theta: Matrix<<E as Pairing>::G1Affine>,
}

impl<E: Pairing> Proof<E> {
    pub fn new<R: Rng>(
        rng: &mut R,
        cks: &CommitmentKeys<E>,
        equ: &Equation<E>,
        x: &[Variable<<E as Pairing>::G1>],
        y: &[Variable<<E as Pairing>::G2>],
    ) -> Self {
        assert_eq!(equ.a.len(), y.len());
        assert_eq!(equ.b.len(), x.len());

        let z = Matrix::<E::ScalarField>::rand(rng, 2, 2);

        let z_u = z_u(cks, &z);
        let z_v = z_v(cks, &z);

        let (t11, t12, t21, t22) = t11_t12_t21_t22::<E>(
            &x.iter().map(|x_i| x_i.rand.clone()).collect::<Vec<_>>(),
            &y.iter().map(|y_j| y_j.rand.clone()).collect::<Vec<_>>(),
            &equ.gamma,
        );

        let phi11 = cks.v.0 .0.mul(t11) + cks.v.1 .0.mul(t12);

        let phi12 = {
            let b_product = equ
                .b
                .iter()
                .zip(x.iter())
                .fold(<E as Pairing>::G2::zero(), |acc, (b_i, x_i)| {
                    acc + b_i.mul(x_i.rand.0)
                });
            let y_product =
                y.iter()
                    .enumerate()
                    .fold(<E as Pairing>::G2::zero(), |acc, (j, y_j)| {
                        let exp = x
                            .iter()
                            .enumerate()
                            .fold(E::ScalarField::zero(), |acc, (i, x_i)| {
                                acc + equ.gamma[(i, j)].mul(x_i.rand.0)
                            });
                        acc + y_j.value.mul(exp)
                    });
            cks.v.0 .1.mul(t11) + cks.v.1 .1.mul(t12) + b_product + y_product
        };

        let phi21 = cks.v.0 .0.mul(t21) + cks.v.1 .0.mul(t22);

        let phi22 = {
            let b_product = equ
                .b
                .iter()
                .zip(x.iter())
                .fold(<E as Pairing>::G2::zero(), |acc, (b_i, x_i)| {
                    acc + b_i.mul(x_i.rand.1)
                });
            let y_product =
                y.iter()
                    .enumerate()
                    .fold(<E as Pairing>::G2::zero(), |acc, (j, y_j)| {
                        let exp = x
                            .iter()
                            .enumerate()
                            .fold(E::ScalarField::zero(), |acc, (i, x_i)| {
                                acc + equ.gamma[(i, j)].mul(x_i.rand.1)
                            });
                        acc + y_j.value.mul(exp)
                    });
            cks.v.0 .1.mul(t21) + cks.v.1 .1.mul(t22) + b_product + y_product
        };

        let phi = Matrix::new(&[[phi11, phi12], [phi21, phi22]]) + z_v;

        let theta11 = <E as Pairing>::G1::zero();

        let theta12 = {
            let a_product = equ
                .a
                .iter()
                .zip(y.iter())
                .fold(<E as Pairing>::G1::zero(), |acc, (a_j, y_j)| {
                    acc + a_j.mul(y_j.rand.0)
                });
            let x_product =
                x.iter()
                    .enumerate()
                    .fold(<E as Pairing>::G1::zero(), |acc, (i, x_i)| {
                        let exp = y
                            .iter()
                            .enumerate()
                            .fold(E::ScalarField::zero(), |acc, (j, y_j)| {
                                acc + equ.gamma[(i, j)].mul(y_j.rand.0)
                            });
                        acc + x_i.value.mul(exp)
                    });
            a_product + x_product
        };

        let theta21 = <E as Pairing>::G1::zero();

        let theta22 = {
            let a_product = equ
                .a
                .iter()
                .zip(y.iter())
                .fold(<E as Pairing>::G1::zero(), |acc, (a_j, y_j)| {
                    acc + a_j.mul(y_j.rand.1)
                });
            let x_product =
                x.iter()
                    .enumerate()
                    .fold(<E as Pairing>::G1::zero(), |acc, (i, x_i)| {
                        let exp = y
                            .iter()
                            .enumerate()
                            .fold(E::ScalarField::zero(), |acc, (j, y_j)| {
                                acc + equ.gamma[(i, j)].mul(y_j.rand.1)
                            });
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

    pub fn randomize<R: Rng>(
        &mut self,
        rng: &mut R,
        cks: &CommitmentKeys<E>,
        equ: &Equation<E>,
        cr: &[ComRandomness<<E as Pairing>::G1>],
        ds: &[ComRandomness<<E as Pairing>::G2>],
    ) {
        assert_eq!(equ.a.len(), ds.len());
        assert_eq!(equ.b.len(), cr.len());

        let z = Matrix::<E::ScalarField>::rand(rng, 2, 2);
        let z_u = z_u(cks, &z);
        let z_v = z_v(cks, &z);

        let c = cr.iter().map(|(c_i, _)| *c_i).collect::<Vec<_>>();
        let d = ds.iter().map(|(d_j, _)| *d_j).collect::<Vec<_>>();
        let r = cr.iter().map(|(_, r_i)| r_i.clone()).collect::<Vec<_>>();
        let s = ds.iter().map(|(_, s_j)| s_j.clone()).collect::<Vec<_>>();

        let (t11, t12, t21, t22) = t11_t12_t21_t22::<E>(&r, &s, &equ.gamma);

        self.phi = {
            let phi11 = {
                let d_product =
                    d.iter()
                        .enumerate()
                        .fold(<E as Pairing>::G2::zero(), |acc, (j, d_j)| {
                            let exp = r
                                .iter()
                                .enumerate()
                                .fold(E::ScalarField::zero(), |acc, (i, r_i)| {
                                    acc + equ.gamma[(i, j)].mul(r_i.0)
                                });
                            acc + d_j.0.mul(exp)
                        });
                let v_product = cks.v.0 .0.mul(&t11) + cks.v.1 .0.mul(&t12);
                d_product + v_product
            };

            let phi12 = {
                let b_product = equ
                    .b
                    .iter()
                    .zip(r.iter())
                    .fold(<E as Pairing>::G2::zero(), |acc, (b_i, r_i)| {
                        acc + b_i.mul(r_i.0)
                    });
                let d_product =
                    d.iter()
                        .enumerate()
                        .fold(<E as Pairing>::G2::zero(), |acc, (j, d_j)| {
                            let exp = r
                                .iter()
                                .enumerate()
                                .fold(E::ScalarField::zero(), |acc, (i, r_i)| {
                                    acc + equ.gamma[(i, j)].mul(r_i.0)
                                });
                            acc + d_j.1.mul(exp)
                        });
                let v_product = cks.v.0 .1.mul(&t11) + cks.v.1 .1.mul(&t12);
                b_product + d_product + v_product
            };

            let phi21 = {
                let d_product =
                    d.iter()
                        .enumerate()
                        .fold(<E as Pairing>::G2::zero(), |acc, (j, d_j)| {
                            let exp = r
                                .iter()
                                .enumerate()
                                .fold(E::ScalarField::zero(), |acc, (i, r_i)| {
                                    acc + equ.gamma[(i, j)].mul(r_i.1)
                                });
                            acc + d_j.0.mul(exp)
                        });
                let v_product = cks.v.0 .0.mul(&t21) + cks.v.1 .0.mul(&t22);
                d_product + v_product
            };

            let phi22 = {
                let b_product = equ
                    .b
                    .iter()
                    .zip(r.iter())
                    .fold(<E as Pairing>::G2::zero(), |acc, (b_i, r_i)| {
                        acc + b_i.mul(r_i.1)
                    });
                let d_product =
                    d.iter()
                        .enumerate()
                        .fold(<E as Pairing>::G2::zero(), |acc, (j, d_j)| {
                            let exp = r
                                .iter()
                                .enumerate()
                                .fold(E::ScalarField::zero(), |acc, (i, r_i)| {
                                    acc + equ.gamma[(i, j)].mul(r_i.1)
                                });
                            acc + d_j.1.mul(exp)
                        });
                let v_product = cks.v.0 .1.mul(&t21) + cks.v.1 .1.mul(&t22);
                b_product + d_product + v_product
            };

            (self.phi.clone().into::<<E as Pairing>::G2>()
                + Matrix::new(&[[phi11, phi12], [phi21, phi22]])
                + z_v)
                .into()
        };

        self.theta = {
            let theta11 = c
                .iter()
                .enumerate()
                .fold(<E as Pairing>::G1::zero(), |acc, (i, c_i)| {
                    let exp = s
                        .iter()
                        .enumerate()
                        .fold(E::ScalarField::zero(), |acc, (j, s_j)| {
                            acc + equ.gamma[(i, j)].mul(s_j.0)
                        });
                    acc + c_i.0.mul(exp)
                });

            let theta12 = {
                let a_product = equ
                    .a
                    .iter()
                    .zip(s.iter())
                    .fold(<E as Pairing>::G1::zero(), |acc, (a_j, s_j)| {
                        acc + a_j.mul(s_j.0)
                    });
                let c_product =
                    c.iter()
                        .enumerate()
                        .fold(<E as Pairing>::G1::zero(), |acc, (i, c_i)| {
                            let exp = s
                                .iter()
                                .enumerate()
                                .fold(E::ScalarField::zero(), |acc, (j, s_j)| {
                                    acc + equ.gamma[(i, j)].mul(s_j.0)
                                });
                            acc + c_i.1.mul(exp)
                        });
                a_product + c_product
            };

            let theta21 = c
                .iter()
                .enumerate()
                .fold(<E as Pairing>::G1::zero(), |acc, (i, c_i)| {
                    let exp = s
                        .iter()
                        .enumerate()
                        .fold(E::ScalarField::zero(), |acc, (j, s_j)| {
                            acc + equ.gamma[(i, j)].mul(s_j.1)
                        });
                    acc + c_i.0.mul(exp)
                });

            let theta22 = {
                let a_product = equ
                    .a
                    .iter()
                    .zip(s.iter())
                    .fold(<E as Pairing>::G1::zero(), |acc, (a_j, s_j)| {
                        acc + a_j.mul(s_j.1)
                    });
                let c_product =
                    c.iter()
                        .enumerate()
                        .fold(<E as Pairing>::G1::zero(), |acc, (i, c_i)| {
                            let exp = s
                                .iter()
                                .enumerate()
                                .fold(E::ScalarField::zero(), |acc, (j, s_j)| {
                                    acc + equ.gamma[(i, j)].mul(s_j.1)
                                });
                            acc + c_i.1.mul(exp)
                        });
                a_product + c_product
            };

            (self.theta.clone().into::<<E as Pairing>::G1>()
                + Matrix::new(&[[theta11, theta12], [theta21, theta22]])
                + z_u)
                .into()
        };
    }
}

fn z_u<E: Pairing>(
    cks: &CommitmentKeys<E>,
    z: &Matrix<E::ScalarField>,
) -> Matrix<<E as Pairing>::G1> {
    Matrix::new(&[
        [
            cks.u.0 .0.mul(z[(0, 0)]) + cks.u.1 .0.mul(z[(0, 1)]),
            cks.u.0 .1.mul(z[(0, 0)]) + cks.u.1 .1.mul(z[(0, 1)]),
        ],
        [
            cks.u.0 .0.mul(z[(1, 0)]) + cks.u.1 .0.mul(z[(1, 1)]),
            cks.u.0 .1.mul(z[(1, 0)]) + cks.u.1 .1.mul(z[(1, 1)]),
        ],
    ])
}

fn z_v<E: Pairing>(
    cks: &CommitmentKeys<E>,
    z: &Matrix<E::ScalarField>,
) -> Matrix<<E as Pairing>::G2> {
    Matrix::new(&[
        [
            cks.v.0 .0.mul(z[(0, 0)].neg()) + cks.v.1 .0.mul(z[(1, 0)].neg()),
            cks.v.0 .1.mul(z[(0, 0)].neg()) + cks.v.1 .1.mul(z[(1, 0)].neg()),
        ],
        [
            cks.v.0 .0.mul(z[(0, 1)].neg()) + cks.v.1 .0.mul(z[(1, 1)].neg()),
            cks.v.0 .1.mul(z[(0, 1)].neg()) + cks.v.1 .1.mul(z[(1, 1)].neg()),
        ],
    ])
}

fn t11_t12_t21_t22<E: Pairing>(
    r: &[Randomness<<E as Pairing>::G1>],
    s: &[Randomness<<E as Pairing>::G2>],
    gamma: &Matrix<E::ScalarField>,
) -> (
    E::ScalarField,
    E::ScalarField,
    E::ScalarField,
    E::ScalarField,
) {
    let mut t11 = E::ScalarField::zero();
    let mut t12 = E::ScalarField::zero();
    let mut t21 = E::ScalarField::zero();
    let mut t22 = E::ScalarField::zero();

    for i in 0..r.len() {
        for j in 0..s.len() {
            let r_i = &r[i];
            let s_j = &s[j];
            t11 += gamma[(i, j)].mul(r_i.0).mul(s_j.0);
            t12 += gamma[(i, j)].mul(r_i.0).mul(s_j.1);
            t21 += gamma[(i, j)].mul(r_i.1).mul(s_j.0);
            t22 += gamma[(i, j)].mul(r_i.1).mul(s_j.1);
        }
    }
    (t11, t12, t21, t22)
}
