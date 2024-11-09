//! Defines the struct [Proof] which implements the functions `Prove` and `RdProof` notated in section 6.3
//! in the paper [Fuc10](https://eprint.iacr.org/2010/233.pdf).

use ark_ec::pairing::Pairing;
use ark_std::{rand::Rng, Zero};
use std::ops::{Add, Div, Mul, Neg};

use crate::{
    com::ComRandomness, commit::CommitmentKey, CommitmentKeys, Equation, Matrix, Randomness,
    Variable,
};

/// Contains the components `φ` and `θ` as a Groth-Sahai proof (without internal randomness `Z`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Proof<E: Pairing> {
    pub(crate) phi: Matrix<<E as Pairing>::G2Affine>,
    pub(crate) theta: Matrix<<E as Pairing>::G1Affine>,
}

impl<E: Pairing> Proof<E> {
    /// Implements the `Prove(ck, E, (X, r), (Y, s))` function defined in the paper. Generates a proof `π` = (`φ`, `θ`)
    /// for the equation `E` with the commitment keys `ck` and the variables `X`, `Y` (and their internal
    /// randomness `r`, `s` respectively).
    ///
    /// ## Panics
    /// Panics if 'a.len() != x.len()' or 'b.len() != y.len()', where `a` and `b` are the constants in the equation `E`.
    ///
    /// ## Example
    ///
    /// ```
    /// use ark_bls12_381::Bls12_381 as E;
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::{test_rng, UniformRand};
    /// use gs_ppe::{setup, CommitmentKeys, Equation, Matrix, ProofSystem, Variable};
    ///
    /// type G1 = <E as Pairing>::G1;
    /// type G2 = <E as Pairing>::G2;
    /// type G1Affine = <E as Pairing>::G1Affine;
    /// type G2Affine = <E as Pairing>::G2Affine;
    /// type Fr = <E as Pairing>::ScalarField;
    ///
    /// let rng = &mut test_rng();
    /// let (a, b) = (G1Affine::rand(rng), G2Affine::rand(rng));
    /// let (x_value, y_value) = (G1Affine::rand(rng), G2Affine::rand(rng));
    /// let (x, y) = (
    ///     Variable::<G1>::new(rng, x_value),
    ///     Variable::<G2>::new(rng, y_value),
    /// );
    /// let gamma = Matrix::<Fr>::rand(rng, 1, 1);
    ///
    /// let cks = CommitmentKeys::<E>::rand(rng);
    ///
    /// // Setup Proof System over Pairing Product Equation:
    /// // e(a, y) + e(x, b) + e(x, y)^gamma = T
    /// let ProofSystem {
    ///     equation,
    ///     c,
    ///     d,
    ///     proof,
    /// } = setup(rng, &cks, &[(a, y)], &[(x, b)], &gamma);
    ///
    /// assert!(equation.verify(&cks, &c, &d, &proof));
    /// ```
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

        let z_u = z_u::<E>(&z, &cks.u);
        let z_v = z_v::<E>(&z, &cks.v);

        let (t11, t12, t21, t22) = t11_t12_t21_t22::<E>(
            &x.iter().map(|x_i| x_i.rand).collect::<Vec<_>>(),
            &y.iter().map(|y_j| y_j.rand).collect::<Vec<_>>(),
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

        // Compute φ as in (7).
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

        // Compute θ as in (7).
        let theta = Matrix::new(&[[theta11, theta12], [theta21, theta22]]) + z_u;

        // π = (φ, θ)
        Proof {
            phi: phi.into(),
            theta: theta.into(),
        }
    }

    /// Implements the Proof Randomization function (proof adaption) `RdProof(ck, E, (c, r), (d, s)), π)` defined in the paper.
    /// Randomized the components (`φ`, `θ`) in this proof for the equation `E` with the commitment keys `ck` and the Commitments `c`, `d`
    /// (and their internal randomness `r`, `s` respectively).
    ///
    /// ## Panics
    /// Panics if 'a.len() != ds.len()' or 'b.len() != cr.len()', where `a` and `b` are the constants in the equation `E`.
    ///
    /// ## Example
    ///
    /// ```
    /// use ark_bls12_381::Bls12_381 as E;
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::{test_rng, UniformRand};
    /// use gs_ppe::{setup, CommitmentKeys, Equation, Matrix, ProofSystem, Variable};
    ///
    /// type G1 = <E as Pairing>::G1;
    /// type G2 = <E as Pairing>::G2;
    /// type G1Affine = <E as Pairing>::G1Affine;
    /// type G2Affine = <E as Pairing>::G2Affine;
    /// type Fr = <E as Pairing>::ScalarField;
    ///
    /// let rng = &mut test_rng();
    /// let (a, b) = (G1Affine::rand(rng), G2Affine::rand(rng));
    /// let (x_value, y_value) = (G1Affine::rand(rng), G2Affine::rand(rng));
    /// let (x, y) = (
    ///     Variable::<G1>::new(rng, x_value),
    ///     Variable::<G2>::new(rng, y_value),
    /// );
    /// let gamma = Matrix::<Fr>::rand(rng, 1, 1);
    ///
    /// let cks = CommitmentKeys::<E>::rand(rng);
    ///
    /// // Setup Proof System over Pairing Product Equation:
    /// // e(a, y) + e(x, b) + e(x, y)^gamma = T
    /// let proof_system = setup(rng, &cks, &[(a, y)], &[(x, b)], &gamma);
    /// let ProofSystem {
    ///     equation,
    ///     c,
    ///     d,
    ///     proof,
    /// } = proof_system.randomize(rng, &cks);
    ///
    /// assert!(equation.verify(&cks, &c, &d, &proof));
    /// ```
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
        let z_u = z_u::<E>(&z, &cks.u);
        let z_v = z_v::<E>(&z, &cks.v);

        let c = cr.iter().map(|(c_i, _)| *c_i).collect::<Vec<_>>();
        let d = ds.iter().map(|(d_j, _)| *d_j).collect::<Vec<_>>();
        let r = cr.iter().map(|(_, r_i)| *r_i).collect::<Vec<_>>();
        let s = ds.iter().map(|(_, s_j)| *s_j).collect::<Vec<_>>();

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

impl<E: Pairing> Add for Proof<E> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let phi = self.phi.into::<<E as Pairing>::G2>() + other.phi.into::<<E as Pairing>::G2>();
        let theta =
            self.theta.into::<<E as Pairing>::G1>() + other.theta.into::<<E as Pairing>::G1>();
        Proof {
            phi: phi.into(),
            theta: theta.into(),
        }
    }
}

impl<E: Pairing> Div for Proof<E> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        let phi =
            self.phi.into::<<E as Pairing>::G2>() + other.phi.into::<<E as Pairing>::G2>().neg();
        let theta = self.theta.into::<<E as Pairing>::G1>()
            + other.theta.into::<<E as Pairing>::G1>().neg();
        Proof {
            phi: phi.into(),
            theta: theta.into(),
        }
    }
}

/// Computes the matrix `Z (x) u` defined in (5).
fn z_u<E: Pairing>(
    z: &Matrix<E::ScalarField>,
    u: &CommitmentKey<<E as Pairing>::G1>,
) -> Matrix<<E as Pairing>::G1> {
    Matrix::new(&[
        [
            u.0 .0.mul(z[(0, 0)]) + u.1 .0.mul(z[(0, 1)]),
            u.0 .1.mul(z[(0, 0)]) + u.1 .1.mul(z[(0, 1)]),
        ],
        [
            u.0 .0.mul(z[(1, 0)]) + u.1 .0.mul(z[(1, 1)]),
            u.0 .1.mul(z[(1, 0)]) + u.1 .1.mul(z[(1, 1)]),
        ],
    ])
}

/// Computes the matrix `Z (x) v` defined in (5).
fn z_v<E: Pairing>(
    z: &Matrix<E::ScalarField>,
    v: &CommitmentKey<<E as Pairing>::G2>,
) -> Matrix<<E as Pairing>::G2> {
    Matrix::new(&[
        [
            v.0 .0.mul(z[(0, 0)].neg()) + v.1 .0.mul(z[(1, 0)].neg()),
            v.0 .1.mul(z[(0, 0)].neg()) + v.1 .1.mul(z[(1, 0)].neg()),
        ],
        [
            v.0 .0.mul(z[(0, 1)].neg()) + v.1 .0.mul(z[(1, 1)].neg()),
            v.0 .1.mul(z[(0, 1)].neg()) + v.1 .1.mul(z[(1, 1)].neg()),
        ],
    ])
}

/// Computes the values defined in (6).
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
