//! Defines the struct [Equation], the equation `E` notated in section 6.3 in the paper [Fuc10](https://eprint.iacr.org/2010/233.pdf).

use ark_ec::pairing::{Pairing, PairingOutput};
use ark_std::Zero;
use std::ops::{Add, Mul};

use crate::{Com, CommitmentKeys, Matrix, Proof};

/// The pairing product equation `E`, represented by:
/// - the constant `a` in a vector of size `n`
/// - the constant `b` in a vector of size `m`
/// - the constant `gamma` in a matrix of dimension `(m, n)`
/// - the target value `target` in the pairing output
///
/// where the values satisfies equation:
///
/// Π e(a_i, y_i) Π e(x_i, b_i) ΠΠ e(x_i, y_i)^gamma_ij = target
///
/// for some x and y.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Equation<E: Pairing> {
    pub(crate) a: Vec<<E as Pairing>::G1>,    // size = n
    pub(crate) b: Vec<<E as Pairing>::G2>,    // size = m
    pub(crate) gamma: Matrix<E::ScalarField>, // dim = (m, n)
    pub(crate) target: PairingOutput<E>,
}

impl<E: Pairing> Equation<E> {
    /// Constructs an equation `E` with the given constants `a`, `b`, `gamma`, and `target`.
    ///
    /// ## Panics
    /// Panics if the dimension of `gamma` != (m, n), where m = b.len() and n = a.len().
    pub fn new(
        a: Vec<<E as Pairing>::G1Affine>,
        b: Vec<<E as Pairing>::G2Affine>,
        gamma: Matrix<E::ScalarField>,
        target: PairingOutput<E>,
    ) -> Self {
        assert_eq!(gamma.dim(), (b.len(), a.len()));
        Self {
            a: a.into_iter().map(|a| a.into()).collect(),
            b: b.into_iter().map(|b| b.into()).collect(),
            gamma,
            target,
        }
    }

    // TODO:
    // "Remark 5. Blazy et al. [BFI+10] show that by using techniques of batch verification, the number of pairing
    // computations can be reduced from 4m + n + 16 to 2m+n+8".

    /// The Verification function `Verify(ck, E, c, d, (φ, θ))`. Verifies the equation `E` with
    /// the given commitments `c`, `d`, and `proof`. Returns false if the verification fails or
    /// the dimensions of the inputs are incorrect.
    ///
    /// ## Example
    ///
    /// ```
    /// use ark_bls12_381::Bls12_381 as E;
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::{test_rng, UniformRand};
    /// use std::ops::Mul;
    /// use gs_ppe::{setup, CommitmentKeys, Equation, Matrix, Variable};
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
    ///
    /// // Explicitly create the equation:
    /// // e(a, y) + e(x, b) + e(x, y)^gamma = T
    /// let target = E::pairing(a, y_value)
    ///     + E::pairing(x_value, b)
    ///     + E::pairing(x_value, y_value).mul(gamma[(0, 0)]);
    /// let equation = Equation::<E>::new(vec![a], vec![b], gamma, target);
    ///
    /// assert_eq!(equation, proof_system.equation);
    /// assert!(equation.verify(&cks, &proof_system.c, &proof_system.d, &proof_system.proof));
    /// ```
    pub fn verify(
        &self,
        cks: &CommitmentKeys<E>,
        c: &[Com<<E as Pairing>::G1>],
        d: &[Com<<E as Pairing>::G2>],
        proof: &Proof<E>,
    ) -> bool {
        let (m, n) = self.gamma.dim();
        if self.a.len() != n
            || self.b.len() != m
            || c.len() != m
            || d.len() != n
            || proof.phi.dim() != (2, 2)
            || proof.theta.dim() != (2, 2)
        {
            return false;
        }
        let u = &cks.u;
        let v = &cks.v;

        // Check Equation 1:
        // Π e(c_i1, Π d_j1^gamma_ij) = e(u11, φ11) e(u21, φ21) e(θ11, v11) e(θ21, v21)
        let lhs = c
            .iter()
            .enumerate()
            .fold(PairingOutput::zero(), |acc, (i, c_i)| {
                let d_product = d
                    .iter()
                    .enumerate()
                    .fold(<E as Pairing>::G2::zero(), |acc, (j, d_j)| {
                        acc + d_j.0.mul(self.gamma[(i, j)])
                    });

                acc + E::pairing(c_i.0, d_product)
            });
        let rhs = E::pairing(u.0 .0, proof.phi[(0, 0)])
            + E::pairing(u.1 .0, proof.phi[(1, 0)])
            + E::pairing(proof.theta[(0, 0)], v.0 .0)
            + E::pairing(proof.theta[(1, 0)], v.1 .0);

        if lhs != rhs {
            return false;
        }

        // create pre-calculated value b_i Π d_j2^gamma_ij for equation 2 and 4 for efficiency.
        let b_d = c.iter().enumerate().fold(Vec::new(), |mut acc, (i, _)| {
            let d_product = d
                .iter()
                .enumerate()
                .fold(<E as Pairing>::G2::zero(), |acc, (j, d_j)| {
                    acc + d_j.1.mul(self.gamma[(i, j)])
                });
            acc.push(self.b[i] + d_product);
            acc
        });

        // Check Equation 2:
        // Π e(c_i1, b_i Π d_j2^gamma_ij) = e(u11, φ12) e(u21, φ22) e(θ11, v12) e(θ21, v22)
        let lhs = c
            .iter()
            .enumerate()
            .fold(PairingOutput::zero(), |acc, (i, c_i)| {
                acc + E::pairing(c_i.0, b_d[i])
            });
        let rhs = E::pairing(u.0 .0, proof.phi[(0, 1)])
            + E::pairing(u.1 .0, proof.phi[(1, 1)])
            + E::pairing(proof.theta[(0, 0)], v.0 .1)
            + E::pairing(proof.theta[(1, 0)], v.1 .1);
        if lhs != rhs {
            return false;
        }

        // Check Equation 3:
        // Π e(a_j Π c_i2^gamma_ij, d_j1) = e(u12, φ11) e(u22, φ21) e(θ12, v11) e(θ22, v21)
        let lhs = d
            .iter()
            .enumerate()
            .fold(PairingOutput::zero(), |acc, (j, d_j)| {
                let c_product = c
                    .iter()
                    .enumerate()
                    .fold(<E as Pairing>::G1::zero(), |acc, (i, c_i)| {
                        acc + c_i.1.mul(self.gamma[(i, j)])
                    });

                acc + E::pairing(self.a[j] + c_product, d_j.0)
            });
        let rhs = E::pairing(u.0 .1, proof.phi[(0, 0)])
            + E::pairing(u.1 .1, proof.phi[(1, 0)])
            + E::pairing(proof.theta[(0, 1)], v.0 .0)
            + E::pairing(proof.theta[(1, 1)], v.1 .0);
        if lhs != rhs {
            return false;
        }

        // Check Equation 4:
        // Π e(a_j, d_j2) Π e(c_i2, b_i Π d_j2^gamma_ij) = t_T e(u12, φ12) e(u22, φ22) e(θ12, v12) e(θ22, v22)
        let lhs = {
            let a_d = self
                .a
                .iter()
                .zip(d.iter())
                .fold(PairingOutput::zero(), |acc, (a_j, d_j)| {
                    acc + E::pairing(a_j, d_j.1)
                });
            let c_bd = c
                .iter()
                .enumerate()
                .fold(PairingOutput::zero(), |acc, (i, c_i)| {
                    acc + E::pairing(c_i.1, b_d[i])
                });
            a_d + c_bd
        };
        let rhs = self.target
            + E::pairing(u.0 .1, proof.phi[(0, 1)])
            + E::pairing(u.1 .1, proof.phi[(1, 1)])
            + E::pairing(proof.theta[(0, 1)], v.0 .1)
            + E::pairing(proof.theta[(1, 1)], v.1 .1);

        lhs == rhs
    }
}

impl<E: Pairing> Add for Equation<E> {
    type Output = Self;

    fn add(self, mut rhs: Self) -> Self {
        let Equation {
            mut a,
            mut b,
            gamma,
            mut target,
        } = self;

        a.append(&mut rhs.a);
        b.append(&mut rhs.b);

        target += rhs.target;

        // Compute [[ gamma1, 0], [0, gamma2]]
        let gamma = {
            let (m, n) = gamma.dim();
            let (m_prime, n_prime) = rhs.gamma.dim();
            let mut gamma1 = gamma.take();
            let gamma22 = rhs.gamma.take();

            let zeros12 = ndarray::Array2::from_elem((m_prime, n_prime), E::ScalarField::zero());
            let mut gamma2 = ndarray::Array2::from_elem((m, n), E::ScalarField::zero());

            gamma1.append(ndarray::Axis(1), zeros12.view()).unwrap();
            gamma2.append(ndarray::Axis(1), gamma22.view()).unwrap();

            gamma1.append(ndarray::Axis(0), gamma2.view()).unwrap();

            Matrix::from(gamma1)
        };

        Self {
            a,
            b,
            gamma,
            target,
        }
    }
}
