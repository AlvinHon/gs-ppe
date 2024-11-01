//! Defines the struct [ExtractKey], the key `ek` for extracting `SXDH Commitments`` defined in section 6.2 in
//! the paper [Fuc10](https://eprint.iacr.org/2010/233.pdf).

use ark_ec::pairing::Pairing;
use std::ops::{Mul, Neg};

use crate::com::Com;

/// The key `ek` for extracting `SXDH Commitments`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct ExtractKey<E: Pairing>(pub E::ScalarField, pub E::ScalarField);

impl<E: Pairing> ExtractKey<E> {
    /// The extract function `Extr(ek, c)` defined in the paper. Extracts the commitment
    /// and returns the value of the committed variable in Group 1.
    ///
    /// ## Example
    ///
    /// ```
    /// use ark_bls12_381::Bls12_381 as E;
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::{test_rng, UniformRand};
    /// use gs_ppe::{CommitmentKeys, Variable};
    ///
    /// type G1Affine = <E as Pairing>::G1Affine;
    ///
    /// let rng = &mut test_rng();
    ///
    /// let (cks, ek) = CommitmentKeys::<E>::rand_ex(rng);
    ///
    /// let x_value = G1Affine::rand(rng);
    /// let x = Variable::new(rng, x_value);
    ///
    /// let c = cks.u.commit(&x);
    /// let x_prime = ek.extract_1(&c);
    ///
    /// assert_eq!(x_prime, x_value);
    /// ```
    pub fn extract_1(&self, c: &Com<<E as Pairing>::G1>) -> E::G1Affine {
        (c.0.mul(&self.0.neg()) + c.1).into()
    }
    /// The extract function `Extr(ek, c)` defined in the paper. Extracts the commitment
    /// and returns the value of the committed variable in Group 2.
    ///
    /// ## Example
    ///
    /// ```rust
    /// use ark_bls12_381::Bls12_381 as E;
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::{test_rng, UniformRand};
    /// use gs_ppe::{CommitmentKeys, Variable};
    ///
    /// type G2Affine = <E as Pairing>::G2Affine;
    ///
    /// let rng = &mut test_rng();
    ///
    /// let (cks, ek) = CommitmentKeys::<E>::rand_ex(rng);
    ///
    /// let y_value = G2Affine::rand(rng);
    /// let y = Variable::new(rng, y_value);
    ///
    /// let d = cks.v.commit(&y);
    /// let y_prime = ek.extract_2(&d);
    ///
    /// assert_eq!(y_prime, y_value);
    /// ```
    pub fn extract_2(&self, c: &Com<<E as Pairing>::G2>) -> E::G2Affine {
        (c.0.mul(&self.1.neg()) + c.1).into()
    }
}
