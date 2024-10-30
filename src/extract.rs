use ark_ec::pairing::Pairing;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use std::ops::{Mul, Neg};

use crate::com::Com;

#[derive(Copy, Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct ExtractKey<E: Pairing>(pub E::ScalarField, pub E::ScalarField);

impl<E: Pairing> ExtractKey<E> {
    pub fn extract_1(&self, c: &Com<<E as Pairing>::G1>) -> E::G1Affine {
        (c.0.mul(&self.0.neg()) + c.1).into()
    }

    pub fn extract_2(&self, c: &Com<<E as Pairing>::G2>) -> E::G2Affine {
        (c.0.mul(&self.1.neg()) + c.1).into()
    }
}

#[cfg(test)]
mod test {
    use ark_bls12_381::Bls12_381 as F;
    use ark_ec::pairing::Pairing;
    use ark_std::{test_rng, UniformRand};

    type G1Affine = <F as Pairing>::G1Affine;
    type G2Affine = <F as Pairing>::G2Affine;

    use crate::{CommitmentKeys, Variable};

    #[test]
    fn test_extract() {
        let rng = &mut test_rng();
        let x_value = G1Affine::rand(rng);
        let y_value = G2Affine::rand(rng);

        let x = Variable::new(rng, x_value);
        let y = Variable::new(rng, y_value);

        let (cks, ek) = CommitmentKeys::<F>::rand_ex(rng);

        let c = cks.u.commit(&x);
        let c_prime = ek.extract_1(&c);
        let d = cks.v.commit(&y);
        let d_prime = ek.extract_2(&d);

        assert_eq!(c_prime, x_value);
        assert_eq!(d_prime, y_value);
    }
}
