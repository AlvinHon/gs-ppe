use ark_ec::CurveGroup;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

#[derive(Copy, Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct Com<G: CurveGroup>(pub G::Affine, pub G::Affine);
