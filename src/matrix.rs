use std::ops::{Add, Index};

use ark_ff::{UniformRand, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Valid};
use ark_std::rand::Rng;
use ndarray::{Array, Ix2};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix<F>
where
    F: Clone,
{
    inner: Array<F, Ix2>,
}

impl<F> Matrix<F>
where
    F: Clone,
{
    pub fn new<const N: usize>(xs: &[[F; N]]) -> Self {
        Self {
            inner: ndarray::arr2(xs),
        }
    }

    /// An array with zero rows and n columns of zeros.
    pub fn zeros_column(n: usize) -> Self
    where
        F: Zero,
    {
        Self {
            inner: Array::zeros((0, n)),
        }
    }

    pub fn rand<R: Rng>(rng: &mut R, rows: usize, cols: usize) -> Self
    where
        F: UniformRand,
    {
        Self {
            inner: Array::from_shape_fn((rows, cols), |_| F::rand(rng)),
        }
    }

    pub fn to_vecs(&self) -> Vec<Vec<F>> {
        self.inner
            .outer_iter()
            .map(|row| row.iter().cloned().collect())
            .collect()
    }

    pub fn from_vecs(vecs: Vec<Vec<F>>) -> Self {
        Self {
            inner: Array::from_shape_vec(
                (vecs.len(), vecs[0].len()),
                vecs.into_iter().flatten().collect(),
            )
            .unwrap(),
        }
    }

    pub fn into<G>(self) -> Matrix<G>
    where
        G: Clone + From<F>,
    {
        Matrix {
            inner: self.inner.mapv(G::from),
        }
    }

    #[inline]
    pub fn dim(&self) -> (usize, usize) {
        self.inner.dim()
    }
}

impl<F, G> From<Array<G, Ix2>> for Matrix<F>
where
    G: Clone,
    F: Clone + From<G>,
{
    fn from(inner: Array<G, Ix2>) -> Self {
        Self {
            inner: inner.mapv(F::from),
        }
    }
}

impl<F> AsRef<Array<F, Ix2>> for Matrix<F>
where
    F: Clone,
{
    fn as_ref(&self) -> &Array<F, Ix2> {
        &self.inner
    }
}

impl<F> Index<(usize, usize)> for Matrix<F>
where
    F: Clone,
{
    type Output = F;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.inner[index]
    }
}

impl<F> Valid for Matrix<F>
where
    F: Clone + Valid,
{
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        Ok(())
    }
}

impl<F> CanonicalSerialize for Matrix<F>
where
    F: Clone + CanonicalSerialize,
{
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        Vec::<Vec<F>>::serialize_with_mode(&self.to_vecs(), writer, compress)
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        Vec::<Vec<F>>::serialized_size(&self.to_vecs(), compress)
    }
}

impl<F> CanonicalDeserialize for Matrix<F>
where
    F: Clone + CanonicalDeserialize,
{
    fn deserialize_with_mode<R: ark_serialize::Read>(
        reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        Vec::<Vec<F>>::deserialize_with_mode(reader, compress, validate).map(Self::from_vecs)
    }
}

impl<F> Add for Matrix<F>
where
    F: Clone + Add<Output = F>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            inner: self.inner + rhs.inner,
        }
    }
}