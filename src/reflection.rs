//! Geometric (reflection) representation of Coxeter groups.
//!
//! Represents group elements as matrices acting on a real vector space.

use crate::matrix::CoxeterMatrix;

/// Reflection representation of a Coxeter group.
///
/// Each generator s_i acts as a reflection in the vector space R^n.
#[derive(Debug, Clone)]
pub struct ReflectionRepresentation {
    /// The Coxeter matrix defining the group.
    pub coxeter: CoxeterMatrix,
    /// Reflection matrices for each generator.
    pub reflections: Vec<Vec<Vec<f64>>>,
}

impl ReflectionRepresentation {
    /// Create the canonical reflection representation for the given Coxeter matrix.
    ///
    /// Uses the Tits representation: the bilinear form B with
    /// B(e_i, e_j) = -cos(π/m_ij) (or -1 if m_ij = ∞).
    pub fn new(coxeter: &CoxeterMatrix) -> Self {
        let n = coxeter.rank;
        let reflections = (0..n).map(|i| compute_reflection(coxeter, i)).collect();
        Self {
            coxeter: coxeter.clone(),
            reflections,
        }
    }

    /// Get the bilinear form matrix B.
    pub fn bilinear_form(&self) -> Vec<Vec<f64>> {
        let n = self.coxeter.rank;
        let mut b = vec![vec![0.0; n]; n];
        for i in 0..n {
            b[i][i] = 1.0;
            for j in (i + 1)..n {
                let m = self.coxeter.entries[i][j];
                let val = if m == 0 {
                    1.0 // corresponds to cos(0) = 1 for infinite order
                } else if m == 2 {
                    0.0 // cos(π/2) = 0
                } else {
                    -(std::f64::consts::PI / (m as f64)).cos()
                };
                b[i][j] = val;
                b[j][i] = val;
            }
        }
        b
    }

    /// Apply a word (sequence of generator indices) to a vector.
    pub fn apply_word(&self, word: &[usize], v: &[f64]) -> Vec<f64> {
        let mut result = v.to_vec();
        for &s in word {
            result = mat_vec(&self.reflections[s], &result);
        }
        result
    }

    /// Compose reflections: get the matrix for a word.
    pub fn word_matrix(&self, word: &[usize]) -> Vec<Vec<f64>> {
        let n = self.coxeter.rank;
        let mut result = identity(n);
        for &s in word {
            result = mat_mul(&self.reflections[s], &result);
        }
        result
    }

    /// Compute the determinant of the matrix for a word.
    /// For reflections (Coxeter group), this is always ±1.
    pub fn word_determinant(&self, word: &[usize]) -> f64 {
        let m = self.word_matrix(word);
        det(&m)
    }

    /// Check if a word represents the identity element.
    pub fn is_identity(&self, word: &[usize]) -> bool {
        let m = self.word_matrix(word);
        is_identity_matrix(&m)
    }
}

/// Compute the reflection matrix for generator `i`.
fn compute_reflection(coxeter: &CoxeterMatrix, i: usize) -> Vec<Vec<f64>> {
    let n = coxeter.rank;
    let mut r = identity(n);
    // s_i(e_j) = e_j - 2*B(e_i,e_j)/B(e_i,e_i) * e_i
    // With B(e_i,e_i) = 1:
    for j in 0..n {
        let m_ij = coxeter.entries[i][j];
        let b_ij = if i == j {
            1.0
        } else if m_ij == 0 {
            1.0
        } else if m_ij == 2 {
            0.0
        } else {
            -(std::f64::consts::PI / (m_ij as f64)).cos()
        };
        r[j][i] -= 2.0 * b_ij;
    }
    r
}

fn identity(n: usize) -> Vec<Vec<f64>> {
    let mut m = vec![vec![0.0; n]; n];
    for i in 0..n { m[i][i] = 1.0; }
    m
}

fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut c = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

fn mat_vec(m: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    m.iter().map(|row| row.iter().zip(v.iter()).map(|(a, b)| a * b).sum()).collect()
}

fn det(m: &[Vec<f64>]) -> f64 {
    let n = m.len();
    if n == 1 { return m[0][0]; }
    if n == 2 { return m[0][0] * m[1][1] - m[0][1] * m[1][0]; }
    let mut d = 0.0;
    for j in 0..n {
        let minor = compute_minor(m, 0, j);
        d += if j % 2 == 0 { 1.0 } else { -1.0 } * m[0][j] * det(&minor);
    }
    d
}

fn compute_minor(m: &[Vec<f64>], row: usize, col: usize) -> Vec<Vec<f64>> {
    let n = m.len();
    (0..n).filter(|&i| i != row).map(|i| {
        (0..n).filter(|&j| j != col).map(|j| m[i][j]).collect()
    }).collect()
}

fn is_identity_matrix(m: &[Vec<f64>]) -> bool {
    let n = m.len();
    for i in 0..n {
        for j in 0..n {
            let expected = if i == j { 1.0 } else { 0.0 };
            if (m[i][j] - expected).abs() > 1e-10 {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reflection_a2() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        assert_eq!(repr.reflections.len(), 2);
    }

    #[test]
    fn test_reflection_determinant() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        // Each reflection has determinant -1
        for r in &repr.reflections {
            let d = det(r);
            assert!((d + 1.0).abs() < 1e-10, "Reflection det should be -1, got {}", d);
        }
    }

    #[test]
    fn test_bilinear_form_a2() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        let b = repr.bilinear_form();
        // B(e_0, e_0) = 1, B(e_0, e_1) = -cos(π/3) = -0.5
        assert!((b[0][0] - 1.0).abs() < 1e-10);
        assert!((b[0][1] + 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_identity_word() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        assert!(repr.is_identity(&[]));
        // s_0 * s_0 = identity
        assert!(repr.is_identity(&[0, 0]));
    }

    #[test]
    fn test_apply_word() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        let v = vec![1.0, 0.0];
        let result = repr.apply_word(&[0], &v);
        assert_eq!(result.len(), 2);
        // Should change the vector
        assert!((result[0] - v[0]).abs() > 0.01 || (result[1] - v[1]).abs() > 0.01);
    }

    #[test]
    fn test_word_determinant_identity() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        let d = repr.word_determinant(&[]);
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_word_determinant_reflection() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        let d = repr.word_determinant(&[0]);
        assert!((d + 1.0).abs() < 1e-10);
    }
}
