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

    /// Generate all reflections in the group (conjugates of simple reflections).
    ///
    /// For a Coxeter group of rank n, returns all reflections
    /// w * s_i * w^{-1} for all group elements w and generators s_i.
    pub fn all_reflections(&self) -> Vec<Vec<Vec<f64>>> {
        let mut reflections = Vec::new();
        let mut seen = std::collections::HashSet::new();

        // Start with simple reflections
        for r in &self.reflections {
            let key = matrix_key_rounded(r);
            if seen.insert(key) {
                reflections.push(r.clone());
            }
        }

        // Conjugate by each known reflection to find new ones
        let mut changed = true;
        while changed {
            changed = false;
            let current = reflections.clone();
            for r1 in &current {
                for r2 in &current {
                    let conjugated = mat_mul(r1, &mat_mul(r2, r1));
                    let key = matrix_key_rounded(&conjugated);
                    if seen.insert(key) {
                        reflections.push(conjugated);
                        changed = true;
                    }
                }
            }
        }

        reflections
    }

    /// Compose two group elements (given as words) into a single word.
    /// Returns the matrix product for the concatenated word w1·w2.
    pub fn compose(&self, w1: &[usize], w2: &[usize]) -> Vec<Vec<f64>> {
        let m2 = self.word_matrix(w2);
        let m1 = self.word_matrix(w1);
        // Group composition: first apply w1, then w2 → mat_mul(m2, m1)
        mat_mul(&m2, &m1)
    }

    /// Compute the Coxeter presentation string.
    ///
    /// Returns the group presentation in terms of generators and relations.
    pub fn coxeter_presentation(&self) -> CoxeterPresentation {
        let generators: Vec<String> = (0..self.coxeter.rank)
            .map(|i| format!("s{}", i))
            .collect();

        let mut relations = vec![];
        // s_i^2 = e
        for i in 0..self.coxeter.rank {
            relations.push(CoxeterRelation::Order2(i));
        }
        // Braid relations
        for i in 0..self.coxeter.rank {
            for j in (i + 1)..self.coxeter.rank {
                let m = self.coxeter.entries[i][j];
                if m > 2 {
                    relations.push(CoxeterRelation::Braid { i, j, m });
                }
            }
        }

        CoxeterPresentation { generators, relations }
    }

    /// Compute the order of a group element represented by a word.
    /// Returns the smallest k > 0 such that w^k = identity, or None.
    pub fn element_order(&self, word: &[usize], max_k: usize) -> Option<usize> {
        let m = self.word_matrix(word);
        if is_identity_matrix(&m) {
            return Some(1);
        }
        let mut current = m.clone();
        for k in 2..=max_k {
            current = mat_mul(&m, &current);
            if is_identity_matrix(&current) {
                return Some(k);
            }
        }
        None
    }

    /// Compute the trace of the matrix for a word.
    pub fn word_trace(&self, word: &[usize]) -> f64 {
        let m = self.word_matrix(word);
        (0..m.len()).map(|i| m[i][i]).sum()
    }

    /// Check if a word represents a reflection (element conjugate to a simple reflection).
    pub fn is_reflection(&self, word: &[usize]) -> bool {
        let m = self.word_matrix(word);
        if !is_identity_matrix(&m) {
            let d = det(&m);
            (d + 1.0).abs() < 1e-8
        } else {
            false
        }
    }

    /// Compute the eigenvalues of a word's matrix.
    pub fn word_eigenvalues(&self, word: &[usize]) -> Vec<f64> {
        let m = self.word_matrix(word);
        let n = m.len();
        if n == 0 {
            return vec![];
        }
        if n == 2 {
            let trace = m[0][0] + m[1][1];
            let determinant = m[0][0] * m[1][1] - m[0][1] * m[1][0];
            let disc = (trace * trace - 4.0 * determinant).max(0.0);
            return vec![(trace + disc.sqrt()) / 2.0, (trace - disc.sqrt()) / 2.0];
        }
        // For larger matrices, use power iteration for dominant eigenvalue
        let mut v = vec![1.0; n];
        for _ in 0..200 {
            let mut w = vec![0.0; n];
            for i in 0..n {
                for j in 0..n {
                    w[i] += m[i][j] * v[j];
                }
            }
            let norm: f64 = w.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm < 1e-15 {
                break;
            }
            v = w.iter().map(|x| x / norm).collect();
        }
        // Rayleigh quotient
        let mut mv = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                mv[i] += m[i][j] * v[j];
            }
        }
        let lambda: f64 = v.iter().zip(mv.iter()).map(|(a, b)| a * b).sum();
        vec![lambda]
    }
}

/// A Coxeter presentation: generators and relations.
#[derive(Debug, Clone)]
pub struct CoxeterPresentation {
    /// Generator names.
    pub generators: Vec<String>,
    /// Relations.
    pub relations: Vec<CoxeterRelation>,
}

/// A relation in the Coxeter presentation.
#[derive(Debug, Clone, PartialEq)]
pub enum CoxeterRelation {
    /// s_i^2 = e
    Order2(usize),
    /// (s_i s_j)^m = e
    Braid { i: usize, j: usize, m: u32 },
}

impl std::fmt::Display for CoxeterPresentation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "<")?;
        for (i, g) in self.generators.iter().enumerate() {
            if i > 0 { write!(f, ", ")?; }
            write!(f, "{}", g)?;
        }
        write!(f, " | ")?;
        for (i, rel) in self.relations.iter().enumerate() {
            if i > 0 { write!(f, ", ")?; }
            match rel {
                CoxeterRelation::Order2(s) => write!(f, "s{}^2", s)?,
                CoxeterRelation::Braid { i: a, j: b, m } => {
                    write!(f, "(s{}s{})^{}", a, b, m)?;
                }
            }
        }
        write!(f, ">")
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

fn matrix_key_rounded(m: &[Vec<f64>]) -> Vec<i64> {
    m.iter().flat_map(|row| {
        row.iter().map(|&x| (x * 1e8).round() as i64)
    }).collect()
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
        assert!((b[0][0] - 1.0).abs() < 1e-10);
        assert!((b[0][1] + 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_identity_word() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        assert!(repr.is_identity(&[]));
        assert!(repr.is_identity(&[0, 0]));
    }

    #[test]
    fn test_apply_word() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        let v = vec![1.0, 0.0];
        let result = repr.apply_word(&[0], &v);
        assert_eq!(result.len(), 2);
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

    #[test]
    fn test_all_reflections_a2() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        let all_ref = repr.all_reflections();
        // A2 has 6 elements, each non-identity element is either a reflection or rotation
        // Simple reflections + their conjugates
        assert!(all_ref.len() >= 2);
        // Each should have det = -1
        for r in &all_ref {
            assert!((det(r) + 1.0).abs() < 1e-8);
        }
    }

    #[test]
    fn test_compose_words() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        let composed = repr.compose(&[0], &[1]);
        let expected = repr.word_matrix(&[0, 1]);
        for i in 0..2 {
            for j in 0..2 {
                assert!((composed[i][j] - expected[i][j]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_coxeter_presentation_a2() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        let pres = repr.coxeter_presentation();
        assert_eq!(pres.generators.len(), 2);
        // s0^2 = e, s1^2 = e, (s0*s1)^3 = e
        assert_eq!(pres.relations.len(), 3);
    }

    #[test]
    fn test_coxeter_presentation_display() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        let pres = repr.coxeter_presentation();
        let s = format!("{}", pres);
        assert!(s.contains("s0"));
        assert!(s.contains("s1"));
    }

    #[test]
    fn test_element_order_reflection() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        // A reflection has order 2
        assert_eq!(repr.element_order(&[0], 10), Some(2));
    }

    #[test]
    fn test_element_order_identity() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        assert_eq!(repr.element_order(&[], 10), Some(1));
        assert_eq!(repr.element_order(&[0, 0], 10), Some(1));
    }

    #[test]
    fn test_word_trace() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        // Identity has trace = rank
        let t = repr.word_trace(&[]);
        assert!((t - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_is_reflection() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        assert!(repr.is_reflection(&[0]));
        assert!(repr.is_reflection(&[1]));
        assert!(!repr.is_reflection(&[]));
    }

    #[test]
    fn test_word_eigenvalues() {
        let coxeter = CoxeterMatrix::type_a(2);
        let repr = ReflectionRepresentation::new(&coxeter);
        let eigs = repr.word_eigenvalues(&[]);
        // Identity has eigenvalues all 1
        assert!(eigs.iter().all(|&e| (e - 1.0).abs() < 1e-8));
    }
}
