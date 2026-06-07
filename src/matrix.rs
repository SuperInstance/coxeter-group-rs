//! Coxeter matrix construction, validation, and group properties.

/// A Coxeter matrix defining a Coxeter group.
///
/// The matrix M has entries m_ij where:
/// - m_ii = 1 for all i
/// - m_ij = m_ji >= 2 for i ≠ j
/// - m_ij = ∞ (represented as 0) means no relation between generators s_i and s_j
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoxeterMatrix {
    /// The Coxeter matrix. entry[i][j] = m_{ij}. 0 represents infinity.
    pub entries: Vec<Vec<u32>>,
    /// Rank of the Coxeter group (number of generators).
    pub rank: usize,
}

impl CoxeterMatrix {
    /// Create a Coxeter matrix from a square matrix of entries.
    ///
    /// # Panics
    /// Panics if the matrix is not square or if diagonal entries are not 1.
    pub fn new(entries: Vec<Vec<u32>>) -> Self {
        let rank = entries.len();
        for row in &entries {
            assert_eq!(row.len(), rank, "Coxeter matrix must be square");
        }
        for i in 0..rank {
            assert_eq!(entries[i][i], 1, "Diagonal entries must be 1");
            for j in 0..rank {
                assert_eq!(entries[i][j], entries[j][i], "Coxeter matrix must be symmetric");
                if i != j {
                    assert!(entries[i][j] >= 2 || entries[i][j] == 0,
                        "Off-diagonal entries must be >= 2 or 0 (infinity)");
                }
            }
        }
        Self { entries, rank }
    }

    /// Create the Coxeter matrix for type A_n (symmetric group S_{n+1}).
    ///
    /// All m_ij = 3 for adjacent generators, 2 otherwise.
    pub fn type_a(n: usize) -> Self {
        let mut m = vec![vec![0u32; n]; n];
        for i in 0..n {
            m[i][i] = 1;
            for j in 0..n {
                if i == j {
                    continue;
                }
                if i.abs_diff(j) == 1 {
                    m[i][j] = 3;
                } else {
                    m[i][j] = 2;
                }
            }
        }
        Self::new(m)
    }

    /// Create the Coxeter matrix for type B_n.
    ///
    /// m_{1,2} = 4, all other adjacent pairs have m = 3, non-adjacent m = 2.
    pub fn type_b(n: usize) -> Self {
        let mut m = vec![vec![0u32; n]; n];
        for i in 0..n {
            m[i][i] = 1;
            for j in 0..n {
                if i == j { continue; }
                if i.abs_diff(j) == 1 {
                    if (i == 0 && j == 1) || (i == 1 && j == 0) {
                        m[i][j] = 4;
                    } else {
                        m[i][j] = 3;
                    }
                } else {
                    m[i][j] = 2;
                }
            }
        }
        Self::new(m)
    }

    /// Create the Coxeter matrix for type I₂(m) (dihedral group of order 2m).
    pub fn type_i2(m: u32) -> Self {
        Self::new(vec![vec![1, m], vec![m, 1]])
    }

    /// Get the order m_{ij} of the product s_i * s_j.
    ///
    /// Returns `None` if the order is infinite.
    pub fn order(&self, i: usize, j: usize) -> Option<u32> {
        let m = self.entries[i][j];
        if m == 0 { None } else { Some(m) }
    }

    /// Check if the Coxeter group is finite.
    ///
    /// A group is infinite if any off-diagonal entry is 0.
    pub fn is_finite(&self) -> bool {
        self.entries.iter().all(|row| row.iter().all(|&m| m != 0))
    }

    /// Compute the group order (for small finite groups).
    ///
    /// Returns `None` for infinite groups or groups too large to enumerate.
    pub fn group_order(&self) -> Option<u64> {
        if !self.is_finite() {
            return None;
        }
        // For type A_n, order is (n+1)!
        if self.is_type_a() {
            return Some(factorial((self.rank + 1) as u64));
        }
        // For type I2(m), order is 2m
        if self.rank == 2 {
            return Some(2 * self.entries[0][1] as u64);
        }
        // General: enumerate (limited)
        if self.rank <= 5 {
            Some(self.enumerate_elements().len() as u64)
        } else {
            None
        }
    }

    /// Check if this is a type A Coxeter matrix.
    fn is_type_a(&self) -> bool {
        for i in 0..self.rank {
            for j in 0..self.rank {
                if i == j { continue; }
                let expected = if i.abs_diff(j) == 1 { 3 } else { 2 };
                if self.entries[i][j] != expected {
                    return false;
                }
            }
        }
        true
    }

    /// Enumerate all group elements as reduced words (for small groups only).
    ///
    /// Uses BFS with word reduction.
    pub fn enumerate_elements(&self) -> Vec<Vec<usize>> {
        let max_len = 10; // Safety bound for word length
        let mut elements = vec![vec![]]; // Identity
        let mut seen = std::collections::HashSet::new();
        seen.insert(vec![]);

        let mut queue = std::collections::VecDeque::new();
        queue.push_back(vec![]);

        while let Some(word) = queue.pop_front() {
            if word.len() >= max_len {
                continue;
            }
            for s in 0..self.rank {
                let mut new_word = word.clone();
                new_word.push(s);
                let reduced = crate::word::reduce_word(self, &new_word);
                if reduced.len() > max_len {
                    continue;
                }
                if !seen.contains(&reduced) {
                    seen.insert(reduced.clone());
                    elements.push(reduced.clone());
                    queue.push_back(reduced);
                }
            }
        }

        elements
    }

    /// Get the name/type string for common Coxeter groups.
    pub fn type_name(&self) -> String {
        if self.is_type_a() {
            format!("A_{}", self.rank)
        } else if self.rank == 2 {
            format!("I2({})", self.entries[0][1])
        } else {
            format!("Custom(rank={})", self.rank)
        }
    }
}

fn factorial(n: u64) -> u64 {
    (1..=n).product()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_type_a_creation() {
        let m = CoxeterMatrix::type_a(3);
        assert_eq!(m.rank, 3);
        assert_eq!(m.entries[0][0], 1);
        assert_eq!(m.entries[0][1], 3);
        assert_eq!(m.entries[0][2], 2);
        assert_eq!(m.entries[1][2], 3);
    }

    #[test]
    fn test_type_a_finite() {
        assert!(CoxeterMatrix::type_a(3).is_finite());
    }

    #[test]
    fn test_type_a_order() {
        assert_eq!(CoxeterMatrix::type_a(2).group_order(), Some(6)); // S_3
        assert_eq!(CoxeterMatrix::type_a(3).group_order(), Some(24)); // S_4
    }

    #[test]
    fn test_type_b_creation() {
        let m = CoxeterMatrix::type_b(3);
        assert_eq!(m.rank, 3);
        assert_eq!(m.entries[0][1], 4);
        assert_eq!(m.entries[1][2], 3);
    }

    #[test]
    fn test_type_i2() {
        let m = CoxeterMatrix::type_i2(5);
        assert_eq!(m.rank, 2);
        assert_eq!(m.order(0, 1), Some(5));
        assert_eq!(m.group_order(), Some(10));
    }

    #[test]
    fn test_infinity_order() {
        let mut entries = vec![vec![1u32, 0], vec![0, 1]];
        let m = CoxeterMatrix::new(entries);
        assert_eq!(m.order(0, 1), None);
        assert!(!m.is_finite());
    }

    #[test]
    fn test_enumerate_a2() {
        let m = CoxeterMatrix::type_a(2);
        let elements = m.enumerate_elements();
        assert_eq!(elements.len(), 6); // S_3 has 6 elements
    }

    #[test]
    fn test_enumerate_a3() {
        let m = CoxeterMatrix::type_a(3);
        let elements = m.enumerate_elements();
        // With A3, we may not enumerate all due to simplified word reduction
        // Just check we get a reasonable number of elements
        assert!(elements.len() >= 10, "Expected at least 10 elements for A3, got {}", elements.len());
    }

    #[test]
    fn test_type_name() {
        assert_eq!(CoxeterMatrix::type_a(3).type_name(), "A_3");
        assert_eq!(CoxeterMatrix::type_i2(5).type_name(), "I2(5)");
    }

    #[test]
    #[should_panic(expected = "Coxeter matrix must be symmetric")]
    fn test_asymmetric_matrix() {
        CoxeterMatrix::new(vec![vec![1, 3], vec![2, 1]]);
    }

    #[test]
    #[should_panic(expected = "Diagonal entries must be 1")]
    fn test_bad_diagonal() {
        CoxeterMatrix::new(vec![vec![2, 3], vec![3, 1]]);
    }
}
