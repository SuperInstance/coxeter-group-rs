//! Coxeter graph / Dynkin diagram utilities.

use crate::matrix::CoxeterMatrix;

/// An edge in the Coxeter graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoxeterEdge {
    pub from: usize,
    pub to: usize,
    /// The Coxeter matrix entry m_{ij}. 0 means infinity.
    pub order: u32,
}

/// A Coxeter graph (dual to the Dynkin diagram).
#[derive(Debug, Clone)]
pub struct CoxeterGraph {
    /// Number of vertices (generators).
    pub num_vertices: usize,
    /// Edges in the graph (only edges with m > 2).
    pub edges: Vec<CoxeterEdge>,
}

impl CoxeterGraph {
    /// Build the Coxeter graph from a Coxeter matrix.
    pub fn from_matrix(coxeter: &CoxeterMatrix) -> Self {
        let mut edges = vec![];
        for i in 0..coxeter.rank {
            for j in (i + 1)..coxeter.rank {
                let m = coxeter.entries[i][j];
                if m > 2 {
                    edges.push(CoxeterEdge { from: i, to: j, order: m });
                }
            }
        }
        Self {
            num_vertices: coxeter.rank,
            edges,
        }
    }

    /// Get the edge between vertices i and j, if any.
    pub fn edge(&self, i: usize, j: usize) -> Option<&CoxeterEdge> {
        self.edges.iter().find(|e| (e.from == i && e.to == j) || (e.from == j && e.to == i))
    }

    /// Get the neighbors of vertex v.
    pub fn neighbors(&self, v: usize) -> Vec<usize> {
        self.edges
            .iter()
            .filter(|e| e.from == v || e.to == v)
            .map(|e| if e.from == v { e.to } else { e.from })
            .collect()
    }

    /// Number of edges.
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    /// Check if the graph is connected.
    pub fn is_connected(&self) -> bool {
        if self.num_vertices <= 1 {
            return true;
        }
        let mut visited = vec![false; self.num_vertices];
        let mut stack = vec![0];
        visited[0] = true;
        let mut count = 1;

        while let Some(v) = stack.pop() {
            for &n in &self.neighbors(v) {
                if !visited[n] {
                    visited[n] = true;
                    count += 1;
                    stack.push(n);
                }
            }
        }

        count == self.num_vertices
    }

    /// Compute the connected components.
    pub fn components(&self) -> Vec<Vec<usize>> {
        let mut visited = vec![false; self.num_vertices];
        let mut components = vec![];

        for start in 0..self.num_vertices {
            if visited[start] {
                continue;
            }
            let mut component = vec![];
            let mut stack = vec![start];
            visited[start] = true;

            while let Some(v) = stack.pop() {
                component.push(v);
                for &n in &self.neighbors(v) {
                    if !visited[n] {
                        visited[n] = true;
                        stack.push(n);
                    }
                }
            }
            components.push(component);
        }
        components
    }

    /// Format as a Graphviz DOT string.
    pub fn to_dot(&self) -> String {
        let mut dot = String::from("graph Coxeter {\n");
        for e in &self.edges {
            let label = if e.order == 0 {
                "∞".to_string()
            } else if e.order == 3 {
                String::new() // Default edge label for m=3
            } else {
                e.order.to_string()
            };
            if label.is_empty() {
                dot.push_str(&format!("  {} -- {};\n", e.from, e.to));
            } else {
                dot.push_str(&format!("  {} -- {} [label=\"{}\"];\n", e.from, e.to, label));
            }
        }
        // Isolated vertices
        let connected: std::collections::HashSet<usize> = self.edges
            .iter()
            .flat_map(|e| [e.from, e.to])
            .collect();
        for v in 0..self.num_vertices {
            if !connected.contains(&v) {
                dot.push_str(&format!("  {};\n", v));
            }
        }
        dot.push_str("}\n");
        dot
    }
}

/// Classify the Coxeter graph type (finite irreducible).
pub fn classify(coxeter: &CoxeterMatrix) -> String {
    let graph = CoxeterGraph::from_matrix(coxeter);
    let components = graph.components();

    if components.len() == 1 {
        classify_irreducible(coxeter, &graph)
    } else {
        components
            .iter()
            .map(|_| "A".to_string())
            .collect::<Vec<_>>()
            .join(" × ")
    }
}

fn classify_irreducible(coxeter: &CoxeterMatrix, graph: &CoxeterGraph) -> String {
    let n = coxeter.rank;

    // Check for A_n pattern: path graph with all edges labeled 3
    if graph.edges.iter().all(|e| e.order == 3) && graph.edges.len() == n - 1 && n > 1 {
        let mut degrees = vec![0usize; n];
        for e in &graph.edges {
            degrees[e.from] += 1;
            degrees[e.to] += 1;
        }
        let leaves = degrees.iter().filter(|&&d| d == 1).count();
        if leaves == 2 && degrees.iter().all(|&d| d <= 2) {
            return format!("A_{}", n);
        }
    }

    // Check for I2(m)
    if n == 2 {
        let m = coxeter.entries[0][1];
        return if m == 3 {
            "A_2".to_string()
        } else if m == 4 {
            "B_2".to_string()
        } else {
            format!("I2({})", m)
        };
    }

    // Check for B_n: path with one edge labeled 4
    if n > 2 {
        let edges_4: Vec<_> = graph.edges.iter().filter(|e| e.order == 4).collect();
        let edges_3: Vec<_> = graph.edges.iter().filter(|e| e.order == 3).collect();
        if edges_4.len() == 1 && edges_3.len() == n - 2 && graph.edges.len() == n - 1 {
            return format!("B_{}", n);
        }
    }

    format!("Custom(rank={})", n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_graph_from_a3() {
        let coxeter = CoxeterMatrix::type_a(3);
        let graph = CoxeterGraph::from_matrix(&coxeter);
        assert_eq!(graph.num_vertices, 3);
        assert_eq!(graph.num_edges(), 2); // 0-1 and 1-2
    }

    #[test]
    fn test_graph_connected() {
        let coxeter = CoxeterMatrix::type_a(3);
        let graph = CoxeterGraph::from_matrix(&coxeter);
        assert!(graph.is_connected());
    }

    #[test]
    fn test_graph_neighbors() {
        let coxeter = CoxeterMatrix::type_a(3);
        let graph = CoxeterGraph::from_matrix(&coxeter);
        let neighbors = graph.neighbors(1);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&0));
        assert!(neighbors.contains(&2));
    }

    #[test]
    fn test_graph_edge() {
        let coxeter = CoxeterMatrix::type_a(3);
        let graph = CoxeterGraph::from_matrix(&coxeter);
        let e = graph.edge(0, 1).unwrap();
        assert_eq!(e.order, 3);
        assert!(graph.edge(0, 2).is_none()); // m_{0,2} = 2, no edge
    }

    #[test]
    fn test_classify_a3() {
        let coxeter = CoxeterMatrix::type_a(3);
        assert_eq!(classify(&coxeter), "A_3");
    }

    #[test]
    fn test_classify_b3() {
        let coxeter = CoxeterMatrix::type_b(3);
        assert_eq!(classify(&coxeter), "B_3");
    }

    #[test]
    fn test_classify_i2() {
        let coxeter = CoxeterMatrix::type_i2(5);
        assert_eq!(classify(&coxeter), "I2(5)");
    }

    #[test]
    fn test_to_dot() {
        let coxeter = CoxeterMatrix::type_a(2);
        let graph = CoxeterGraph::from_matrix(&coxeter);
        let dot = graph.to_dot();
        assert!(dot.contains("graph Coxeter"));
        assert!(dot.contains("0 -- 1"));
    }
}
