# coxeter-group-rs

Coxeter group computations: Coxeter matrix, reflection representation, word reduction, and Bruhat order.

## Features

- **Matrix**: Coxeter matrix construction (type A, B, I₂), group enumeration
- **Reflection**: Geometric reflection representation with bilinear forms
- **Word**: Word reduction with braid relations and commuting generators
- **Bruhat**: Bruhat order with subword property
- **Graph**: Coxeter graph, connected components, DOT output

Pure Rust, no external dependencies.

## Usage

```rust
use coxeter_group_rs::CoxeterMatrix;

let m = CoxeterMatrix::type_a(3);
assert!(m.is_finite());
assert_eq!(m.group_order(), Some(24)); // S_4
```

License: MIT OR Apache-2.0
