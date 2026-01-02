# kreep

A minimal Rust library for abstract algebraic structures and prime field arithmetic.

## Overview

kreep provides trait definitions for common algebraic structures (rings, groups, fields) and a concrete implementation of prime fields `Fp<P>` using const generics.

## Features

- **Zero dependencies** - only uses Rust standard library
- **Const generics** - prime modulus is a compile-time constant for type safety
- **Safe arithmetic** - uses `u128` internally to prevent overflow
- **Full operator support** - `+`, `-`, `*`, `/`, and negation
- **Efficient inverse** - extended Euclidean algorithm for multiplicative inverses

## Quick Start

```rust
use kreep::{Field, Fp};

// Define a prime field with modulus 17
type F17 = Fp<17>;

fn main() {
    let a = F17::new(5);
    let b = F17::new(9);

    println!("a + b = {:?}", a + b);           // Fp<17>(14)
    println!("a * b = {:?}", a * b);           // Fp<17>(11)
    println!("a^-1 = {:?}", a.inverse());      // Some(Fp<17>(7))
    println!("a / b = {:?}", a / b);           // Fp<17>(3)
}
```

## Running Examples

```bash
cargo run --example f17
cargo run --example playground
```

## Running Tests

```bash
cargo test
```

## Notes

- The modulus `P` should be prime for `Fp<P>` to be a valid field
- Maximum modulus size is `u64::MAX`
- Division by zero will panic

## License

MIT
