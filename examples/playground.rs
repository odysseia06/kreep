use kreep::Fp;
use std::collections::HashSet;

type F17 = Fp<17>;

fn main() {
    // Create a set containing all elements of Fp<17>
    let field_elements: HashSet<F17> = (0..17).map(|i| F17::new(i)).collect();

    println!("Number of elements in Fp<17>: {}", field_elements.len());

    // Print all elements
    let mut elements: Vec<_> = field_elements.iter().collect();
    elements.sort_by_key(|e| e.value());

    println!("All elements:");
    for elem in elements {
        println!("  {:?}", elem);
    }
}
