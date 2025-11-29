use kreep::{Field, Fp, Ring};

type F17 = Fp<17>;

fn main() {
    let a = F17::new(5);
    let b = F17::new(9);

    let sum = a + b;
    let prod = a * b;
    let inv_a = a.inverse().unwrap();

    println!("a = {:?}", a);
    println!("b = {:?}", b);
    println!("a + b = {:?}", sum);
    println!("a * b = {:?}", prod);
    println!("a^-1 = {:?}", inv_a);
    println!("a * a^-1 = {:?}", a * inv_a); // should be ONE
}
