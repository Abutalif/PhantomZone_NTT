use app::dft::{ntt::NttBuilder, DFT};
use num_traits::PrimInt;
use rand::{thread_rng, Rng};

fn main() {
    // known value polynomials
    let q = 7681;
    let n = 4;
    let ntt = NttBuilder::new(q, n).build().expect("Cannot build");
    let mut a = vec![1u64,2,3,4];
    let mut b = vec![5u64,6,7,8];
    
    ntt.forward_inplace(&mut a);
    ntt.forward_inplace(&mut b);
    
    let mut c = vec!(0; n);
    for i in 0..n {
        c[i] = ((a[i] as u128 * b[i] as u128) % q as u128) as u64;
    }
    ntt.backward_inplace(&mut c);

    println!("mine: {c:?}");


    let q = 0x1fffffffffe00001u64;
    let n = 16;
    let ntt = NttBuilder::new(q, n).build().expect("Cannot build");

    // known element slice.
    let mut a:Vec<u64> = (1..=n as u64).collect();
    let a_clone = a.clone();
    println!("a = {a:?}");
    ntt.forward_inplace(&mut a);
    ntt.backward_inplace(&mut a);
    println!("INTT(NTT(a)) = {a:?}"); // Ok, So far, this works.
    assert_eq!(a, a_clone);

    // random filled slice
    let mut rng = thread_rng();
    let mut b = vec![0u64; n].into_boxed_slice();

    for el in b.iter_mut() {
        *el = rng.gen_range(0..q);
    }

    let b_clone = b.clone();
    println!("b = {b:?}");
    ntt.forward_inplace(&mut b);
    ntt.backward_inplace(&mut b);
    println!("INTT(NTT(b)) = {b:?}"); // Ok, So far, this works.

    assert_eq!(b, b_clone);

    // Special prime
    let spec_q:u64 = ((((1 << 63) - (1 << 31)))<<1) + 1;
    let spec_bit = 0xffffffff00000001u64;

    assert_eq!(spec_q, spec_bit);

    println!("Special Q: {spec_q}; spocal Bit: {spec_bit}");


}