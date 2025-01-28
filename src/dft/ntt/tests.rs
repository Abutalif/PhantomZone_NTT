use super::*;
use rand::thread_rng;

/// Computes element-wise multiplication of slices *a* and *b* in the ring Z<sub>q</sub>(x)/(x<sup>n</sup>+1).
/// Saves the result in-place of *a*.
fn elementwise_mod_mul(a: &mut [u64], b: &mut [u64], q: u64) {
    for i in 0..a.len() {
        a[i] = a[i].mod_mul(b[i], q);
    }
}

/// Negacyclic convolution of polynomials *a* and *b*. 
fn manual_negacyclic_convolution(a: &[u64], b: &[u64], q: u64) -> Vec<u64> {
    let n = a.len();
    let mut res = vec![0; n];
    for k in 0..n {
        for i in 0..=k {
            res[k] = res[k].mod_add(a[i].mod_mul(b[k-i], q), q);
        } 

        for i in k+1..n {
            res[k] = res[k].mod_sub(a[i].mod_mul(b[k+n-i], q), q);
        }
    }

    res
}

/// Tests if `a == INTT(NTT(a))` for a lazy calculation with known numbers.
#[test]
fn forward_reverse_known_lazy() {
    let q = 0x1fffffffffe00001u64;
    let n = 16;
    let ntt = NttBuilder::new(q, n).build().expect("Cannot build");

    // known element slice.
    let mut a:Vec<u64> = (1..=n as u64).collect();
    let a_clone = a.clone();

    ntt.forward_inplace_lazy(&mut a);
    ntt.backward_inplace_lazy(&mut a);

    assert_eq!(a, a_clone);
}

/// Tests if `a == INTT(NTT(a))` for a greedy calculation with known numbers.
#[test]
fn forward_reverse_known_greedy() {
    let q = 0x1fffffffffe00001u64;
    let n = 16;
    let ntt = NttBuilder::new(q, n).build().expect("Cannot build");

    // known element slice.
    let mut a:Vec<u64> = (1..=n as u64).collect();
    let a_clone = a.clone();

    ntt.forward_inplace(&mut a);
    ntt.backward_inplace(&mut a);

    assert_eq!(a, a_clone);
}

/// Tests if `a == INTT(NTT(a))` for a greedy calculation with random numbers.
#[test]
fn forward_reverse_random_greedy() {
    let q = 0x1fffffffffe00001u64;
    let n = 16;
    let ntt = NttBuilder::new(q, n).build().expect("Cannot build");

    let mut rng = thread_rng();
    let mut b = vec![0u64; n].into_boxed_slice();

    for el in b.iter_mut() {
        *el = rng.gen_range(0..q);
    }

    let b_clone = b.clone();
    ntt.forward_inplace(&mut b);
    ntt.backward_inplace(&mut b);

    assert_eq!(b, b_clone);
}

/// Tests if `a == INTT(NTT(a))` for a lazy calculation with random numbers.
#[test]
fn forward_reverse_random_lazy() {
    let q = 0x1fffffffffe00001u64;
    let n = 16;
    let ntt = NttBuilder::new(q, n).build().expect("Cannot build");

    let mut rng = thread_rng();
    let mut b = vec![0u64; n].into_boxed_slice();

    for el in b.iter_mut() {
        *el = rng.gen_range(0..q);
    }

    let b_clone = b.clone();
    ntt.forward_inplace_lazy(&mut b);
    ntt.backward_inplace_lazy(&mut b);

    assert_eq!(b, b_clone);
}

/// Checks if lazy negacyclic polynomial multiplication result matches with manual calculation with known numbers.
#[test]
fn polynomial_mul_known_lazy() {
    let q = 7681;
    let n = 4;
    let ntt = NttBuilder::new(q, n).build().expect("Cannot build");
    let mut a = vec![1u64,2,3,4];
    let mut b = vec![5u64,6,7,8];
    let manual = vec![7625u64, 7645, 2, 60];

    
    ntt.forward_inplace_lazy(&mut a);
    ntt.forward_inplace_lazy(&mut b);

    elementwise_mod_mul(&mut a, &mut b, q);
    
    ntt.backward_inplace_lazy(&mut a);

    assert_eq!(manual, a);
}

/// Checks if greedy negacyclic polynomial multiplication result matches with manual calculation with known numbers.
#[test]
fn polynomial_mul_known_greedy() {
    let q = 7681;
    let n = 4;
    let ntt = NttBuilder::new(q, n).build().expect("Cannot build");
    let mut a = vec![1u64,2,3,4];
    let mut b = vec![5u64,6,7,8];
    let manual = vec![7625u64, 7645, 2, 60];
    
    ntt.forward_inplace(&mut a);
    ntt.forward_inplace(&mut b);
    
    elementwise_mod_mul(&mut a, &mut b, q);

    ntt.backward_inplace(&mut a);

    assert_eq!(manual, a);
}

/// Checks if greedy negacyclic polynomial multiplication result matches with school-book formula computation with random numbers.
#[test]
fn polynomial_mul_random_greedy() {
    let q = 0x1fffffffffe00001u64;
    let n = 32;
    let ntt = NttBuilder::new(q, n).build().expect("Cannot build");
    let mut rng = thread_rng();
    let mut a = vec![0u64; n].into_boxed_slice();
    let mut b = vec![0u64; n].into_boxed_slice();
    for i in 0..n {
        let a_i = rng.gen_range(0..q);
        let b_i =rng.gen_range(0..q);
        a[i] = a_i;
        b[i] = b_i;
    }

    let c = manual_negacyclic_convolution(&a, &b, q).into_boxed_slice();
    
    ntt.forward_inplace(&mut a);
    ntt.forward_inplace(&mut b);
    
    elementwise_mod_mul(&mut a, &mut b, q);

    ntt.backward_inplace(&mut a);
    assert_eq!(c, a);
}

/// Checks if lazy negacyclic polynomial multiplication result matches with school-book formula computation with random numbers.
#[test]
fn polynomial_mul_random_lazy() {
    let q = 0x1fffffffffe00001u64;
    let n = 32;
    let ntt = NttBuilder::new(q, n).build().expect("Cannot build");
    
    let mut rng = thread_rng();
    let mut a = vec![0u64; n].into_boxed_slice();
    let mut b = vec![0u64; n].into_boxed_slice();
    for i in 0..n {
        let a_i = rng.gen_range(0..q);
        let b_i =rng.gen_range(0..q);
        a[i] = a_i;
        b[i] = b_i;
    }

    let c = manual_negacyclic_convolution(&a, &b, q).into_boxed_slice();
    
    ntt.forward_inplace_lazy(&mut a);
    ntt.forward_inplace_lazy(&mut b);
    
    elementwise_mod_mul(&mut a, &mut b, q);

    ntt.backward_inplace_lazy(&mut a);
    assert_eq!(c, a);
}

/// Test for the provided NTT friendly prime and its 2<sup>17</sup>-th primitive root to uphold `a == INTT(NTT(a))`.
/// Uses a randomly generated polynomial.
#[test]
fn forward_reverse_default_ntt() {
    let ntt = Table::default();
    let q = ntt.q;
    let n = ntt.n;

    let mut rng = thread_rng();
    let mut b = vec![0u64; n].into_boxed_slice();

    for el in b.iter_mut() {
        *el = rng.gen_range(0..q);
    }

    let b_clone = b.clone();
    ntt.forward_inplace_lazy(&mut b);
    ntt.backward_inplace_lazy(&mut b);
    assert_eq!(b, b_clone);
}

/// Test for polynomial multiplication using the provided NTT friendly prime and its 2<sup>17</sup>-th primitive root.
/// Uses randomly generated polynomials of *(2<sup>16</sup> - 1)*-th degree.
/// Due to the size of the polynomail the test has been marked as `#[ignore]`
#[test]
#[ignore]
fn polynomial_mul_default_ntt() {
    let ntt = Table::default();
    let n = ntt.n;
    let q = ntt.q;

    let mut rng = thread_rng();
    let mut a = vec![0u64; n].into_boxed_slice();
    let mut b = vec![0u64; n].into_boxed_slice();
    for i in 0..n {
        let a_i = rng.gen_range(0..q);
        let b_i =rng.gen_range(0..q);
        a[i] = a_i;
        b[i] = b_i;
    }

    let c = manual_negacyclic_convolution(&a, &b, q).into_boxed_slice();
    
    ntt.forward_inplace(&mut a);
    ntt.forward_inplace(&mut b);
    
    elementwise_mod_mul(&mut a, &mut b, q);

    ntt.backward_inplace(&mut a);
    assert_eq!(c, a);
}