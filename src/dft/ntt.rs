use rand::Rng;

use crate::dft::DFT;
use crate::mod_ops::ModOps;

pub struct NttBuilder<O> {
    q: O,
    n: usize,
}

#[derive(Debug)]
pub enum Error {
    NoPrimaryRoot,
    NotPowerOfTwoN,
}

impl<O> NttBuilder<O> {
    pub fn new(q: O, n: usize) -> Self {
        Self { q, n }
    }
}

impl NttBuilder<u64> {
    pub fn build(self) -> Result<Table<u64>, Error> {
        let q = self.q;
        let n = self.n;
        
        let psi = self.find_nth_root( 2*self.n as u64).expect("cannot find psi");
        let inv_psi = psi.mod_inv(self.q);

        let (psi_pow, psi_inv_pow) = calculate_wiggle(psi, inv_psi, q, n);
        
        Ok(Table {
            q,
            n,
            psi_pow,
            psi_inv_pow,
        })
    }

    fn find_nth_root(&self, n: u64) -> Result<u64, Error> {
        if !n.is_power_of_two() {
            return Err(Error::NotPowerOfTwoN);
        }

        if (self.q - 1) % n != 0 {
            return Err(Error::NoPrimaryRoot);
        }

        let mut rng = rand::thread_rng();

        let exp = (self.q - 1) / n;

        for _ in 0..n {
            let num = rng.gen_range(1..self.q).mod_exp(exp, self.q);
            
            if num.mod_exp(n/2, self.q) != 1 {
                return Ok(num);
            }
        }

        Err(Error::NoPrimaryRoot)
    }
}

pub struct Table<O> {
    /// NTT friendly prime modulus
    q: O,
    /// number of terms in polynomial
    n: usize,
    /// powers of 2n-th root of unity in bit reversed order
    psi_pow: Box<[O]>,
    /// powers of inverse 2n-th root of unity root in bit reversed order
    psi_inv_pow: Box<[O]>,
}

impl Default for Table<u64> {
    fn default() -> Self {
        let q = 0x1fffffffffe00001u64;
        let psi = 0x15eb043c7aa2b01fu64;
        let n = 16;

        let inv_psi = psi.mod_inv(q);

        let (psi_pow, psi_inv_pow) = calculate_wiggle(psi, inv_psi, q, n);

        Self {
            q,
            n,
            psi_pow,
            psi_inv_pow,
        }
    }
}

/// Calculates wiggle factor values.
/// Returns a tupple: (boxed slice of powers of psi, boxed slice of powers of inverse psi)
fn calculate_wiggle(psi:u64, inv_psi: u64, q: u64, n: usize) -> (Box<[u64]>,Box<[u64]>) {
    let mut psi_pow = vec![0; n];
    let mut psi_inv_pow = vec![0; n];

    let mut temp_psi = 1;
    let mut temp_inv_psi = 1;
    let bit_start = n.leading_zeros()+1;

    for i in 0..n {
        let index = i.reverse_bits() >> bit_start;
        psi_pow[index] = temp_psi;
        psi_inv_pow[index] = temp_inv_psi;
        
        temp_psi = temp_psi.mod_mul(psi, q);
        temp_inv_psi = temp_inv_psi.mod_mul(inv_psi, q);
    }

    (psi_pow.into_boxed_slice(), psi_inv_pow.into_boxed_slice())
}

impl DFT<u64> for Table<u64> {
    /// NTT forward routine
    ///
    /// - `a`: vector with each element in range `[0, q)`
    fn forward_inplace(&self, a: &mut [u64]) {
        self.forward_inplace_core::<false>(a)
    }

    /// NTT forward lazy routine
    ///
    /// - `a`: vector with each element in range `[0, 2q)`
    fn forward_inplace_lazy(&self, a: &mut [u64]) {
        self.forward_inplace_core::<true>(a)
    }

    /// NTT backward routine
    ///
    /// - `a`: vector with each element in range `[0, q)`
    fn backward_inplace(&self, a: &mut [u64]) {
        self.backward_inplace_core::<false>(a)
    }

    /// NTT backward lazy routine
    ///
    /// - `a`: vector with each element in range `[0, 2q)`
    fn backward_inplace_lazy(&self, a: &mut [u64]) {
        self.backward_inplace_core::<true>(a)
    }
}

impl Table<u64> {
    fn forward_inplace_core<const LAZY: bool>(&self, a: &mut [u64]) {
        let mut t = self.n;
        let q = self.q;
        let mut m = 1;
        while m < self.n {
            t = t / 2;
            for i in 0..m {
                let j1 = 2 * i * t;
                let j2 = j1 + t - 1;
                let s = self.psi_pow[m + i];
                for j in j1..=j2 {
                    let u = a[j];
                    let v = a[j + t].mod_mul(s, q);
                    a[j] = u.mod_add(v, q);
                    a[j + t] = u.mod_sub(v, q);
                }
            }
            m = m * 2;
        }
    }

    fn backward_inplace_core<const LAZY: bool>(&self, a: &mut [u64]) {
        let q = self.q;
        let mut t = 1;
        let mut m = self.n;

        while m > 1 {
            let mut j1 = 0;
            let h = m/2;
            for i in 0..h {
                let j2 = j1 + t -1;
                let s = self.psi_inv_pow[h+i];
                for j in j1..=j2 {
                    let u = a[j];
                    let v = a[j+t];
                    a[j] = u.mod_add(v,q);
                    a[j+t] = u.mod_sub(v, q).mod_mul(s, q);
                }
                j1 += 2*t;
            }
            t *=2;
            m /= 2;
        }

        let n_inv = (self.n as u64).mod_inv(q);
        for i in 0..self.n {
            a[i] = a[i].mod_mul(n_inv, q);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    /// Tests if `a == INTT(NTT(a))` for a lazy calculation with known numbers.
    #[test]
    fn self_inverse_known_lazy() {
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
    fn self_inverse_known_greedy() {
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
    fn self_inverse_random_greedy() {
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
    fn self_inverse_random_lazy() {
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

    /// Checks if lazy polynomial multiplication result matches with manual calculation.
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

    /// Checks if greedy polynomial multiplication result matches with manual calculation.
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

    fn elementwise_mod_mul(a: &mut [u64], b: &mut [u64], q: u64) {
        for i in 0..a.len() {
            a[i] = ((a[i] as u128 * b[i] as u128) % q as u128) as u64;
        }
    }

    fn _manual_polynomial_mod_mul(){
        // TODO
    }

    /// Test for the provided NTT friendly prime and its 2^17-th primitive root to uphold `a == INTT(NTT(a))`.
    /// Uses a randomly generated polynomial.
    #[test]
    fn self_inverse_default_ntt() {
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

    /// Test for polynomial multiplication using the provided NTT friendly prime and its 2^17-th primitive root.
    /// Uses randomly generated polynomials of 15-th degree.
    #[test]
    #[ignore]
    fn polynomial_mul_default_ntt() {
        let ntt = Table::default();
        let q = ntt.q;
        let _n = ntt.n;

        let mut a = vec![1u64,2,3,4];
        let mut b = vec![5u64,6,7,8];
        let manual = vec![7625u64, 7645, 2, 60];
        
        ntt.forward_inplace(&mut a);
        ntt.forward_inplace(&mut b);
        
        elementwise_mod_mul(&mut a, &mut b, q);

        ntt.backward_inplace(&mut a);

        assert_eq!(manual, a);
    }
}