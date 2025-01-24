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
        let psi = self.find_nth_root( 2*self.n as u64).expect("cannot find psi");
        let psi_inv = psi.mod_inv(self.q);
        let mut psi_pow = vec![0; self.n];
        let mut psi_inv_pow = vec![0; self.n];

        assert_eq!(psi.mod_exp(2*self.n as u64, self.q), 1);
        assert_eq!((psi as u128 * psi_inv as u128) % self.q as u128, 1);

        let mut temp_psi = 1;
        let mut temp_inv_psi = 1;
        let bit_start = self.n.leading_zeros()+1;

        for i in 0..self.n {
            let index = i.reverse_bits() >> bit_start;
            psi_pow[index] = temp_psi;
            psi_inv_pow[index] = temp_inv_psi;
            
            temp_psi = temp_psi.mod_mul(psi, self.q);
            temp_inv_psi = temp_inv_psi.mod_mul(psi_inv, self.q);
        }
        
        Ok(Table {
            q: self.q,
            n: self.n,
            psi_pow: psi_pow.into_boxed_slice(),
            psi_inv_pow: psi_inv_pow.into_boxed_slice(),
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
            let num = rng.gen_range(1..self.q);
            
            let num = num.mod_exp(exp, self.q);

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
        let _psi = 0x15eb043c7aa2b01fu64;
        let n = 16;
        Self {
            q,
            n,
            psi_pow: todo!(),
            psi_inv_pow: todo!(),
        }
    }
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
