use crate::dft::DFT;

pub struct NttBuilder<O> {
    q: O,
    n: usize,
}

#[derive(Debug)]
pub enum Error {
    NoPrimaryRoot,
    NotEvenN,
}

impl<O> NttBuilder<O> {
    pub fn new(q: O, n: usize)-> Self{
        Self {
            q,
            n,
        }
    }
}


impl NttBuilder<u64> {
    pub fn build(self) -> Result<Table<u64>, Error> {
        let psi = self.find_primitive_root();
        Ok(Table {
            q: self.q,
            n: self.n,
            psi,
            psi_pow: todo!(),
            psi_inv: todo!(),
            psi_inv_pow: todo!(),
        })
    }

    fn find_primitive_root(&self) -> u64 {
        todo!()
    }

}

pub struct Table<O> {
    /// NTT friendly prime modulus
    q: O,
    /// 2n-th root of unity
    psi: O,
    /// reversed 2n-th root of unity
    psi_inv: O,
    /// number of terms in polynomial
    n: usize,
    /// powers of 2n-th root of unity in bit reversed order
    psi_pow: Box<[O]>,
    /// powers of inverse 2n-th root of unity root in bit reversed order
    psi_inv_pow: Box<[O]>,
}

impl Default for Table<u64> {
    fn default() -> Self {
        Self { 
            q: 0x1fffffffffe00001u64,
            psi: 0x15eb043c7aa2b01fu64,
            n: 16,
            psi_pow: todo!(),
            psi_inv: todo!(),
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
        while m <t {
            t = t/2; 
            for i in 0..m {
                let j1 = 2*i*m;
                let j2 = j1 + t  - 1;
                let s = self.psi_pow[m+i];
                // TODO: overflow and modulus handling 
                for j in j1..=j2 {
                    let u = a[j];
                    let v = a[j+t]*s;
                    a[j] = (u+v)%q;
                    a[j+t] = (u-v)%q;
                }
            }
            m = m*2;
        }
    }
    fn backward_inplace_core<const LAZY: bool>(&self, a: &mut [u64]) {}
}
