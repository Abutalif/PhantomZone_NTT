use std::ops::{Add, Mul, Sub, Div};

pub trait ModOps: Add<Output = Self> + Mul<Output = Self> + Sub<Output = Self> + Div<Output = Self> + Copy + Sized {
    fn mod_add(self, rhs: Self, modulus: Self) -> Self;
    fn mod_mul(self, rhs: Self, modulus: Self) -> Self;
    fn mod_sub(self, rhs: Self, modulus: Self) -> Self;
    fn mod_exp(self, exp: Self, modulus: Self) -> Self;
    fn mod_inv(self, modulus: Self) -> Self;
}

impl ModOps for u64 {
    fn mod_add(self, rhs: Self, modulus: Self) -> Self {
        assert!(
            modulus <= 1 << 61,
            "prime q must be smaller or equal to 2^61"
        );
        // assert!(
        //     self < modulus && rhs < modulus,
        //     "Elements in the ring Z_q must be smaller than q"
        // );

        ((self as u128 + rhs as u128) % modulus as u128) as u64
    }

    fn mod_mul(self, rhs: Self, modulus: Self) -> Self {
        assert!(
            modulus <= 1 << 61,
            "prime q must be smaller or equal to 2^61"
        );
        // assert!(
        //     self < modulus && rhs < modulus,
        //     "Elements in Zq ring must be smaller than q"
        // );
        if self <= 1<<32 && rhs <= 1<<32 {
            return (self * rhs) % modulus;
        }

        // Russian peasant multiplication
        let mut a = self as u128;
        let mut b = rhs as u128;
        let q = modulus as u128;
        let mut res = 0;

        while b > 0 {
            if b & 1 == 1 {
                res += a % q;
            }
            a *= 2;
            b /= 2;
        }
        (res%q) as u64
    }

    fn mod_sub(self, rhs: Self, modulus: Self) -> Self {
        assert!(
            modulus <= 1 << 61,
            "prime q must be smaller or equal to 2^61"
        );
        // assert!(
        //     self < modulus && rhs < modulus,
        //     "Elements in the ring Z_q must be smaller than q"
        // );
        (self + modulus - rhs) % modulus
        // let temp = self.wrapping_sub(rhs);
        // (temp.min(self.wrapping_add(modulus))) % modulus
    }

    fn mod_exp(self, exp: Self, modulus: Self) -> Self {
        assert!(
            modulus <= 1 << 61,
            "prime q must be smaller or equal to 2^61"
        );
        assert!(
            self < modulus,
            "Elements in the ring Z_q must be smaller than q"
        );

        // exponentiation by squaring
        let mut res = 1u128;
        let mut a = self as u128;
        let mut b = exp;
        let q = modulus as u128;

        while b > 0 {
            if b % 2 == 1 {
                res = res * a % q;
            }
            a = a * a % q;
            b /= 2;
        }

        res as u64
    }

    fn mod_inv(self, modulus: Self) -> Self {
        assert!(
            modulus <= 1 << 61,
            "prime q must be smaller or equal to 2^61"
        );
        assert!(
            self < modulus,
            "Elements in the ring Z_q must be smaller than q"
        );

        // Fermat's little theorem
        self.mod_exp(modulus - 2, modulus)
    }
}
