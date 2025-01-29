use std::ops::{Add, Div, Mul, Sub};

pub trait ModOps:
    Add<Output = Self> + Mul<Output = Self> + Sub<Output = Self> + Div<Output = Self> + Copy + Sized
{
    fn mod_add(self, rhs: Self, modulus: Self) -> Self;
    fn mod_mul(self, rhs: Self, modulus: Self) -> Self;
    fn mod_sub(self, rhs: Self, modulus: Self) -> Self;
    fn mod_exp(self, exp: Self, modulus: Self) -> Self;
    fn mod_inv(self, modulus: Self) -> Self;
}

// for all modular ops we will assume that self is smaller than modulus.
impl ModOps for u64 {
    fn mod_add(self, rhs: Self, modulus: Self) -> Self {
        debug_assert!(self < modulus && rhs < modulus);

        let mut sum = self + rhs;
        if sum >= modulus {
            sum -= modulus;
        }
        sum
    }

    fn mod_mul(self, rhs: Self, modulus: Self) -> Self {
        // Russian peasant multiplication
        let mut a = self;
        let mut b = rhs;
        let q = modulus;
        let mut res = 0;

        while b > 0 {
            if (b & 1) != 0 {
                res += a;
                if res >= q {
                    res -= q;
                }
            }
            a = a << 1;
            if a >= q {
                a -= q;
            }
            b = b >> 1;
        }
        res
    }

    fn mod_sub(self, rhs: Self, modulus: Self) -> Self {
        if rhs > self {
            self + modulus - rhs
        } else {
            self - rhs
        }
    }

    fn mod_exp(self, exp: Self, modulus: Self) -> Self {
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
        // Fermat's little theorem
        self.mod_exp(modulus - 2, modulus)
    }
}

#[cfg(test)]
mod tests {}
