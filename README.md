# NWC-based NTT in Rust

This project represents a Negative Wrapped Convolution based Number Theoretic Transform (NWC-NTT) implemented in Rust. It operates over the ring $Z_{Q}[X]/X^N+1$ where $N$ is a power of two and $Q$ is an NTT friendly prime satisfying $Q\equiv 1\mod 2N$[^1].

The NWC-NTT is used in cryptography, specifically post-quantum cryptography algorithms, such as the homomorphic encryption and Ring Learning with Errors (RLwE). The main use-case is the polynomial multiplication in the quotent ring.

## Implementation details

The NTT forward and reverse transforms were implemented according to Longa and Naehrig[^2] (Algorithms 1 & 2). The solution is present in two forms:
1. **Lazy** - elements of the input vector are allowed to be in the range $[0, 2q)$.
2. **Greedy** - elements of the input vector must be in the range $[0, q)$.

The default trait was implemented for a 61-bit prime `0x1fffffffffe00001` with a $2^{17}$ root of unity equal to `0x15eb043c7aa2b01f`. Although, users are not limited to this given prime. They can provide any NTT-friendlty prime and its $2^{17}$ root of unity will calculated according to method described in the post[^2].

The unit test are provided.

### Mod ops

For robust execution and ease of development a special trait for modular arithmetics was introduced - `ModMul`.
```
pub trait ModOps: Add<Output = Self> + Mul<Output = Self> + Sub<Output = Self> + Div<Output = Self> + Copy + Sized {
    fn mod_add(self, rhs: Self, modulus: Self) -> Self;
    fn mod_mul(self, rhs: Self, modulus: Self) -> Self;
    fn mod_sub(self, rhs: Self, modulus: Self) -> Self;
    fn mod_exp(self, exp: Self, modulus: Self) -> Self;
    fn mod_inv(self, modulus: Self) -> Self;
}
```
Its trait bounds are restricful enough to ensure that only primitive types can implement it. The trait was implemented for **u64** for demonstration purposes.

The **u64** version of the code used following algorithms for more efficient operations.
* `mod_mul` used *Russian peasant multiplication*.
* `mod_exp` used *Exponentiation by squaring*.
* `mod_inv` used *Fermat's little theorem*.


## Future plans

There are two future optimizations are planned.
1. Use instruction sets from **Intel CPUs** with **AVX2** or **AVX512** extensions for speeding up mathicatical calculations. These will utilize *unsafe* Rust.
2. Add implementation optimazied for a Goldilock prime: `0xffffffff00000001` = $2^{64} - 2^{32} + 1$. It has $2^17$-th root of unity equal to `0xabd0a6e8aa3d8a0e`. This prime has a reduction process that has no branches in assembly code. 

## References
[^1]: [Number Theoretic Transform and Its Applications in Lattice-based Cryptosystems: A Survey](https://arxiv.org/pdf/2211.13546)
[^2]: [Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography](https://eprint.iacr.org/2016/504)
[^3]: [Finding the n-th root of unity in a finite field](https://crypto.stackexchange.com/a/63616)
[^4]: [A Complete Beginner Guide to the Number Theoretic Transform (NTT)](https://eprint.iacr.org/2024/585.pdf)