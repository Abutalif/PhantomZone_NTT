use app::dft::ntt::NttBuilder;

fn main() {
    let q = 1u64;
    let n = 16;
    let ntt = NttBuilder::new(q, n).build().expect("Cannot build");
}