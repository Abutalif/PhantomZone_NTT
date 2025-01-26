fn main() {
    // Special prime
    let spec_q:u64 = ((((1 << 63) - (1 << 31)))<<1) + 1;
    let spec_bit = 0xffffffff00000001u64;

    assert_eq!(spec_q, spec_bit);

    println!("Special Q: {spec_q}; special Bit: {spec_bit}");
}