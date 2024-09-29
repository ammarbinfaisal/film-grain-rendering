use std::f32::consts::PI;

#[derive(Debug, Clone, Copy)]
pub struct NoisePrng {
    state: u32,
}

impl NoisePrng {
    pub fn new() -> Self {
        Self { state: 0 }
    }
}

fn wang_hash(mut seed: u32) -> u32 {
    seed = (seed ^ 61) ^ (seed.wrapping_shr(16));
    seed = seed.wrapping_mul(9);
    seed = seed ^ (seed.wrapping_shr(4));
    seed = seed.wrapping_mul(668265261);
    seed = seed ^ (seed.wrapping_shr(15));
    seed
}

pub fn cellseed(x: u32, y: u32, offset: u32) -> u32 {
    const PERIOD: u32 = 65536; // 2^16
    let mut s = ((y % PERIOD) * PERIOD + (x % PERIOD)) + offset;
    if s == 0 {
        s = 1;
    }
    s
}

pub fn mysrand(p: &mut NoisePrng, seed: u32) {
    p.state = wang_hash(seed);
}

pub fn myrand(p: &mut NoisePrng) -> u32 {
    p.state ^= p.state.wrapping_shl(13);
    p.state ^= p.state.wrapping_shr(17);
    p.state ^= p.state.wrapping_shl(5);
    p.state
}

pub fn myrand_uniform_0_1(p: &mut NoisePrng) -> f32 {
    let a = myrand(p) as f32 / u32::MAX as f32;
    return a;
}

pub fn myrand_gaussian_0_1(p: &mut NoisePrng) -> f32 {
    let u = myrand_uniform_0_1(p);
    let v = myrand_uniform_0_1(p);
    f32::sqrt(-2.0 * f32::ln(u)) * f32::cos(2.0 * PI * v)
}

pub fn my_rand_poisson(p: &mut NoisePrng, lambda: f32, prod_in: f32) -> u32 {
    let u = myrand_uniform_0_1(p);
    let mut x = 0;
    let mut prod = if prod_in <= 0.0 {
        f32::exp(-lambda)
    } else {
        prod_in
    };

    let mut sum = prod;
    while u > sum && x < (10000.0 * lambda).floor() as u32 {
        x += 1;
        prod *= lambda / x as f32;
        sum += prod;
    }

    x
}
