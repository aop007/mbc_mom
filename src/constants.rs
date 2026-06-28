use std::f64::consts::PI;
use num_complex::Complex64;

pub const C: f64 = 299_792_458.0_f64;
pub const MU_0: f64 = 4.0e-7 * PI;
pub const EPS_0: f64 = 1.0 / (MU_0 / (C * C));
pub const ETA: f64 = MU_0 * C;
pub const J: Complex64 = Complex64::new(0.0, 1.0);