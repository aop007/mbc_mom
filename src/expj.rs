use num_complex::Complex64;
use std::f64::consts::PI;

const EULER_GAMMA: f64 = 0.5772156649015328;

/// Evaluates the exponential integral E_1(z) for a complex argument z.
pub fn c_e1(z: Complex64) -> Complex64 {
    let ab = z.norm();
    
    // Handle singularity
    if ab == 0.0 {
        return Complex64::new(f64::INFINITY, 0.0);
    }

    // Use Richmond's Gauss-Laguerre rational approximation for large fields
    if (z.re >= 0.0 && ab > 10.0) || (z.re < 0.0 && z.im.abs() > 10.0) {
        return richmond_large_e1(z);
    }
    
    // Default to power series expansion near the origin
    power_series_e1(z)
}

fn power_series_e1(z: Complex64) -> Complex64 {
    let mut sum = Complex64::new(0.0, 0.0);
    let mut term = Complex64::new(1.0, 0.0); // T_0 = 1

    // Taylor series: \sum_{n=1}^\infty \frac{(-z)^n}{n \cdot n!}
    // Recurrence: T_n = T_{n-1} * (-z) / n
    for n in 1..100 {
        let n_f64 = n as f64;
        term = term * (-z) / n_f64;
        let delta = term / n_f64;
        sum += delta;
        
        if delta.norm() < 1e-15 {
            break;
        }
    }

    -EULER_GAMMA - z.ln() - sum
}

fn richmond_large_e1(z: Complex64) -> Complex64 {
    // Gauss-Laguerre quadrature constants extracted directly from 
    // J.H. Richmond's 1974 Fortran EXPJ subroutine.
    let w = [
        0.409319, 0.421831, 0.147126, 0.0206335,
        0.00107401, 0.0000158654, 0.0000000317031
    ];
    let x = [
        0.193044, 1.02666, 2.56788, 4.90035,
        8.18215, 12.7342, 19.3957
    ];
    
    let mut sum = Complex64::new(0.0, 0.0);
    for i in 0..7 {
        sum += w[i] / (z + x[i]);
    }
    sum * (-z).exp()
}

/// Evaluates the path integral tracking negative real branch crossings.
/// W12 = E_1(v1) - E_1(v2) + j * 2n * PI
pub fn expj_path(v1: Complex64, v2: Complex64) -> Complex64 {
    let e1 = c_e1(v1);
    let e2 = c_e1(v2);
    
    // Richmond's branch cut tracking logic
    let z = v2 / v1;
    let mut th = z.arg() - v2.arg() + v1.arg();
    
    let ab = th.abs();
    if ab < 1.0 {
        th = 0.0;
    } else if th > 1.0 {
        th = 2.0 * PI;
    } else if th < -1.0 {
        th = -2.0 * PI;
    }
    
    e1 - e2 + Complex64::new(0.0, th)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_e1_real_positive() {
        // E_1(1.0) is approx 0.2193839
        let z = Complex64::new(1.0, 0.0);
        let res = c_e1(z);
        assert!((res.re - 0.2193839).abs() < 1e-5);
        assert!((res.im - 0.0).abs() < 1e-5);
    }

    #[test]
    fn test_expj_path_no_crossing() {
        let v1 = Complex64::new(1.0, 1.0);
        let v2 = Complex64::new(2.0, 2.0);
        let res = expj_path(v1, v2);
        let expected = c_e1(v1) - c_e1(v2);
        assert!((res - expected).norm() < 1e-10);
    }

    #[test]
    fn test_expj_path_crossing() {
        // Path from quadrant 2 to quadrant 3 crosses the negative real axis
        let v1 = Complex64::new(-1.0, 1.0);
        let v2 = Complex64::new(-1.0, -1.0);
        let res = expj_path(v1, v2);
        
        // Ensure the imaginary jump of roughly 2*PI occurred
        assert!(res.im.abs() > PI); 
    }
}