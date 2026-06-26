use num_complex::Complex64;
use std::f64::consts::PI;
use rayon::prelude::*;


/// Holds the pre-calculated physical constants for the ground environment
#[derive(Clone, Copy, Debug)]
pub struct GroundPhysics {
    pub freq_hz: f64,
    pub k0: f64,              // Wavenumber in free space
    pub eps_c: Complex64,     // Complex relative permittivity of the ground
    pub k1: Complex64,        // Wavenumber in the ground
    pub pole: Complex64,      // The Sommerfeld Pole location (lambda_p)
}

impl GroundPhysics {
    pub fn new(freq_hz: f64, eps_r: f64, sigma: f64) -> Self {
        let c = 299_792_458.0;
        let eps_0 = 8.8541878128e-12;
        let omega = 2.0 * PI * freq_hz;
        
        let k0 = omega / c;
        let eps_c = Complex64::new(eps_r, -sigma / (omega * eps_0));
        
        // k1 = k0 * sqrt(eps_c)
        let k1 = Complex64::new(k0, 0.0) * eps_c.sqrt();
        
        // Locate the Sommerfeld Pole: lambda_p = k0 * sqrt(eps_c / (eps_c + 1))
        let eps_plus_1 = eps_c + Complex64::new(1.0, 0.0);
        let pole = Complex64::new(k0, 0.0) * (eps_c / eps_plus_1).sqrt();
        
        Self { freq_hz, k0, eps_c, k1, pole }
    }

    /// Evaluates the integration variable lambda along an elliptical contour.
    /// 't' parameter goes from 0.0 to 1.0 (the length of the contour).
    /// The contour dips under the real axis to avoid the pole.
    pub fn contour_path(&self, t: f64, t_max: f64) -> (Complex64, Complex64) {
        // We will stretch the contour slightly past k0
        let a = self.k0 * 1.2; 
        
        // The depth of the contour depends on how close the pole is to the real axis
        // We ensure we dip at least far enough to clear it safely
        let pole_imag = self.pole.im.abs();
        let b = (self.k0 * 0.2).max(pole_imag * 2.0); 
        
        // Map t into an angle from 0 to PI
        let theta = (t / t_max) * PI;
        
        // Parametric ellipse equation: x = a/2 * (1 - cos(theta)), y = -b * sin(theta)
        let lambda_real = (a / 2.0) * (1.0 - theta.cos());
        let lambda_imag = -b * theta.sin();
        
        let lambda = Complex64::new(lambda_real, lambda_imag);
        
        // We also need the derivative of lambda with respect to t (d_lambda / dt)
        // to correctly scale the numerical integration steps
        let d_real_dt = (a / 2.0) * (PI / t_max) * theta.sin();
        let d_imag_dt = -b * (PI / t_max) * theta.cos();
        let d_lambda = Complex64::new(d_real_dt, d_imag_dt);
        
        (lambda, d_lambda)
    }

    /// Evaluates the TM (Vertical) Sommerfeld integral for a specific radial distance (rho) 
    /// and combined height (z_sum = z + z').
    pub fn evaluate_tm(&self, rho: f64, z_sum: f64) -> Complex64 {
        let mut sum = Complex64::new(0.0, 0.0);
        
        // 1. Define the Integrand Function F(lambda)
        let integrand = |lambda: Complex64| -> Complex64 {
            let l2 = lambda * lambda;
            let k0_2 = Complex64::new(self.k0 * self.k0, 0.0);
            
            // u0 = sqrt(lambda^2 - k0^2)
            let mut u0 = (l2 - k0_2).sqrt();
            // Enforce Radiation Condition (Real part must be positive)
            if u0.re < 0.0 { u0 = -u0; } 
            
            // u1 = sqrt(lambda^2 - k1^2)
            let k1_2 = self.k1 * self.k1;
            let mut u1 = (l2 - k1_2).sqrt();
            if u1.re < 0.0 { u1 = -u1; }
            
            // Fresnel TM Reflection Coefficient
            let r_v = (self.eps_c * u0 - u1) / (self.eps_c * u0 + u1);
            
            // Exponential decay away from the boundary
            let exp_term = (-u0 * z_sum).exp();
            
            // Combine the fundamental physics terms
            let term = r_v * (lambda / u0) * exp_term;
            
            // Multiply by our fast Abramowitz & Stegun Bessel function.
            // Since our contour is very shallow, lambda.re accurately dominates the oscillation.
            let bessel = crate::sommerfeld::j0(lambda.re * rho); 
            
            term * Complex64::new(bessel, 0.0)
        };

        // 2. Part 1: Contour Integration (Simpson's Rule)
        // We evaluate along our safe elliptical path to dodge the Sommerfeld Pole
        let contour_steps = 100;
        let t_max = 1.0;
        let dt = t_max / (contour_steps as f64);
        
        for i in 0..contour_steps {
            let t1 = (i as f64) * dt;
            let t2 = (i as f64 + 0.5) * dt;
            let t3 = (i as f64 + 1.0) * dt;
            
            let (l1, dl1) = self.contour_path(t1, t_max);
            let (l2, dl2) = self.contour_path(t2, t_max);
            let (l3, dl3) = self.contour_path(t3, t_max);
            
            let val = (integrand(l1)*dl1 + integrand(l2)*dl2*4.0 + integrand(l3)*dl3) * (dt / 6.0);
            sum += val;
        }

        // 3. Part 2: Real Axis Tail Integration
        // Start exactly where the contour left off (lambda = 1.2 k0)
        let mut lambda_real = self.k0 * 1.2;
        let step_size = self.k0 * 0.05; 
        
        // March towards infinity until exponential decay kills the integrand
        for _ in 0..1000 {
            let l_cplx = Complex64::new(lambda_real, 0.0);
            
            // Basic trapezoidal numerical step
            let val = integrand(l_cplx) * step_size;
            sum += val;
            
            lambda_real += step_size;
            
            // Convergence Check: If the contribution drops below 10^-8, we are done
            if val.norm() < 1e-8 {
                break;
            }
        }
        
        sum
    }
}

/// Computes the Bessel function of the first kind, order zero, J0(x).
/// Uses the optimized Abramowitz & Stegun polynomial approximations.
pub fn j0(x: f64) -> f64 {
    let ax = x.abs();
    
    // Domain 1: Small arguments (Highly accurate polynomial series)
    if ax <= 3.0 {
        let y = (x / 3.0) * (x / 3.0);
        return 1.0 
            - 2.2499997 * y 
            + 1.2656208 * y.powi(2) 
            - 0.3163866 * y.powi(3) 
            + 0.0444479 * y.powi(4) 
            - 0.0039444 * y.powi(5) 
            + 0.0002100 * y.powi(6);
    } 
    // Domain 2: Large arguments (Asymptotic expansion with phase shift)
    else {
        let y = 3.0 / ax;
        let f0 = 0.79788456 
            - 0.00000077 * y 
            - 0.00552740 * y.powi(2) 
            - 0.00009512 * y.powi(3) 
            + 0.00137237 * y.powi(4) 
            - 0.00072805 * y.powi(5) 
            + 0.00014476 * y.powi(6);
            
        let theta0 = ax 
            - 0.78539816 
            - 0.04166397 * y 
            - 0.00003954 * y.powi(2) 
            + 0.00262573 * y.powi(3) 
            - 0.00054125 * y.powi(4) 
            - 0.00029333 * y.powi(5) 
            + 0.00013558 * y.powi(6);
            
        return f0 * theta0.cos() / ax.sqrt();
    }
}

/// A highly optimized 2D Lookup Table for Sommerfeld Integrals.
/// Uses logarithmic spacing for both radial distance (rho) and height (z_sum).
#[derive(Clone)]
pub struct SommerfeldTable {
    rho_min: f64,
    rho_max: f64,
    z_min: f64,
    z_max: f64,
    rho_pts: usize,
    z_pts: usize,
    grid: Vec<Complex64>,
}

impl SommerfeldTable {
    /// Pre-computes the Sommerfeld integrals over a 2D logarithmic grid.
    pub fn generate(physics: &GroundPhysics) -> Self {
        let rho_min = 1e-4;   // 0.1 mm (Near-field safety limit)
        let rho_max = 1000.0; // 1 km (Far-field limit)
        let z_min = 1e-4;     // 0.1 mm
        let z_max = 1000.0;
        
        let rho_pts = 100;
        let z_pts = 100;

        let get_val = |min: f64, max: f64, pts: usize, i: usize| -> f64 {
            let log_min = min.ln();
            let log_max = max.ln();
            let frac = i as f64 / (pts - 1) as f64;
            (log_min + frac * (log_max - log_min)).exp()
        };

        let mut grid = vec![Complex64::new(0.0, 0.0); rho_pts * z_pts];

        // Fire up all CPU cores to crunch the 10,000 complex contour integrals!
        grid.par_iter_mut().enumerate().for_each(|(idx, val)| {
            let r_idx = idx / z_pts;
            let z_idx = idx % z_pts;
            
            let rho = get_val(rho_min, rho_max, rho_pts, r_idx);
            let z_sum = get_val(z_min, z_max, z_pts, z_idx);
            
            *val = physics.evaluate_tm(rho, z_sum);
        });

        Self { rho_min, rho_max, z_min, z_max, rho_pts, z_pts, grid }
    }

    /// Performs a fast 2D Bilinear Interpolation on the logarithmic grid.
    pub fn interpolate(&self, mut rho: f64, mut z_sum: f64) -> Complex64 {
        // Enforce safety limits to prevent out-of-bounds indexing or singularity blowups
        rho = rho.clamp(self.rho_min, self.rho_max);
        z_sum = z_sum.clamp(self.z_min, self.z_max);

        // Convert the requested physical coordinates to grid indices
        let r_frac = (rho.ln() - self.rho_min.ln()) / (self.rho_max.ln() - self.rho_min.ln()) * (self.rho_pts - 1) as f64;
        let z_frac = (z_sum.ln() - self.z_min.ln()) / (self.z_max.ln() - self.z_min.ln()) * (self.z_pts - 1) as f64;

        let r0 = (r_frac.floor() as usize).min(self.rho_pts - 2);
        let r1 = r0 + 1;
        let z0 = (z_frac.floor() as usize).min(self.z_pts - 2);
        let z1 = z0 + 1;

        let dr = r_frac - r0 as f64;
        let dz = z_frac - z0 as f64;

        // Retrieve the 4 surrounding points from the flattened 1D array
        let v00 = self.grid[r0 * self.z_pts + z0];
        let v10 = self.grid[r1 * self.z_pts + z0];
        let v01 = self.grid[r0 * self.z_pts + z1];
        let v11 = self.grid[r1 * self.z_pts + z1];

        // Bilinear blend in complex space
        let c0 = v00 * (1.0 - dr) + v10 * dr;
        let c1 = v01 * (1.0 - dr) + v11 * dr;
        
        c0 * (1.0 - dz) + c1 * dz
    }
}