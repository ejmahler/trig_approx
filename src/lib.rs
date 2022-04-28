use polynomial::OddPolynomial;
use roots;
use optimize::{self, NelderMeadBuilder, Minimizer};
use ndarray::{self, ArrayView1, Array};
use float_next_after::NextAfter;

mod polynomial;

const HALF_PI : f64 = std::f64::consts::PI * 0.5;
const HALF_PI_SQUARED : f64 = HALF_PI * HALF_PI;

pub struct ApproximationResult<const D: usize> {
    pub minimized_absolute_polynomial : OddPolynomial<D>,
    pub minimized_absolute_error : f64,
    pub minimized_relative_polynomial : OddPolynomial<D>,
    pub minimized_relative_error : f64,
}

fn iter_float_range(min: f64, max: f64, steps: usize) -> impl Iterator<Item=f64> {
    let width = max - min;
    let divisor = (steps - 1) as f64;
    (0..steps).map(move |i| (i as f64 / divisor) * width + min)
}

fn compute_maximum_error<ErrorFn: Fn(f64) -> f64, ErrorDerivativeFn: Copy + Fn(f64) -> f64>(min: f64, max: f64, error_fn: ErrorFn, derivative_fn: ErrorDerivativeFn) -> f64 {
    let mut prev_sample = min;
    let mut prev_derivative = derivative_fn(min);
    let mut largest_absolute_error = error_fn(min).abs();
    for sample in iter_float_range(min, max, 50).skip(1) {
        let new_derivative = derivative_fn(sample);

        if new_derivative.signum() != prev_derivative.signum() {
            let solved_sample = roots::find_root_regula_falsi(prev_sample, sample, derivative_fn, &mut 1e-15f64);
            let error_at_solved = error_fn(solved_sample.unwrap()).abs();
            largest_absolute_error = largest_absolute_error.max(error_at_solved);
        }
        prev_sample = sample;
        prev_derivative = new_derivative;

    }
    largest_absolute_error
}

fn compute_maximum_absolute_error<const D: usize>(min: f64, max: f64, sin_fn: OddPolynomial<D>) -> f64 {
    let sin_derivative_fn = sin_fn.derivative();
    let absolute_error_fn = |x| {
        sin_fn.eval(x) - x.sin()
    };
    let absolute_error_derivative_fn = |x| {
        sin_derivative_fn.eval(x) - x.cos()
    };
    compute_maximum_error(min, max, absolute_error_fn, absolute_error_derivative_fn)
}
fn compute_maximum_relative_error<const D: usize>(min: f64, max: f64, sin_fn: OddPolynomial<D>) -> f64 {
    let sin_derivative_fn = sin_fn.derivative();
    let relative_error_fn = |x| {
        sin_fn.eval(x) / x.sin() - 1.0
    };
    let relative_error_derivative_fn = |x: f64| {
        let (sin, cos) = x.sin_cos();

        (sin_derivative_fn.eval(x) - sin_fn.eval(x) * cos / sin) / sin
    };
    compute_maximum_error(min, max, relative_error_fn, relative_error_derivative_fn)
}

pub fn approximate_deg5() -> polynomial::OddPolynomial<3> {
    let build_sin_coeffs_deg5 = |a| {
        let raw_polynomial = polynomial::EvenPolynomial([
            -1.0 / HALF_PI_SQUARED,
            a,
            0.0
        ]);
        let derivative_coeffs = raw_polynomial.polynomial_multiply_plusminus(HALF_PI);
        let coeffs = derivative_coeffs.integral();
        coeffs
    };

    let deg5_sample_fn = |sample| {
        let poly = build_sin_coeffs_deg5(sample);
        let result = poly.eval(HALF_PI);
        result - 1.0
    };

    let solved_a = roots::find_root_regula_falsi(-1.0, 1.0, deg5_sample_fn, &mut 1e-50f64).unwrap();
    let solved_coeffs = build_sin_coeffs_deg5(solved_a);
    solved_coeffs
}

fn build_sin_coeffs<const C: usize, const P: usize>(constant_params: [f64; C], dependent_param: f64) -> OddPolynomial<P> {
    let mut raw_polynomial = polynomial::EvenPolynomial([0.0; P]);
    raw_polynomial.0[0] = -1.0 / HALF_PI_SQUARED;
    for i in 0..C {
        raw_polynomial.0[i+1] = constant_params[i];
    }
    raw_polynomial.0[raw_polynomial.0.len() - 2] = dependent_param;

    let derivative_coeffs = raw_polynomial.polynomial_multiply_plusminus(HALF_PI);
    let coeffs = derivative_coeffs.integral();
    coeffs
}

fn solve_final_coefficient<const C: usize, const P: usize>(params: ArrayView1<f64>) -> OddPolynomial<P> {
    let constants = {
        let mut constants_builder = [0.0; C];
        for i in 0..C {
            constants_builder[i] = params[i];
        }
        constants_builder
    };
    let sample_fn = |dependent| {
        let poly = build_sin_coeffs::<C, P>(constants, dependent);
        let result = poly.eval(HALF_PI);
        result - 1.0
    };

    let mut solved_dependent = roots::find_root_regula_falsi(-1.0, 1.0, sample_fn, &mut 1e-50f64).unwrap();
    let mut coeffs = build_sin_coeffs(constants, solved_dependent);
    while coeffs.eval(HALF_PI) < 1.0 {
        solved_dependent = solved_dependent.next_after(std::f64::INFINITY);
        coeffs = build_sin_coeffs(constants, solved_dependent);
    }
    while coeffs.eval(HALF_PI) > 1.0 {
        solved_dependent = solved_dependent.next_after(std::f64::NEG_INFINITY);
        coeffs = build_sin_coeffs(constants, solved_dependent);
    }
    coeffs
}

pub fn approximate_sin<const C: usize, const P: usize>() -> ApproximationResult<P> {
    let compute_absolute_error = |params: ArrayView1<f64>| {
        compute_maximum_absolute_error(0.0, HALF_PI, solve_final_coefficient::<C, P>(params))
    };
    let compute_relative_error = |params: ArrayView1<f64>| {
        compute_maximum_relative_error(0.0001, HALF_PI, solve_final_coefficient::<C, P>(params))
    };

    let minimizer = NelderMeadBuilder::default()
        .xtol(1e-20f64)
        .ftol(1e-20f64)
        .maxiter(5000000) 
        .build()
        .unwrap();
    let args = Array::from_vec(vec![0.0; C]);
    let result_absolute = minimizer.minimize(compute_absolute_error, args.view());
    let result_relative = minimizer.minimize(compute_relative_error, result_absolute.view());
    let result_absolute = minimizer.minimize(compute_absolute_error, result_relative.view());

    ApproximationResult {
        minimized_absolute_polynomial: solve_final_coefficient::<C, P>(result_absolute.view()),
        minimized_absolute_error: compute_absolute_error(result_absolute.view()),
        minimized_relative_polynomial: solve_final_coefficient::<C, P>(result_relative.view()),
        minimized_relative_error: compute_relative_error(result_relative.view()),
    }
}
pub fn approximate_deg7() -> ApproximationResult<4> {
    approximate_sin::<1, 4>()
}
pub fn approximate_deg9() -> ApproximationResult<5> {
    approximate_sin::<2, 5>()
}
pub fn approximate_deg11() -> ApproximationResult<6> {
    approximate_sin::<3, 6>()
}
pub fn approximate_deg13() -> ApproximationResult<7> {
    approximate_sin::<4, 7>()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn print_results<const D: usize>(result: ApproximationResult<D>) {
        println!("Degree {} polynomial:", result.minimized_absolute_polynomial.get_degree());
        println!("Minimized absolute error: [max error: {}, {}]", result.minimized_absolute_error, result.minimized_absolute_polynomial.to_string());
        println!("Minimized relative error: [max error: {}, {}]", result.minimized_relative_error, result.minimized_relative_polynomial.to_string());
    }

    #[test]
    fn test_approximation() {
        let coeffs_deg5 = approximate_deg5();
        println!("Degree 5 polynomial:\n{}", coeffs_deg5.to_string());
        println!("Max absolute error: {}", compute_maximum_absolute_error( 0.0, HALF_PI, coeffs_deg5));
        println!("Max relative error: {}", compute_maximum_relative_error( 0.0001, HALF_PI, coeffs_deg5));

        
        print_results(approximate_deg7());
        print_results(approximate_deg9());
        print_results(approximate_deg11());
        print_results(approximate_deg13());
    }
}
