use core::arch::x86_64::*;

#[derive(Copy, Clone, Debug)]
pub struct OddPolynomial<const D: usize>(pub [f64; D]);

#[derive(Copy, Clone, Debug)]
pub struct EvenPolynomial<const D: usize>(pub [f64; D]);

#[derive(Copy, Clone, Debug)]
pub struct Polynomial<const D: usize>(pub [f64; D]);

impl<const D: usize> OddPolynomial<D> {
    pub fn eval(&self, x: f64) -> f64 {
        unsafe { self.eval_unsafe(x) }
    }

    #[target_feature(enable = "fma")]
    pub unsafe fn eval_unsafe(&self, x: f64) -> f64 {
        let x_vec = _mm_set_sd(x);
        let xsquared = _mm_mul_sd(x_vec, x_vec);
        let mut result = _mm_set_sd(self.0[self.0.len() - 1]);

        for coeff in self.0[0..(self.0.len() - 1)].iter().rev() {
            let coeff_vec = _mm_set_sd(*coeff);
            result = _mm_fmadd_sd(result, xsquared, coeff_vec);
        }

        result = _mm_mul_sd(result, x_vec);
        
        let mut result_scalar = 0.0;
        _mm_store_sd(&mut result_scalar, result);
        result_scalar
    }
    pub fn derivative(&self) -> EvenPolynomial<D> {
        let mut result = self.0;
        for (i, e) in result.iter_mut().enumerate() {
            let old_exponent = (i * 2 + 1) as f64;
            *e = *e * old_exponent;
        }
        EvenPolynomial(result)
    }
    pub fn get_degree(&self) -> usize {
        D * 2 - 1
    }
    pub fn to_string(&self) -> String {
        let mut chunks = vec![];
        for i in (0..self.0.len()).rev() {
            let exponent = i * 2 + 1;
            if exponent == 1 {
                chunks.push(format!("{}x", self.0[i]));
            } else {
                chunks.push(format!("{}x^{}", self.0[i], exponent));
            }
        }
        chunks.join(" + ")
    }
}

impl<const D: usize> EvenPolynomial<D> {
    pub fn eval(&self, x: f64) -> f64 {
        unsafe { self.eval_unsafe(x) }
    }

    #[target_feature(enable = "fma")]
    pub unsafe fn eval_unsafe(&self, x: f64) -> f64 {
        let x_vec = _mm_set_sd(x);
        let xsquared = _mm_mul_sd(x_vec, x_vec);
        let mut result = _mm_set_sd(self.0[self.0.len() - 1]);

        for coeff in self.0[0..(self.0.len() - 1)].iter().rev() {
            let coeff_vec = _mm_set_sd(*coeff);
            result = _mm_fmadd_sd(result, xsquared, coeff_vec);
        }
        
        let mut result_scalar = 0.0;
        _mm_store_sd(&mut result_scalar, result);
        result_scalar
    }
    
    pub fn integral(&self) -> OddPolynomial<D> {
        let mut result = self.0;
        for (i, e) in result.iter_mut().enumerate() {
            let new_exponent = (i * 2 + 1) as f64;
            *e = *e / new_exponent;
        }
        OddPolynomial(result)
    }

    // Multiplies this polynomial by (x + c) * (x - c). The most-significant coefficeint of this polynomial should be zero, to make room for the multiplication.
    pub fn polynomial_multiply_plusminus(&self, c: f64) -> EvenPolynomial<D> {
        let mut result = [0.0; D];

        let neg_c_squared = -(c * c);
        for i in (1..result.len()).rev() {
            result[i] = self.0[i - 1] + self.0[i] * neg_c_squared;
        }
        result[0] = self.0[0] * neg_c_squared;

        EvenPolynomial(result)
    }
}

macro_rules! impl_full_polynomial {
    ($even_degree: expr, $full_degree: expr) => {
        impl EvenPolynomial<$even_degree> {
            #[allow(dead_code)]
            pub fn to_full_polynomial(&self) -> Polynomial<$full_degree> {
                let mut result = [0.0; $full_degree];
                for i in 0..$even_degree {
                    result[i*2] = self.0[i];
                }
                Polynomial(result)
            }
        }
    };
}
impl_full_polynomial!(1, 1);
impl_full_polynomial!(2, 3);
impl_full_polynomial!(3, 5);
impl_full_polynomial!(4, 7);
impl_full_polynomial!(5, 9);
impl_full_polynomial!(6, 11);
impl_full_polynomial!(7, 13);
impl_full_polynomial!(8, 15);
impl_full_polynomial!(9, 17);


#[cfg(test)]
mod tests {
    use super::*;

    #[track_caller]
    fn assert_nearly_equal(expected: f64, observed: f64) {
        assert!((expected - observed).abs() < 0.001);
    }

    #[test]
    fn test_polynomial_multiply() {
        // Test a "polynomial" of zero
        for c in 0..10 {
            let poly = EvenPolynomial([0.0]);
            let multiplied = poly.polynomial_multiply_plusminus(c as f64);
            assert_nearly_equal(0.0, multiplied.0[0]);
        }

        // Test with a "polynomial" of just a constant
        for i in 0..3
        {
            let poly = EvenPolynomial([i as f64, 0.0]);
            for c in 0..4 {
                let multiplied = poly.polynomial_multiply_plusminus(c as f64);
                assert_nearly_equal(-(c * c) as f64 * i as f64, multiplied.0[0]);
                assert_nearly_equal(i as f64, multiplied.0[1]);
            }
        }

        struct TestData<const D: usize> {
            polynomial: EvenPolynomial<D>,
            c: f64,
            expected: EvenPolynomial<D>,
        }

        // Test with a quadratic polynomial
        let quadratic_tests = vec![
            TestData { polynomial: EvenPolynomial([0.0, 0.0, 0.0]), c: 0.0, expected: EvenPolynomial([0.0, 0.0, 0.0]) },
            TestData { polynomial: EvenPolynomial([1.0, 0.0, 0.0]), c: 0.0, expected: EvenPolynomial([0.0, 1.0, 0.0]) },
            TestData { polynomial: EvenPolynomial([3.0, 2.0, 0.0]), c: 0.0, expected: EvenPolynomial([0.0, 3.0, 2.0]) },
            TestData { polynomial: EvenPolynomial([0.0, 0.0, 0.0]), c: 0.0, expected: EvenPolynomial([0.0, 0.0, 0.0]) },
            TestData { polynomial: EvenPolynomial([0.0, 0.0, 0.0]), c: 1.0, expected: EvenPolynomial([0.0, 0.0, 0.0]) },
            TestData { polynomial: EvenPolynomial([0.0, 0.0, 0.0]), c: 2.0, expected: EvenPolynomial([0.0, 0.0, 0.0]) },

            TestData { polynomial: EvenPolynomial([1.0, 0.0, 0.0]), c: 1.5, expected: EvenPolynomial([-2.25, 1.0, 0.0]) },
            TestData { polynomial: EvenPolynomial([1.0, 0.0, 0.0]), c: -1.5, expected: EvenPolynomial([-2.25, 1.0, 0.0]) },

            TestData { polynomial: EvenPolynomial([2.0, 3.0, 0.0]), c: -1.5, expected: EvenPolynomial([-4.5, -4.75, 3.0]) },
            TestData { polynomial: EvenPolynomial([2.0, 3.0, 0.0]), c: 1.5, expected: EvenPolynomial([-4.5, -4.75, 3.0]) },

            TestData { polynomial: EvenPolynomial([-2.0, 3.0, 0.0]), c: 1.5, expected: EvenPolynomial([4.5, -8.75, 3.0]) },
            TestData { polynomial: EvenPolynomial([-2.0, 3.0, 0.0]), c: -1.5, expected: EvenPolynomial([4.5, -8.75, 3.0]) },

            TestData { polynomial: EvenPolynomial([1.0, 0.0, 0.0]), c: 10.5, expected: EvenPolynomial([-110.25, 1.0, 0.0]) },
            TestData { polynomial: EvenPolynomial([1.0, 0.0, 0.0]), c: -10.5, expected: EvenPolynomial([-110.25, 1.0, 0.0]) },

            TestData { polynomial: EvenPolynomial([2.0, 3.0, 0.0]), c: 10.5, expected: EvenPolynomial([-220.5, -328.75, 3.0]) },
            TestData { polynomial: EvenPolynomial([2.0, 3.0, 0.0]), c: -10.5, expected: EvenPolynomial([-220.5, -328.75, 3.0]) },

            TestData { polynomial: EvenPolynomial([-2.0, 3.0, 0.0]), c: 10.5, expected: EvenPolynomial([220.5, -332.75, 3.0]) },
            TestData { polynomial: EvenPolynomial([-2.0, 3.0, 0.0]), c: -10.5, expected: EvenPolynomial([220.5, -332.75, 3.0]) },
        ];
        for test in quadratic_tests {
            let multiplied = test.polynomial.polynomial_multiply_plusminus(test.c);
            for i in 0..test.expected.0.len() {
                assert_nearly_equal(test.expected.0[i], multiplied.0[i]);
            }
        }
    }

    #[test]
    fn test_polynomial_calculus() {

        // Test the even polynomial's integral function 
        struct IntegralTestData<const D: usize> {
            polynomial: EvenPolynomial<D>,
            expected: OddPolynomial<D>,
        }

        let integral_tests = vec![
            IntegralTestData { polynomial: EvenPolynomial([0.0, 0.0, 0.0]), expected: OddPolynomial([0.0, 0.0, 0.0]) },
            IntegralTestData { polynomial: EvenPolynomial([0.0, 0.0, 1.0]), expected: OddPolynomial([0.0, 0.0, 0.2]) },
            IntegralTestData { polynomial: EvenPolynomial([0.0, 1.0, 0.0]), expected: OddPolynomial([0.0, 1.0 / 3.0, 0.0]) },
            IntegralTestData { polynomial: EvenPolynomial([1.0, 0.0, 0.0]), expected: OddPolynomial([1.0, 0.0, 0.0]) },
            IntegralTestData { polynomial: EvenPolynomial([3.0, 6.0, 5.0]), expected: OddPolynomial([3.0, 2.0, 1.0]) },
        ];
        for test in integral_tests {
            let integral = test.polynomial.integral();
            println!("polynomial: {:?}, integral: {:?}, expectd: {:?}", test.polynomial, integral, test.expected);
            for i in 0..test.expected.0.len() {
                assert_nearly_equal(test.expected.0[i], integral.0[i]);
            }
        }

        // Test the odd polynomial's derivative function 
        struct DerivativeTestData<const D: usize> {
            polynomial: OddPolynomial<D>,
            expected: EvenPolynomial<D>,
        }

        let integral_tests = vec![
            DerivativeTestData { polynomial: OddPolynomial([0.0, 0.0, 0.0]), expected: EvenPolynomial([0.0, 0.0, 0.0]) },
            DerivativeTestData { polynomial: OddPolynomial([0.0, 0.0, 1.0]), expected: EvenPolynomial([0.0, 0.0, 5.0]) },
            DerivativeTestData { polynomial: OddPolynomial([0.0, 1.0, 0.0]), expected: EvenPolynomial([0.0, 3.0, 0.0]) },
            DerivativeTestData { polynomial: OddPolynomial([1.0, 0.0, 0.0]), expected: EvenPolynomial([1.0, 0.0, 0.0]) },
            DerivativeTestData { polynomial: OddPolynomial([3.0, 6.0, 5.0]), expected: EvenPolynomial([3.0, 18.0, 25.0]) },
        ];
        for test in integral_tests {
            let derivative = test.polynomial.derivative();
            for i in 0..test.expected.0.len() {
                assert_nearly_equal(test.expected.0[i], derivative.0[i]);
            }
        }

        // Test that derivative and integral are opposites of each other
        let inverse_tests = vec![
            EvenPolynomial([0.0, 0.0, 0.0, 0.0]),
            EvenPolynomial([1.0, 0.0, 0.0, 0.0]),
            EvenPolynomial([0.0, 1.0, 0.0, 0.0]),
            EvenPolynomial([0.0, 0.0, 1.0, 0.0]),
            EvenPolynomial([0.0, 0.0, 0.0, 1.0]),
            EvenPolynomial([234.0, 25.0, -345.0, 1.0]),
            EvenPolynomial([0.001, 123.1, 1121.7, 0.0]),
            EvenPolynomial([-100.5, 1.1, 1.1, 2.2]),
        ];
        for test in inverse_tests {
            let integral = test.integral();
            let derivative = integral.derivative();
            for i in 0..test.0.len() {
                assert_nearly_equal(test.0[i], derivative.0[i]);
            }
        }
    }
}
