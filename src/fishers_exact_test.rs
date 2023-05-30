// Functions needed to perform Fisher's exact test
use thiserror::Error;

use crate::constants::CLOSE_TO_ZERO;

// Errors occurring during the Fisher's Exact Test
#[derive(Error, Debug, PartialEq)]
pub enum FishersExactError {
    /// inproper inputs to hypergeometric test function
    #[error("inproper inputs provided")]
    BadInputs,

    /// attempt to take the logarithm of 0 (for any base this is infinite)
    #[error("cannot take the log of 0")]
    LogZero,

    /// taylor series is unable to converge
    #[error("Taylor series failed to converge")]
    ConvergeFailure,

    /// infinite gamma function value
    #[error("Gamma function returned infinity")]
    InfiniteGamma,

    /// NaN
    #[error("Result is returning a NaN")]
    NaNReturn,

    /// value out of bounds
    #[error("Value is out of bounds")]
    OutOfBoundsValue,
}

// Handle attempting to return 0
fn try_return_zero(as_log: bool) -> Result<f64, FishersExactError> {
    if as_log {
        Err(FishersExactError::LogZero)
    } else {
        Ok(0.0)
    }
}

// Evaluates the n-term Chebyshev series g at x
// Based on a Fortran subroutine by Fullerton at Los Alamos Scientific Laboratory
//     https://www.netlib.org/slatec/fnlib/csevl.f
fn chebyshev_eval(x: f64, g: &[f64], n: i64) -> Result<f64, FishersExactError> {
    if n < 1 || n > 1000 {
        return Err(FishersExactError::OutOfBoundsValue);
    }

    if x < -1.1 || x > 1.1 {
        return Err(FishersExactError::OutOfBoundsValue);
    }

    let twox = 2.0 * x;
    let mut b0: f64 = 0.0;
    let mut b1: f64 = 0.0;
    let mut b2: f64 = 0.0;

    for i in 1..n+1 {
        b2 = b1;
        b1 = b0;
        b0 = b1*twox - b2 + g[(n - i) as usize];
    }

    Ok(0.5 * (b0 - b2))
}

// Natural log of the gamma correction factor for x >= 10.0
// Based on a Fortran subroutine by Fullerton at Los Alamos Scientific Laboratory
//     https://www.netlib.org/slatec/fnlib/r9lgmc.f
fn lgammacor(x: f64) -> Result<f64, FishersExactError> {
    const ALGMCS: [f64; 6] = [
         0.1666389480451863247205729650822e+0 ,
        -0.1384948176067563840732986059135e-4 ,
         0.9810825646924729426157171547487e-8 ,
        -0.1809129475572494194263306266719e-10,
         0.6221098041892605227126015543416e-13,
        -0.3399615005417721944303330599666e-15
    ];
    const NALGM: i64 = 5;
    const XBIG: f64 = 94906265.62425156;     // 2 ^ 26.5
    const XMAX: f64 = 3.745194030963158e306; // f64::MAX / 48

    if x < 10.0 {
        return Err(FishersExactError::NaNReturn);
    } else if x >= XMAX {
        return Err(FishersExactError::OutOfBoundsValue);
    } else if x < XBIG {
        let tmp = 10.0 / x;
        let ret = match chebyshev_eval(2.0*tmp*tmp - 1.0, &ALGMCS, NALGM) {
            Ok(f) => f,
            Err(err) => { return Err(err); },
        };
        return Ok(ret / x);
    }

    Ok(1.0 / (12.0*x))
}

// Calculate gamma function
// Based on a Fortran subroutine by Fullerton at Los Alamos Scientific Laboratory
//     http://www.netlib.org/slatec/fnlib/gamma.f
fn gammafn(x: f64) -> Result<f64, FishersExactError> {
    const GAMCS: [f64; 23] = [
         0.8571195590989331421920062399942e-2 ,
         0.4415381324841006757191315771652e-2 ,
         0.5685043681599363378632664588789e-1 ,
        -0.4219835396418560501012500186624e-2 ,
         0.1326808181212460220584006796352e-2 ,
        -0.1893024529798880432523947023886e-3 ,
         0.3606925327441245256578082217225e-4 ,
        -0.6056761904460864218485548290365e-5 ,
         0.1055829546302283344731823509093e-5 ,
        -0.1811967365542384048291855891166e-6 ,
         0.3117724964715322277790254593169e-7 ,
        -0.5354219639019687140874081024347e-8 ,
         0.9193275519859588946887786825940e-9 ,
        -0.1577941280288339761767423273953e-9 ,
         0.2707980622934954543266540433089e-10,
        -0.4646818653825730144081661058933e-11,
         0.7973350192007419656460767175359e-12,
        -0.1368078209830916025799499172309e-12,
         0.2347319486563800657233471771688e-13,
        -0.4027432614949066932766570534699e-14,
         0.6910051747372100912138336975257e-15,
        -0.1185584500221992907052387126192e-15,
         0.2034148542496373955201026051932e-16
    ];
    const NGAM: i64 = 22;
    const XMIN: f64 = -170.5674972726612;
    const XMAX: f64 =  171.61447887182298;
    const XSML: f64 = 1.01005016708_f64 * f64::MIN; // 0.01_f64.exp() * f64::MIN
    
    // Check for problem inputs
    if x.abs() < CLOSE_TO_ZERO || (x < 0.0 && (x - x.round()).abs() < CLOSE_TO_ZERO) {
        return Err(FishersExactError::NaNReturn);
    }

    let mut y = x.abs();

    // Calculate gamma(x) for |x| <= 10
    if y <= 10.0 {
        // Find gamma(1+y) FOR 0 <= y < 1
        let mut n = x as i64;
        if x < 0.0 {
            n -= 1;
        }

        y = x - n as f64;
        n -= 1;

        let value = match chebyshev_eval(2.0*y - 1.0, &GAMCS, NGAM) {
            Ok(f) => f + 0.9375,
            Err(err) => return Err(err),
        };

        if n == 0 {
            return Ok(value);
        }

        if n < 0 {
            // gamma(x) for -10 <= x < 1
            if y < XSML {
                return Err(FishersExactError::InfiniteGamma);
            }

            return Ok( (0..-n).fold(value, |acc, e| acc / (x+e as f64)) );
        } else {
            // gamma(x) for 2 <= x <= 10
            return Ok( (1..n+1).fold(value, |acc, e| acc * (y+e as f64)) );
        }
    } else {
        // gamma(x) for |x| > 10
        if x > XMAX {
            return Err(FishersExactError::InfiniteGamma);
        }
        if x < XMIN {
            return Ok(0.0);
        }

        let value = if y <= 50.0 && (y - y.trunc()).abs() < CLOSE_TO_ZERO {
            (2..y as i64).fold(1.0, |acc, e| acc * e as f64)
        } else {
            let tmp = if (2.0*y - (2.0*y).trunc()).abs() < CLOSE_TO_ZERO {
                match stirling_error(y) {
                    Ok(f) => f,
                    Err(err) => { return Err(err); },
                }
            } else {
                match lgammacor(y) {
                    Ok(f) => f,
                    Err(err) => { return Err(err); }
                }
            };

            ((y-0.5) * y.ln() - y + std::f64::consts::TAU.sqrt().ln() + tmp).exp()
        };

        if x > 0.0 {
            return Ok(value);
        }

        let sinpiy = (y*std::f64::consts::PI).sin().abs();
        if sinpiy < CLOSE_TO_ZERO {
            return Err(FishersExactError::InfiniteGamma);
        }

        return Ok(-std::f64::consts::PI / (y * sinpiy * value));
    }
}

// Compute ln( |gamma(x)| )
// Based on a Fortran subroutine by Fullerton at Los Alamos Scientific Laboratory
//     https://www.netlib.org/slatec/fnlib/alngam.f
fn lgammafn(x: f64) -> Result<f64, FishersExactError> {
    const XMAX: f64  = 2.5327372760800758e+305; // f64::MAX / f64::MAX.ln()

    // Catch problem inputs
    if x <= 0.0 && (x-x.trunc()).abs() < CLOSE_TO_ZERO {
        return Err(FishersExactError::InfiniteGamma);
    }

    let y = x.abs();

    if y < 1e-306 {
        return Ok(-y.ln());
    }

    // |x| <= 10
    if y <= 10.0 {
        let tmp = match gammafn(x) {
            Ok(f) => f,
            Err(err) => { return Err(err); },
        };
        return Ok(tmp.abs().ln());
    }

    if y > XMAX {
        return Err(FishersExactError::OutOfBoundsValue);
    }

    // x > 10
    if x > 0.0 {
        if x > 1.0e17 {
            return Ok(x*x.ln() - 1.0);
        } else if x > 4934720.0 {
            return Ok(std::f64::consts::TAU.sqrt().ln() + (x-0.5)*x.ln() - x);
        } else {
            let tmp = match lgammacor(x) {
                Ok(f) => f,
                Err(err) => { return Err(err); },
            };
            return Ok(std::f64::consts::TAU.sqrt().ln() + (x-0.5)*x.ln() - x + tmp);
        }
    }

    // x < -10
    let sinpiy = (y*std::f64::consts::PI).sin().abs();

    if sinpiy < CLOSE_TO_ZERO {
        return Err(FishersExactError::OutOfBoundsValue);
    }

    let tmp = match lgammacor(y) {
        Ok(f) => f,
        Err(err) => { return Err(err); },
    };
    let ans = std::f64::consts::FRAC_PI_2.sqrt().ln() + (x-0.5)*y.ln() - x - sinpiy.ln() - tmp;

    Ok(ans)
}

// Computes natural log of the error term in Stirling's formula
//     https://cran.r-project.org/web/packages/DPQ/vignettes/log1pmx-etc.pdf
fn stirling_error(n: f64) -> Result<f64, FishersExactError> {
    const S0: f64 = 0.083333333333333333333;        // 1/12
    const S1: f64 = 0.00277777777777777777778;      // 1/360
    const S2: f64 = 0.00079365079365079365079365;   // 1/1260
    const S3: f64 = 0.000595238095238095238095238;  // 1/1680
    const S4: f64 = 0.0008417508417508417508417508; // 1/1188

    //exact values for n/2 for 0 <= n <= 30
    const SFERR_HALVES: [f64; 31] = [
        0.0, /* n=0 - wrong, place holder only */
        0.1534264097200273452913848  , /*  0.5 */
        0.0810614667953272582196702  , /*  1.0 */
        0.0548141210519176538961390  , /*  1.5 */
        0.0413406959554092940938221  , /*  2.0 */
        0.03316287351993628748511048 , /*  2.5 */
        0.02767792568499833914878929 , /*  3.0 */
        0.02374616365629749597132920 , /*  3.5 */
        0.02079067210376509311152277 , /*  4.0 */
        0.01848845053267318523077934 , /*  4.5 */
        0.01664469118982119216319487 , /*  5.0 */
        0.01513497322191737887351255 , /*  5.5 */
        0.01387612882307074799874573 , /*  6.0 */
        0.01281046524292022692424986 , /*  6.5 */
        0.01189670994589177009505572 , /*  7.0 */
        0.01110455975820691732662991 , /*  7.5 */
        0.010411265261972096497478567, /*  8.0 */
        0.009799416126158803298389475, /*  8.5 */
        0.009255462182712732917728637, /*  9.0 */
        0.008768700134139385462952823, /*  9.5 */
        0.008330563433362871256469318, /* 10.0 */
        0.007934114564314020547248100, /* 10.5 */
        0.007573675487951840794972024, /* 11.0 */
        0.007244554301320383179543912, /* 11.5 */
        0.006942840107209529865664152, /* 12.0 */
        0.006665247032707682442354394, /* 12.5 */
        0.006408994188004207068439631, /* 13.0 */
        0.006171712263039457647532867, /* 13.5 */
        0.005951370112758847735624416, /* 14.0 */
        0.005746216513010115682023589, /* 14.5 */
        0.005554733551962801371038690  /* 15.0 */
    ];

    if n <= 15.0 {
        let nn = n + n;
        if (nn - nn.round()).abs() < CLOSE_TO_ZERO {
            return Ok(SFERR_HALVES[nn as usize]);
        }

        let tmp = match lgammafn(n + 1.0) {
            Ok(f) => f,
            Err(err) => { return Err(err); },
        };
        return Ok(tmp - (n + 0.5)*n.ln() + n - std::f64::consts::TAU.sqrt().ln());
    }

    let nn = n*n;
    if n > 500.0 {
        return Ok((S0 - S1/nn) / n);
    }
    if n > 80.0 {
        return Ok((S0 - (S1 - S2/nn) / nn) / n);
    }
    if n > 35.0 {
        return Ok((S0 - (S1 - (S2 - S3/nn) / nn) / nn) / n);
    }

    // 15 < n <= 35
    Ok((S0 - (S1 - (S2 - (S3 - S4/nn) / nn) / nn) / nn) / n)
}

// Unit deviance of a Poisson distribution
//     https://en.wikipedia.org/wiki/Deviance_(statistics)
// Based on Catherine Loader's C implementation for R
fn deviance(y: f64, mean: f64) -> Result<f64, FishersExactError> {
    // Will fail if mean is 0, so break early if mean ~ 0
    if mean < CLOSE_TO_ZERO {
        return Err(FishersExactError::BadInputs);
    }

    // If y approx. equals mean (y/mean ~ 1), then evaluate using a Taylor series
    if (y - mean).abs() < 0.1*(y + mean) {
        let mut x: f64 = (y-mean) / (y+mean);
        let mut s: f64 = (y-mean) * x; // (y-mean)^2 / (y+mean)

        if s.abs() < f64::MIN {
            // Catch case where x underflows to 0.0
            return Ok(s);
        }

        let mut expansion_term: f64 = 2.0*y*x;
        x *= x;

        for j in 1..1000 {
            expansion_term *= x;
            let previous = s;
            s += expansion_term / ((j<<1)+1) as f64;
            if (s - previous).abs() < CLOSE_TO_ZERO {
                // Last expansion term was basically zero, so end here
                return Ok(s);
            }
        }

        return Err(FishersExactError::ConvergeFailure);
    }

    // If y/mean is not around 1, then use the standard equation
    Ok(y * (y/mean).ln() + mean - y)
}

// Calculate the binomial probability
//     https://www.r-project.org/doc/reports/CLoader-dbinom-2002.pdf
// Based on Catherine Loader's C implementation for R
fn binomial_prob(x: f64, n: f64, p: f64, q: f64, as_log: bool) -> Result<f64, FishersExactError> {
    // Check for quick returns
    if p.abs() < CLOSE_TO_ZERO {
        if x.abs() < CLOSE_TO_ZERO {
            if as_log { return Ok(0.0); } else { return Ok(1.0); }
        } else {
            return try_return_zero(as_log);
        }
    }

    if q.abs() < CLOSE_TO_ZERO {
        if (x-n).abs() < CLOSE_TO_ZERO {
            if as_log { return Ok(0.0); } else { return Ok(1.0); }
        } else {
            return try_return_zero(as_log);
        }
    }

    // Check for x ~ 0 or x ~ n
    if x.abs() < CLOSE_TO_ZERO {
        if n.abs() < CLOSE_TO_ZERO {
            if as_log { return Ok(0.0); } else { return Ok(1.0); }
        }

        let tmp = if p < 0.1 { 
            let t = match deviance(n, n*q) {
                Ok(f) => f,
                Err(err) => { return Err(err); },
            };

            -t - n*p
        } else { 
            n*q.ln()
        };

        if as_log {
            return Ok(tmp);
        } else {
            return Ok(tmp.exp());
        }
    }

    if (x - n).abs() < CLOSE_TO_ZERO {
        let tmp = if q < 0.1 {
            let t = match deviance(n, n*p) {
                Ok(f) => f,
                Err(err) => { return Err(err); },
            };

            -t - n*q
        } else {
            n*p.ln()
        };

        if as_log {
            return Ok(tmp);
        } else {
            return Ok(tmp.exp());
        }
    }

    if x < 0.0 || x > n {
        return try_return_zero(as_log);
    }

    let n_err = match stirling_error(n) {
        Ok(f) => f,
        Err(err) => { return Err(err); },
    };
    let x_err = match stirling_error(x) {
        Ok(f) => f,
        Err(err) => { return Err(err); },
    };
    let nx_err = match stirling_error(n-x) {
        Ok(f) => f,
        Err(err) => { return Err(err); },
    };
    let x_dev = match deviance(x, n*p) {
        Ok(f) => f,
        Err(err) => { return Err(err); },
    };
    let nx_dev = match deviance(n-x, n*q) {
        Ok(f) => f,
        Err(err) => { return Err(err); },
    };
    let correction = n_err - x_err - nx_err - x_dev - nx_dev;
    let binom_approx = std::f64::consts::TAU.ln() + x.ln() + (1.0 - x/n).ln(); // TAU = 2*pi

    if as_log {
        Ok(correction - 0.5*binom_approx)
    } else {
        Ok((correction - 0.5*binom_approx).exp())
    }
}

// Calculate hypergeometric probability of drawing d red balls (without replacement) from a bucket
// with r red balls and b blue balls in n draws
//     d <= r and d >= 0
//     n <= r+b, n >= d, and n > b+d
fn hypergeometric_prob(d: f64, r: f64, b: f64, n: f64, as_log: bool) -> Result<f64, FishersExactError> {
    // Check for inputs issues
    if d.is_nan() || r.is_nan() || b.is_nan() || n.is_nan() {
        return Err(FishersExactError::BadInputs);
    }
    if d < 0.0 || r < 0.0 || b < 0.0 || n < 0.0 {
        return Err(FishersExactError::BadInputs);
    }
    if r < d || n < d || r+b < n || b < n-d {
        return Err(FishersExactError::BadInputs);
    }

    if n.abs() < CLOSE_TO_ZERO {
        if d.abs() < CLOSE_TO_ZERO {
            if as_log { return Ok(0.0); } else { return Ok(1.0); }
        } else {
            return try_return_zero(as_log);
        }
    }

    // Probabilities
    let p = n / (r + b);
    let q = (r + b - n) / (r + b);

    let p1: f64 = match binomial_prob(d, r, p, q, as_log) {
        Ok(f) => f,
        Err(err) => { return Err(err); },
    };
    let p2: f64 = match binomial_prob(n-d, b, p, q, as_log) {
        Ok(f) => f,
        Err(err) => { return Err(err); },
    };
    let p3: f64 = match binomial_prob(n, r+b, p, q, as_log) {
        Ok(f) => f,
        Err(err) => { return Err(err); },
    };

    if as_log {
        Ok(p1 + p2 - p3)
    } else {
        Ok(p1 * p2 / p3)
    }
}

pub fn fishers_exact_test(m11: i64, m12: i64, m21: i64, m22: i64) -> f64 {
    let c1 = m11 + m21; // column 1
    let c2 = m12 + m22; // column 2
    let r1 = m11 + m12; // row 1

    let lo = 0_i64.max(r1-c2);
    let hi = r1.min(c1);
    let test_values: Vec<i64> = (lo..hi+1).collect();

    // Find hypergeometric probability for each test value
    let mut log_hypergeo: Vec<f64> = test_values
        .iter()
        .map(
            |tv| match hypergeometric_prob(*tv as f64, c1 as f64, c2 as f64, r1 as f64, true) {
                Ok(f) => f,
                Err(err) => {
                    eprintln!("Error performing fisher's exact test: {}", err);
                    quit::with_code(1);
                },
            }
        ).collect();

    // In the vein of R, do the conditional maximum likelihood estimate
    let lh_max = *log_hypergeo.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
    for x in log_hypergeo.iter_mut() {
        *x = (*x - lh_max).exp();
    }

    let lh_sum: f64 = log_hypergeo.iter().sum();
    for x in log_hypergeo.iter_mut() {
        *x /= lh_sum;
    }

    // Calculate p-value
    let comp: f64 = log_hypergeo[(m11 - lo) as usize];
    let pval = log_hypergeo.iter().filter(|&x| *x <= 1.0000001*comp).fold(0.0, |acc, e| acc + e);

    if pval > 1.0 {
        1.0
    } else {
        pval
    }
}

#[cfg(test)]
mod tests{
    use super::*;
    use crate::constants;

    fn float_equality(x: f64, y: f64) -> bool {
        (x - y).abs() < constants::CLOSE_TO_ZERO
    }

    #[test]
    fn test_try_return_zero_true() {
        let test = try_return_zero(true);
        assert_eq!(test, Err(FishersExactError::LogZero));
    }

    #[test]
    fn test_try_return_zero_false() {
        let test = try_return_zero(false).unwrap();
        assert!(float_equality(test, 0.0));
    }

    #[test]
    fn test_lgammacor_pass() {
        let test = lgammacor(20.0).unwrap();
        assert!(float_equality(test, 0.0041663196919969146));
    }

    #[test]
    fn test_lgammacor_big_pass() {
        let test = lgammacor(94906270.0).unwrap();
        assert!(float_equality(test, 0.0000000008780593034931552));
    }

    #[test]
    fn test_lgammacor_nan_error() {
        let test = lgammacor(1.0);
        assert_eq!(test, Err(FishersExactError::NaNReturn));
    }

    #[test]
    fn test_lgammacor_oob_error() {
        let test = lgammacor(f64::MAX);
        assert_eq!(test, Err(FishersExactError::OutOfBoundsValue));
    }

    #[test]
    fn test_gammafn_int() {
        let test = gammafn(1.0).unwrap();
        assert!(float_equality(test, 1.0));
    }

    #[test]
    fn test_gammafn_half_int() {
        let test = gammafn(1.5).unwrap();
        assert!(float_equality(test, 0.8862269254527580137));
    }

    #[test]
    fn test_gammafn_neg() {
        let test = gammafn(-1.5).unwrap();
        assert!(float_equality(test, 2.3632718012073548053));
    }

    #[test]
    fn test_gammafn_eleven() {
        let test = gammafn(11.0).unwrap();
        assert!(float_equality(test, 3628800.0));
    }

    #[test]
    fn test_lgammafn_zero() {
        let test = lgammafn(0.0);
        assert_eq!(test, Err(FishersExactError::InfiniteGamma));
    }

    #[test]
    fn test_lgammafn_one() {
        let test = lgammafn(1.0).unwrap();
        assert!(float_equality(test, 0.0));
    }

    #[test]
    fn test_lgammafn_minus_half() {
        let test = lgammafn(-0.5).unwrap();
        assert!(float_equality(test, 1.2655121234846453682));
    }

    #[test]
    fn test_lgammafn_eleven() {
        let test = lgammafn(11.0).unwrap();
        assert!(float_equality(test, 15.104412573075515));
    }

    #[test]
    fn test_deviance_pass() {
        let test = deviance(2.0, 1.0).unwrap();
        assert!(float_equality(test, 0.38629436112));
    }

    #[test]
    fn test_deviance_taylor_pass() {
        let test = deviance(1.1, 1.0).unwrap();
        assert!(float_equality(test, 0.00484119778));
    }

    #[test]
    fn test_deviance_taylor_same_inputs() {
        let test = deviance(1.0, 1.0).unwrap();
        assert!(float_equality(test, 0.0));
    }

    #[test]
    // Definitely don't have tests that cover all possible paths through this function...
    fn test_binomial_prob_p_x_zero() {
        let test = binomial_prob(0.0, 2.0, 0.0, 1.0, false).unwrap();
        assert!(float_equality(test, 1.0));
    }

    #[test]
    fn test_binomial_prob_log_false() {
        let test = binomial_prob(5.0, 10.0, 0.5, 0.5, false).unwrap();
        assert!(float_equality(test, 0.24609375));
    }

    #[test]
    fn test_binomial_prob_log_true() {
        let test = binomial_prob(5.0, 10.0, 0.5, 0.5, true).unwrap();
        assert!(float_equality(test, -1.40204271808803));
    }

    #[test]
    // Only testing first case in if-statement
    // Proper testing should cover all cases in statement
    fn test_hypergeometric_prob_nan() {
        let test = hypergeometric_prob(f64::NAN, 1.0, 1.0, 1.0, true);
        assert_eq!(test, Err(FishersExactError::BadInputs));
    }

    #[test]
    // Only testing first case in if-statement
    // Proper testing should cover all cases in statement
    fn test_hypergeometric_prob_negative_input() {
        let test = hypergeometric_prob(-1.0, 1.0, 1.0, 1.0, true);
        assert_eq!(test, Err(FishersExactError::BadInputs));
    }

    #[test]
    // Only testing first case in if-statement
    // Proper testing should cover all cases in statement
    fn test_hypergeometric_prob_bad_input() {
        let test = hypergeometric_prob(2.0, 1.0, 1.0, 1.0, true);
        assert_eq!(test, Err(FishersExactError::BadInputs));
    }

    #[test]
    fn test_fishers_exact_test_9125() {
        let test = fishers_exact_test(9, 1, 2, 5);
        assert!(float_equality(test, 0.0345022624434389));
    }

    #[test]
    fn test_fishers_exact_test_5103() {
        let test = fishers_exact_test(5, 1, 0, 3);
        assert!(float_equality(test, 0.0476190476190476));
    }
}
