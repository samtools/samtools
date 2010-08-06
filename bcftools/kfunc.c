#include <math.h>


/* Log gamma function
 * \log{\Gamma(z)}
 * AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
 */
double kf_lgamma(double z)
{
	double x = 0;
	x += 0.1659470187408462e-06 / (z+7);
	x += 0.9934937113930748e-05 / (z+6);
	x -= 0.1385710331296526     / (z+5);
	x += 12.50734324009056      / (z+4);
	x -= 176.6150291498386      / (z+3);
	x += 771.3234287757674      / (z+2);
	x -= 1259.139216722289      / (z+1);
	x += 676.5203681218835      / z;
	x += 0.9999999999995183;
	return log(x) - 5.58106146679532777 - z + (z-0.5) * log(z+6.5);
}

/* complementary error function
 * \frac{2}{\sqrt{\pi}} \int_x^{\infty} e^{-t^2} dt
 * AS66, 2nd algorithm, http://lib.stat.cmu.edu/apstat/66
 */
double kf_erfc(double x)
{
	const double p0 = 220.2068679123761;
	const double p1 = 221.2135961699311;
	const double p2 = 112.0792914978709;
	const double p3 = 33.912866078383;
	const double p4 = 6.37396220353165;
	const double p5 = .7003830644436881;
	const double p6 = .03526249659989109;
	const double q0 = 440.4137358247522;
	const double q1 = 793.8265125199484;
	const double q2 = 637.3336333788311;
	const double q3 = 296.5642487796737;
	const double q4 = 86.78073220294608;
	const double q5 = 16.06417757920695;
	const double q6 = 1.755667163182642;
	const double q7 = .08838834764831844;
	double expntl, z, p;
	z = fabs(x) * M_SQRT2;
	if (z > 37.) return x > 0.? 0. : 2.;
	expntl = exp(z * z * - .5);
	if (z < 10. / M_SQRT2) // for small z
	    p = expntl * ((((((p6 * z + p5) * z + p4) * z + p3) * z + p2) * z + p1) * z + p0)
			/ (((((((q7 * z + q6) * z + q5) * z + q4) * z + q3) * z + q2) * z + q1) * z + q0);
	else p = expntl / 2.506628274631001 / (z + 1. / (z + 2. / (z + 3. / (z + 4. / (z + .65)))));
	return x > 0.? 2. * p : 2. * (1. - p);
}

/* Regularized (incomplete lower) gamma function
 * \frac{\gamma(p,x)}{\Gamma(p)}=\frac{1}{\Gamma(p)} \int_0^x t^{p-1}e^{-t} dt
 * AS245, http://lib.stat.cmu.edu/apstat/245
 */
double kf_gammap(double p, double x)
{
    double ret_val;
    double a, b, c, an, rn, pn1, pn2, pn3, pn4, pn5, pn6, arg;

	if (x == 0.) return 0.;
	// The following line is not thoroughly tested, so it is commented out.
	if (p > 1e3) return .5 * kf_erfc(-M_SQRT1_2 * sqrt(p) * 3. * (pow(x / p, 1./3.) + 1. / (p * 9.) - 1.));
	if (x > 1e8) return 1.;
	if (x <= 1. || x < p) { // series expansion
		c = 1.;
		arg = p * log(x) - x - kf_lgamma(p + 1.);
		ret_val = 1.;
		a = p;
		while (c > 1e-14) {
			a += 1.;
			c = c * x / a;
			ret_val += c;
		}
		arg += log(ret_val);
		ret_val = 0.;
		if (arg >= -88.) ret_val = exp(arg);
	} else { // continued expansion
		arg = p * log(x) - x - kf_lgamma(p);
		a = 1. - p;
		b = a + x + 1.;
		c = 0.;
		pn1 = 1.;
		pn2 = x;
		pn3 = x + 1.;
		pn4 = x * b;
		ret_val = pn3 / pn4;
		while (1) {
			a += 1.; b += 2.; c += 1.;
			an = a * c;
			pn5 = b * pn3 - an * pn1;
			pn6 = b * pn4 - an * pn2;
			if (fabs(pn6) > 0.) {
				rn = pn5 / pn6;
				if (fabs(ret_val - rn) <= fmin(1e-14, rn * 1e-14)) break;
				ret_val = rn;
			}
			pn1 = pn3; pn2 = pn4; pn3 = pn5; pn4 = pn6;
			if (fabs(pn5) >= 1e37)
				pn1 /= 1e37, pn2 /= 1e37, pn3 /= 1e37, pn4 /= 1e37;
		}
		arg += log(ret_val);
		ret_val = 1.;
		if (arg >= -88.) ret_val = 1. - exp(arg);
    }
	return ret_val;
}

/* Numerical Recipe separates series expansion and continued
 * expansion. This may potentially reduce underflow for some
 * combinations of p and x. Nonetheless, the precision here is good
 * enough for me. I will not spend more time for now.
 */
double kf_gammaq(double p, double x)
{
	return 1. - kf_gammap(p, x);
}

#ifdef KF_MAIN
#include <stdio.h>
int main(int argc, char *argv[])
{
	double x = 10, y = 2.5;
	printf("erfc(%lg): %lg, %lg\n", x, erfc(x), kf_erfc(x));
	printf("lower-gamma(%lg,%lg): %lg\n", x, y, (1.0-kf_gammap(y, x))*tgamma(y));
	return 0;
}
#endif
