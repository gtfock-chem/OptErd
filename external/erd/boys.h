#pragma once

#include <stddef.h>

#define NGRID    920
#define MGRID     10

#pragma offload_attribute(push, target(mic))
static const double tmax = 46.0;
static const double tvstep = 20.0;
static const double tstep = 0.05;

extern double boys_table[NGRID + 1][MGRID + 1];

inline double boys0(double x) {
	const double scale = 0x1.C5BF891B4EF6Bp-1; /* sqrt(pi) / 2 */
	const double limit = 0x1.0000000000000p+0; /* 1 */
	const double sqrt_x = __builtin_sqrt(x);
	const double f = scale * __builtin_erf(sqrt_x) / sqrt_x;
	return f < limit ? f : limit;
}

inline double boys1(double x) {
	const double scale = 0x1.C5BF891B4EF6Bp-2; /* sqrt(pi) / 4 */
	const double limit = 0x1.5555555555555p-2; /* 1/3 */
	const double sqrt_x = __builtin_sqrt(x);
	const double f = (scale * __builtin_erf(sqrt_x) / sqrt_x - 0.5 * __builtin_exp(-x)) / x;
	return f < limit ? f : limit;
}

inline double boys2(double x) {
	const double scale = 0x1.544FA6D47B390p-1; /* 3 * sqrt(pi) / 8 */
	const double limit = 0x1.999999999999Ap-3; /* 1/5 */
	const double sqrt_x = __builtin_sqrt(x);
	const double x_square = x * x;
	const double f = (scale * __builtin_erf(sqrt_x) / sqrt_x - (0.75 + 0.5 * x) * __builtin_exp(-x)) / x_square;
	return f < limit ? f : limit;
}

inline double boys3(double x) {
	const double scale = 0x1.A96390899A074p+0; /* 15 * sqrt(pi) / 16 */
	const double limit = 0x1.2492492492492p-3; /* 1/7 */
	const double sqrt_x = __builtin_sqrt(x);
	const double x_square = x * x;
	const double x_cube = x_square * x;
	const double f = (scale * __builtin_erf(sqrt_x) / sqrt_x - (1.875 + 1.25 * x + 0.5 * x_square) * __builtin_exp(-x)) / x_cube;
	return f < limit ? f : limit;
}

inline double boys4(double x) {
	const double scale = 0x1.74371E7866C65p+2; /* 105 * sqrt(pi) / 32 */
	const double limit = 0x1.C71C71C71C71Cp-4; /* 1/9 */
	const double sqrt_x = __builtin_sqrt(x);
	const double x_square = x * x;
	const double x_pow4 = x_square * x_square;
	const double x_cube = x_square * x;
	const double f = (scale * __builtin_erf(sqrt_x) / sqrt_x - (6.5625 + 4.375 * x + 1.75 * x_square + 0.5 * x_cube) * __builtin_exp(-x)) / x_pow4;
	return f < limit ? f : limit;
}

#pragma offload_attribute(pop)
