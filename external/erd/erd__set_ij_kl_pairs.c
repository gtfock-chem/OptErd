#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <yepPredefines.h>

extern double erd__dsqmin_line_segments_ (double *, double *,
	double *, double *,
	double *, double *,
	double *, double *,
	double *, double *,
	double *, double *);

#define MIN(a,b)    ((a) > (b) ? (b) : (a))

#define TOL 1e-14

extern double mytime[10];

uint64_t erd__set_ij_kl_pairs_ticks = 0;

static YEP_INLINE double pow3o4(double x) {
	return __builtin_sqrt(x * __builtin_sqrt(x));
}

static YEP_INLINE double vector_min(const double *YEP_RESTRICT vector, size_t length) {
	double result = vector[0];
	for (size_t i = 1; i < length; i++) {
		const double element = vector[i];
		result = (element < result) ? element : result;
	}
	return result;
}

static YEP_NOINLINE void set_pairs(int npgtoa, int npgtob,
	int atomab, int equalab,
	int swaprs, double rnabsq,
	double *YEP_RESTRICT alphaa, double *YEP_RESTRICT alphab,
	double *YEP_RESTRICT ftable, int mgrid,
	double tmax, double tstep,
	double tvstep,
	uint32_t *YEP_RESTRICT nij_ptr, int *YEP_RESTRICT prima, int *YEP_RESTRICT primb,
	double *YEP_RESTRICT rho, double qmin, double smaxcd, double rminsq,
	double *YEP_RESTRICT mytable)
{
	YEP_ALIGN(128) double ssss1[128];
	YEP_ALIGN(128) double ssss2[128];

	uint32_t nij = 0;
	if (equalab) {
		for (int i = 0; i < npgtoa; i++) {
			const double a = alphaa[i];
			for (int j = 0; j <= i; j++) {
				const double b = alphab[j];
				const double p = a + b;
				const double ab = a * b;
				const double pinv = 1.0 / p;
				const double pqpinv = 1.0 / (p + qmin);
				const double t = rminsq * p * qmin * pqpinv;
				const double ssssmx = pow3o4(ab) * smaxcd * pinv * sqrt (pqpinv);

				double f0 = (t == 0.0) ? 1.0 : sqrt (M_PI / t) * 0.5;
				if (t <= tmax) {
					const int tgrid = __builtin_lround(t * tvstep);
					//   f0 = mytable[tgrid];
					const double delta = tgrid * tstep - t;
					f0 = (((((  ftable[(mgrid + 1) * tgrid + 6] * delta * 0.166666666666667
						+ ftable[(mgrid + 1) * tgrid + 5]) * delta * 0.2
						+ ftable[(mgrid + 1) * tgrid + 4]) * delta * 0.25
						+ ftable[(mgrid + 1) * tgrid + 3]) * delta * 0.333333333333333
						+ ftable[(mgrid + 1) * tgrid + 2]) * delta * 0.5
						+ ftable[(mgrid + 1) * tgrid + 1]) * delta
						+ ftable[(mgrid + 1) * tgrid + 0];
				}
				if (ssssmx * f0 >= TOL) {
					rho[nij] = 1.0;
					prima[nij] = i + 1;
					primb[nij] = j + 1;
					nij += 1;
				}
			}
		}
	} else {
		double *alphai, *alphaj;
		int *primi, *primj;
		int endi, endj;
		if (swaprs) {
			endi = npgtob;
			endj = npgtoa;
			alphai = alphab;
			alphaj = alphaa;
			primi = primb;
			primj = prima;
		} else {
			endi = npgtoa;
			endj = npgtob;
			alphai = alphaa;
			alphaj = alphab;
			primi = prima;
			primj = primb;
		}

		if (atomab) {
			for (int i = 0; i < endi; i++) {
				const double a = alphai[i];
				#pragma simd
				for (int j = 0; j < endj; j++) {
					__assume_aligned(alphaj, 64);
					const double b = alphaj[j];
					const double p = a + b;
					const double ab = a * b;
					const double pinv = 1.0 / p;
					const double pqpinv = 1.0 / (p + qmin);
					const double t = rminsq * p * qmin * pqpinv;
					const double ssssmx = pow3o4(ab) * smaxcd * pinv * sqrt (pqpinv);

					double f0 = (t == 0.0) ? 1.0 : sqrt (M_PI / t) * 0.5;
					if (t <= tmax) {
						const int tgrid = __builtin_lround(t * tvstep);
						// f0 = mytable[tgrid];
						const double delta = tgrid * tstep - t;
						f0 = (((((ftable[(mgrid + 1) * tgrid + 6] * delta * 0.166666666666667
							+ ftable[(mgrid + 1) * tgrid + 5]) * delta * 0.2
							+ ftable[(mgrid + 1) * tgrid + 4]) * delta * 0.25
							+ ftable[(mgrid + 1) * tgrid + 3]) * delta * 0.333333333333333
							+ ftable[(mgrid + 1) * tgrid + 2]) * delta * 0.5
							+ ftable[(mgrid + 1) * tgrid + 1]) * delta
							+ ftable[(mgrid + 1) * tgrid];
					}
					ssss1[j] = ssssmx * f0;
				}
				for (int j = 0; j < endj; j++) {                
					if (ssss1[j] >= TOL) {
						rho[nij] = 1.0;
						primi[nij] = i + 1;
						primj[nij] = j + 1;
						nij += 1;
					}
				}
			}
		} else {
			for (int i = 0; i < endi; i++) {
				const double a = alphai[i];
				#pragma simd
				for (int j = 0; j < endj; j++) {
					__assume_aligned(alphaj, 64);
					const double b = alphaj[j];
					const double p = a + b;
					const double ab = a * b;
					const double pinv = 1.0 / p;
					const double pqpinv = 1.0 / (p + qmin);
					const double rhoab = exp (-ab * rnabsq * pinv);
					const double t = rminsq * p * qmin * pqpinv;
					const double ssssmx = pow3o4(ab) * rhoab * smaxcd * pinv * sqrt (pqpinv);

					double f0 = (t == 0.0) ? 1.0 : sqrt (M_PI / t) * 0.5;
					if (t <= tmax) {
						const int tgrid = __builtin_lround(t * tvstep);
						f0 = mytable[tgrid];
						const double delta = tgrid * tstep - t;
						f0 = (((((ftable[(mgrid + 1) * tgrid + 6] * delta * 0.166666666666667
							  + ftable[(mgrid + 1) * tgrid + 5]) * delta * 0.2
							  + ftable[(mgrid + 1) * tgrid + 4]) * delta * 0.25
							  + ftable[(mgrid + 1) * tgrid + 3]) * delta * 0.333333333333333
							  + ftable[(mgrid + 1) * tgrid + 2]) * delta * 0.5
							  + ftable[(mgrid + 1) * tgrid + 1]) * delta
							  + ftable[(mgrid + 1) * tgrid];
					}
					ssss1[j] = ssssmx * f0;
					ssss2[j] = rhoab;
				}
				for (int j = 0; j < endj; j++) {                 
					if (ssss1[j] >= TOL) {
						rho[nij] = ssss2[j];
						primi[nij] = i + 1;
						primj[nij] = j + 1;
						nij += 1;
					}
				}
			}
		}
	}
	*nij_ptr = nij;
}

static int erd__set_ij_kl_pairs(int npgtoa, int npgtob, int npgtoc, int npgtod,
	int atomab, int atomcd, int equalab, int equalcd,
	int swaprs, int swaptu,
	double *YEP_RESTRICT xa, double *YEP_RESTRICT ya, double *YEP_RESTRICT za,
	double *YEP_RESTRICT xb, double *YEP_RESTRICT yb, double *YEP_RESTRICT zb,
	double *YEP_RESTRICT xc, double *YEP_RESTRICT yc, double *YEP_RESTRICT zc,
	double *YEP_RESTRICT xd, double *YEP_RESTRICT yd, double *YEP_RESTRICT zd,
	double rnabsq, double rncdsq, double prefact,
	double *YEP_RESTRICT alphaa, double *YEP_RESTRICT alphab,
	double *YEP_RESTRICT alphac, double *YEP_RESTRICT alphad,
	double *YEP_RESTRICT ftable, int mgrid, int ngrid,
	double tmax, double tstep,
	double tvstep, int screen, int *YEP_RESTRICT empty,
	int *YEP_RESTRICT nij_ptr, int *YEP_RESTRICT nkl_ptr, int *YEP_RESTRICT prima,
	int *YEP_RESTRICT primb, int *YEP_RESTRICT primc, int *YEP_RESTRICT primd,
	double *YEP_RESTRICT rho)
{
	/* System generated locals */
	double pmin;
	double qmin;
	double pinv;
	double qinv;
	double abmin;
	double cdmin;   
	double smaxab;
	double smaxcd;
	double rminsq;
	double mytable[64];

	const uint64_t ticks_start = __rdtsc();
	*empty = 0;
#if 0
	for (i = 0; i <= 46; i++) {
		double delta;
		int tgrid;
		double f0;
		 tgrid = (int) (i * tvstep + .5);
		 delta = tgrid * tstep - i;
		 f0 = (((((  ftable[(mgrid + 1) * tgrid + 6] * delta *
					 0.166666666666667
				   + ftable[(mgrid + 1) * tgrid + 5]) * delta * 0.2
				   + ftable[(mgrid + 1) * tgrid + 4]) * delta * 0.25
				   + ftable[(mgrid + 1) * tgrid + 3]) * delta *
				   0.333333333333333
				   + ftable[(mgrid + 1) * tgrid + 2]) * delta * 0.5
				   + ftable[(mgrid + 1) * tgrid + 1]) * delta
				   + ftable[(mgrid + 1) * tgrid + 0];
		 mytable[i] = f0;
	}
#endif

	uint32_t nij = 0;
	uint32_t nkl = 0;

	// if not screening
	if (!(screen)) {
		if (equalab) {
			for (int i = 0; i < npgtoa; i++) {
				for (int j = 0; j <= i; j++) {
					rho[nij] = 1.00;
					prima[nij] = i + 1;
					primb[nij] = j + 1;
					nij += 1;
				}
			}
		} else {
			if (swaprs) {
				if (atomab) {
					for (int j = 0; j < npgtob; j++) {
						for (int i = 0; i < npgtoa; i++) {
							rho[nij] = 1.0;
							prima[nij] = i + 1;
							primb[nij] = j + 1;
							nij += 1;
						}
					}
				} else {
					for (int j = 0; j < npgtob; j++) {
						const double b = alphab[j];
						for (int i = 0; i < npgtoa; i++) {
							const double a = alphaa[i];
							rho[nij] = exp (-a * b * rnabsq / (a + b));
							prima[nij] = i + 1;
							primb[nij] = j + 1;
							nij += 1;
						}
					}
				}
			} else {
				if (atomab) {
					for (int i = 0; i < npgtoa; i++) {
						for (int j = 0; j < npgtob; j++) {
							rho[nij] = 1.0;
							prima[nij] = i + 1;
							primb[nij] = j + 1;
							nij += 1;
						}
					}
				} else {
					for (int i = 0; i < npgtoa; i++) {
						const double a = alphaa[i];
						for (int j = 0; j < npgtob; j++) {
							const double b = alphab[j];
							rho[nij] = exp(-a * b * rnabsq / (a + b));
							prima[nij] = i + 1;
							primb[nij] = j + 1;
							nij += 1;
						}
					}
				}
			}
		}

		nkl = 0;
		if (equalcd) {
			for (int k = 0; k < npgtoc; k++) {
				for (int l = 0; l <= k; l++) {
					rho[nij + nkl] = 1.0;
					primc[nkl] = k + 1;
					primd[nkl] = l + 1;
					nkl += 1;
				}
			}
		} else {
			if (swaptu) {
				if (atomcd) {
					for (int l = 0; l < npgtod; l++) {
						for (int k = 0; k < npgtoc; k++) {
							rho[nij + nkl] = 1.0;
							primc[nkl] = k + 1;
							primd[nkl] = l + 1;
							nkl += 1;
						}
					}
				} else {
					for (int l = 0; l < npgtod; l++) {
						const double d = alphad[l];
						for (int k = 0; k < npgtoc; k++) {
							const double c = alphac[k];
							rho[nij + nkl] = exp (-c * d * rncdsq / (c + d));
							primc[nkl] = k + 1;
							primd[nkl] = l + 1;
							nkl += 1;
						}
					}
				}
			} else {
				if (atomcd) {
					for (int k = 0; k < npgtoc; k++) {
						for (int l = 0; l < npgtod; l++) {
							rho[nij + nkl] = 1.0;
							primc[nkl] = k + 1;
							primd[nkl] = l + 1;
							nkl += 1;
						}
					}
				} else {
					for (int k = 0; k < npgtoc; k++) {
						const double c = alphac[k];
						for (int l = 0; l < npgtod; l++) {
							const double d = alphad[l];
							rho[nij + nkl] = exp (-c * d * rncdsq / (c + d));
							primc[nkl] = k + 1;
							primd[nkl] = l + 1;
							nkl += 1;
						}
					}
				}
			}
		}
	}

	// compute min
	rminsq = erd__dsqmin_line_segments_(xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd);

	const double a = vector_min(alphaa, npgtoa);
	const double b = vector_min(alphab, npgtob);
	const double c = vector_min(alphac, npgtoc);
	const double d = vector_min(alphad, npgtod);

	pmin = a + b;
	qmin = c + d;
	abmin = a * b;
	cdmin = c * d;
	pinv = 1.0 / pmin;
	qinv = 1.0 / qmin;   
	smaxab = prefact * pow3o4(abmin) * exp (-abmin * rnabsq * pinv) * pinv;
	smaxcd = prefact * pow3o4(cdmin) * exp (-cdmin * rncdsq * qinv) * qinv;

	/* ...perform K2 primitive screening on A,B part. */
	set_pairs(npgtoa, npgtob, atomab, equalab, swaprs, rnabsq,
		alphaa, alphab, ftable, mgrid, tmax, tstep, tvstep,
		&nij, prima, primb, rho, qmin, smaxcd, rminsq, mytable);
	if (nij == 0) {
		*empty = 1;
		const uint64_t ticks_elapsed = __rdtsc() - ticks_start;
		erd__set_ij_kl_pairs_ticks += ticks_elapsed;

		*nij_ptr = nij;
		*nkl_ptr = nkl;
		return 0;
	}

	set_pairs(npgtoc, npgtod, atomcd, equalcd, swaptu, rncdsq,
		alphac, alphad, ftable, mgrid, tmax, tstep, tvstep,
		&nkl, primc, primd, &(rho[nij]), pmin, smaxab, rminsq, mytable);
	if (nkl == 0) {
		*empty = 1;
		const uint64_t ticks_elapsed = __rdtsc() - ticks_start;
		erd__set_ij_kl_pairs_ticks += ticks_elapsed;

		*nij_ptr = nij;
		*nkl_ptr = nkl;
		return 0;
	}

	const uint64_t ticks_elapsed = __rdtsc() - ticks_start;
	erd__set_ij_kl_pairs_ticks += ticks_elapsed;

	*nij_ptr = nij;
	*nkl_ptr = nkl;
	return 0;
}

int erd__set_ij_kl_pairs_(int *YEP_RESTRICT npgtoa, int *YEP_RESTRICT npgtob,
	int *YEP_RESTRICT npgtoc, int *YEP_RESTRICT npgtod, int *YEP_RESTRICT npgtoab,
	int *YEP_RESTRICT npgtocd, int *YEP_RESTRICT atomab, int *YEP_RESTRICT atomcd,
	int *YEP_RESTRICT equalab, int *YEP_RESTRICT equalcd, int *YEP_RESTRICT swaprs,
	int *YEP_RESTRICT swaptu, double *YEP_RESTRICT xa, double *YEP_RESTRICT ya,
	double *YEP_RESTRICT za, double *YEP_RESTRICT xb, double *YEP_RESTRICT yb,
	double *YEP_RESTRICT zb, double *YEP_RESTRICT xc, double *YEP_RESTRICT yc,
	double *YEP_RESTRICT zc, double *YEP_RESTRICT xd, double *YEP_RESTRICT yd,
	double *YEP_RESTRICT zd, double *YEP_RESTRICT rnabsq,
	double *YEP_RESTRICT rncdsq, double *YEP_RESTRICT prefact,
	double *YEP_RESTRICT alphaa, double *YEP_RESTRICT alphab,
	double *YEP_RESTRICT alphac, double *YEP_RESTRICT alphad,
	double *YEP_RESTRICT ftable, int *YEP_RESTRICT mgrid, int *YEP_RESTRICT ngrid,
	double *YEP_RESTRICT tmax, double *YEP_RESTRICT tstep,
	double *YEP_RESTRICT tvstep, int *YEP_RESTRICT screen, int *YEP_RESTRICT empty,
	int *YEP_RESTRICT nij, int *YEP_RESTRICT nkl, int *YEP_RESTRICT prima,
	int *YEP_RESTRICT primb, int *YEP_RESTRICT primc, int *YEP_RESTRICT primd,
	double *YEP_RESTRICT rho)
{
	return erd__set_ij_kl_pairs(*npgtoa, *npgtob, *npgtoc, *npgtod,
		*atomab, *atomcd,
		*equalab, *equalcd,
		*swaprs, *swaptu,
		xa, ya, za,
		xb, yb, zb,
		xc, yc, zc,
		xd, yd, zd,
		*rnabsq, *rncdsq,
		*prefact,
		alphaa, alphab, alphac, alphad,
		ftable, *mgrid, *ngrid,
		*tmax, *tstep,
		*tvstep, *screen, empty,
		nij, nkl, prima,
		primb, primc, primd,
		rho);
}

