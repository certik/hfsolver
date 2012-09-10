#include <stdio.h>
#include <algorithm>
#include <cmath>

#include <libint2.h>
extern "C" {
    #include "cints.h"
}

#define SWAP_AB    0x01
#define SWAP_CD    0x02
#define SWAP_AB_CD 0x04

#define MAXFAC 100
#define EPS 1.0E-17

static double *df = NULL; // Double factorial, filled by calc_f()

double* init_array(unsigned long int size) {
  double* result = new double[size];
  for (int i = 0; i < size; i++)
    result[i] = 0.0;
  return result;
}

void free_array(double* array) {
  delete[] array;
}

void calc_f(double *F, int n, double t) {
  int i, m, k;
  int m2;
  double t2;
  double num;
  double sum;
  double term1, term2;
  static double K = 1.0 / M_2_SQRTPI;
  double et;

  if (df == NULL) {
    df = init_array(2 * MAXFAC);
    df[0] = 1.0;
    df[1] = 1.0;
    df[2] = 1.0;
    for (i = 3; i < MAXFAC * 2; i++) {
      df[i] = (i - 1) * df[i - 2];
    }
  }

  if (t > 20.0) { /* For big t's do upward recursion */
    t2 = 2 * t;
    et = exp(-t);
    t = sqrt(t);
    F[0] = K * erf(t) / t;
    for (m = 0; m <= n - 1; m++) {
      F[m + 1] = ((2 * m + 1) * F[m] - et) / (t2);
    }
  } else { /* For smaller t's compute F with highest n using
   asymptotic series (see I. Shavitt in
   Methods in Computational Physics, ed. B. Alder eta l,
   vol 2, 1963, page 8) */
    et = exp(-t);
    t2 = 2 * t;
    m2 = 2 * n;
    num = df[m2];
    i = 0;
    sum = 1.0 / (m2 + 1);
    do {
      i++;
      num = num * t2;
      term1 = num / df[m2 + 2 * i + 2];
      sum += term1;
    } while (fabs(term1) > EPS && i < MAXFAC);
    F[n] = sum * et;
    for (m = n - 1; m >= 0; m--) { /* And then do downward recursion */
      F[m] = (t2 * F[m + 1] + et) / (2 * m + 1);
    }
  }
}

template<typename LibintEval>
void prep_libint2(LibintEval* erieval, unsigned int am1, double alpha1,
                  double A[3], unsigned int am2, double alpha2, double B[3],
                  unsigned int am3, double alpha3, double C[3],
                  unsigned int am4, double alpha4, double D[3])
{

  const unsigned int am = am1 + am2 + am3 + am4;
  if (am + 1 > 21) {
      printf("prep_libint2: am is too large\n");
      abort();
  }
  double* F = init_array(21);

  const double gammap = alpha1 + alpha2;
  const double Px = (alpha1 * A[0] + alpha2 * B[0]) / gammap;
  const double Py = (alpha1 * A[1] + alpha2 * B[1]) / gammap;
  const double Pz = (alpha1 * A[2] + alpha2 * B[2]) / gammap;
  const double PAx = Px - A[0];
  const double PAy = Py - A[1];
  const double PAz = Pz - A[2];
  const double PBx = Px - B[0];
  const double PBy = Py - B[1];
  const double PBz = Pz - B[2];
  const double AB_x = A[0] - B[0];
  const double AB_y = A[1] - B[1];
  const double AB_z = A[2] - B[2];
  const double AB2 = AB_x * AB_x + AB_y * AB_y + AB_z * AB_z;

#if LIBINT2_DEFINED(eri,PA_x)
  erieval->PA_x[0] = PAx;
#endif
#if LIBINT2_DEFINED(eri,PA_y)
  erieval->PA_y[0] = PAy;
#endif
#if LIBINT2_DEFINED(eri,PA_z)
  erieval->PA_z[0] = PAz;
#endif
#if LIBINT2_DEFINED(eri,PB_x)
  erieval->PB_x[0] = PBx;
#endif
#if LIBINT2_DEFINED(eri,PB_y)
  erieval->PB_y[0] = PBy;
#endif
#if LIBINT2_DEFINED(eri,PB_z)
  erieval->PB_z[0] = PBz;
#endif

#if LIBINT2_DEFINED(eri,AB_x)
  erieval->AB_x[0] = AB_x;
#endif
#if LIBINT2_DEFINED(eri,AB_y)
  erieval->AB_y[0] = AB_y;
#endif
#if LIBINT2_DEFINED(eri,AB_z)
  erieval->AB_z[0] = AB_z;
#endif
#if LIBINT2_DEFINED(eri,oo2z)
  erieval->oo2z[0] = 0.5/gammap;
#endif

  const double gammaq = alpha3 + alpha4;
  const double gammapq = gammap * gammaq / (gammap + gammaq);
  const double Qx = (alpha3 * C[0] + alpha4 * D[0]) / gammaq;
  const double Qy = (alpha3 * C[1] + alpha4 * D[1]) / gammaq;
  const double Qz = (alpha3 * C[2] + alpha4 * D[2]) / gammaq;
  const double QCx = Qx - C[0];
  const double QCy = Qy - C[1];
  const double QCz = Qz - C[2];
  const double QDx = Qx - D[0];
  const double QDy = Qy - D[1];
  const double QDz = Qz - D[2];
  const double CD_x = C[0] - D[0];
  const double CD_y = C[1] - D[1];
  const double CD_z = C[2] - D[2];
  const double CD2 = CD_x * CD_x + CD_y * CD_y + CD_z * CD_z;

#if LIBINT2_DEFINED(eri,QC_x)
  erieval->QC_x[0] = QCx;
#endif
#if LIBINT2_DEFINED(eri,QC_y)
  erieval->QC_y[0] = QCy;
#endif
#if LIBINT2_DEFINED(eri,QC_z)
  erieval->QC_z[0] = QCz;
#endif
#if LIBINT2_DEFINED(eri,QD_x)
  erieval->QD_x[0] = QDx;
#endif
#if LIBINT2_DEFINED(eri,QD_y)
  erieval->QD_y[0] = QDy;
#endif
#if LIBINT2_DEFINED(eri,QD_z)
  erieval->QD_z[0] = QDz;
#endif

#if LIBINT2_DEFINED(eri,CD_x)
  erieval->CD_x[0] = CD_x;
#endif
#if LIBINT2_DEFINED(eri,CD_y)
  erieval->CD_y[0] = CD_y;
#endif
#if LIBINT2_DEFINED(eri,CD_z)
  erieval->CD_z[0] = CD_z;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
  erieval->oo2e[0] = 0.5/gammaq;
#endif

  // Prefactors for interelecttron transfer relation
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_x)
  erieval->TwoPRepITR_pfac0_0_x[0] = - (alpha2*AB_x + alpha4*CD_x)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_y)
  erieval->TwoPRepITR_pfac0_0_y[0] = - (alpha2*AB_y + alpha4*CD_y)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_z)
  erieval->TwoPRepITR_pfac0_0_z[0] = - (alpha2*AB_z + alpha4*CD_z)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_x)
  erieval->TwoPRepITR_pfac0_1_x[0] = - (alpha2*AB_x + alpha4*CD_x)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_y)
  erieval->TwoPRepITR_pfac0_1_y[0] = - (alpha2*AB_y + alpha4*CD_y)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_z)
  erieval->TwoPRepITR_pfac0_1_z[0] = - (alpha2*AB_z + alpha4*CD_z)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac1_0)
  erieval->TwoPRepITR_pfac1_0[0] = -gammaq/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac1_1)
  erieval->TwoPRepITR_pfac1_1[0] = -gammap/gammaq;
#endif

  const double PQx = Px - Qx;
  const double PQy = Py - Qy;
  const double PQz = Pz - Qz;
  const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
  const double Wx = (gammap * Px + gammaq * Qx) / (gammap + gammaq);
  const double Wy = (gammap * Py + gammaq * Qy) / (gammap + gammaq);
  const double Wz = (gammap * Pz + gammaq * Qz) / (gammap + gammaq);

#if LIBINT2_DEFINED(eri,WP_x)
  erieval->WP_x[0] = Wx - Px;
#endif
#if LIBINT2_DEFINED(eri,WP_y)
  erieval->WP_y[0] = Wy - Py;
#endif
#if LIBINT2_DEFINED(eri,WP_z)
  erieval->WP_z[0] = Wz - Pz;
#endif
#if LIBINT2_DEFINED(eri,WQ_x)
  erieval->WQ_x[0] = Wx - Qx;
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
  erieval->WQ_y[0] = Wy - Qy;
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
  erieval->WQ_z[0] = Wz - Qz;
#endif
#if LIBINT2_DEFINED(eri,oo2ze)
  erieval->oo2ze[0] = 0.5/(gammap+gammaq);
#endif
#if LIBINT2_DEFINED(eri,roz)
  erieval->roz[0] = gammapq/gammap;
#endif
#if LIBINT2_DEFINED(eri,roe)
  erieval->roe[0] = gammapq/gammaq;
#endif

  double K1 = exp(-alpha1 * alpha2 * AB2 / gammap);
  double K2 = exp(-alpha3 * alpha4 * CD2 / gammaq);
  double pfac = 2 * pow(M_PI, 2.5) * K1 * K2 / (gammap * gammaq
      * sqrt(gammap + gammaq));

  calc_f(F, am, PQ2 * gammapq);

  // using dangerous macros from libint2.h
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
  erieval->LIBINT_T_SS_EREP_SS(0)[0] = pfac*F[0];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
  erieval->LIBINT_T_SS_EREP_SS(1)[0] = pfac*F[1];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
  erieval->LIBINT_T_SS_EREP_SS(2)[0] = pfac*F[2];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
  erieval->LIBINT_T_SS_EREP_SS(3)[0] = pfac*F[3];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
  erieval->LIBINT_T_SS_EREP_SS(4)[0] = pfac*F[4];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
  erieval->LIBINT_T_SS_EREP_SS(5)[0] = pfac*F[5];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
  erieval->LIBINT_T_SS_EREP_SS(6)[0] = pfac*F[6];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
  erieval->LIBINT_T_SS_EREP_SS(7)[0] = pfac*F[7];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
  erieval->LIBINT_T_SS_EREP_SS(8)[0] = pfac*F[8];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
  erieval->LIBINT_T_SS_EREP_SS(9)[0] = pfac*F[9];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
  erieval->LIBINT_T_SS_EREP_SS(10)[0] = pfac*F[10];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
  erieval->LIBINT_T_SS_EREP_SS(11)[0] = pfac*F[11];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
  erieval->LIBINT_T_SS_EREP_SS(12)[0] = pfac*F[12];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
  erieval->LIBINT_T_SS_EREP_SS(13)[0] = pfac*F[13];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
  erieval->LIBINT_T_SS_EREP_SS(14)[0] = pfac*F[14];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
  erieval->LIBINT_T_SS_EREP_SS(15)[0] = pfac*F[15];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
  erieval->LIBINT_T_SS_EREP_SS(16)[0] = pfac*F[16];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
  erieval->LIBINT_T_SS_EREP_SS(17)[0] = pfac*F[17];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
  erieval->LIBINT_T_SS_EREP_SS(18)[0] = pfac*F[18];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
  erieval->LIBINT_T_SS_EREP_SS(19)[0] = pfac*F[19];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
  erieval->LIBINT_T_SS_EREP_SS(20)[0] = pfac*F[20];
#endif

  free_array(F);
}

void swap(int *a, int *b)
{
    int tmp;
    tmp = *a;
    *a = *b;
    *b = tmp;
}

bool can_compute_permut(int i, int j, int k, int l, int permut)
{
    // We need to apply the last swap first
    if (permut & SWAP_AB_CD) {
        swap(&i, &k);
        swap(&j, &l);
    }
    if (permut & SWAP_CD) swap(&k, &l);
    if (permut & SWAP_AB) swap(&i, &j);
    return (i >= j && k >= l && i*(i-1)/2+j >= k*(k-1)/2+l);
}

int max(int a, int b, int c, int d)
{
    return std::max(std::max(a, b), std::max(c, d));
}

extern "C" {

void getInts2(double A[3], double B[3], double C[3], double D[3],
        int nprima, int nprimb, int nprimc, int nprimd,
        double coefa[], double coefb[], double coefc[],
            double coefd[],
        double alphaa[], double alphab[], double alphac[],
            double alphad[],
        int ishell, int jshell, int kshell, int lshell,
        int ilambda, int jlambda, int klambda, int llambda,
        int permut,
        int n2,
        double r[])
{
    if (ilambda < jlambda) {
        getInts2(B, A, C, D,
                nprimb, nprima, nprimc, nprimd,
                coefb, coefa, coefc, coefd,
                alphab, alphaa, alphac, alphad,
                jshell, ishell, kshell, lshell,
                jlambda, ilambda, klambda, llambda,
                permut | SWAP_AB,
                n2, r);
        return;
    }
    if (klambda < llambda) {
        getInts2(A, B, D, C,
                nprima, nprimb, nprimd, nprimc,
                coefa, coefb, coefd, coefc,
                alphaa, alphab, alphad, alphac,
                ishell, jshell, lshell, kshell,
                ilambda, jlambda, llambda, klambda,
                permut | SWAP_CD,
                n2, r);
        return;
    }
    if (ilambda + jlambda > klambda + llambda) {
        getInts2(C, D, A, B,
                nprimc, nprimd, nprima, nprimb,
                coefc, coefd, coefa, coefb,
                alphac, alphad, alphaa, alphab,
                kshell, lshell, ishell, jshell,
                klambda, llambda, ilambda, jlambda,
                permut | SWAP_AB_CD,
                n2, r);
        return;
    }

    bool can_compute;
    can_compute = (ilambda >= jlambda) && (klambda >= llambda)
                          && (ilambda + jlambda <= klambda + llambda);
    if (can_compute == false) {
        printf("internal error: permutation is wrong\n");
        abort();
    }

    LIBINT2_PREFIXED_NAME(libint2_static_init)();
    int MAXprim = nprima*nprimb*nprimc*nprimd;
    Libint_t inteval[MAXprim];
    const unsigned int ammax = max(ilambda, jlambda, klambda, llambda);
    int n;
    for (n=0; n<MAXprim; n++) {
        LIBINT2_PREFIXED_NAME(libint2_init_eri)(&inteval[n], ammax, 0);
    }
    if (LIBINT2_PREFIXED_NAME (libint2_build_eri)
        [ilambda][jlambda][klambda][llambda] == NULL) {
        if (ilambda == 0 && jlambda == 0 && klambda == 0 && llambda == 0) {
            printf("The (ss|ss) case is not handled.\n");
            abort();
        } else {
            printf("The ordering is wrong.\n");
            abort();
        }
    }

    int ip, jp, kp, lp;
    n = 0;
    for (ip=0; ip < nprima; ip++)
    for (jp=0; jp < nprimb; jp++)
    for (kp=0; kp < nprimc; kp++)
    for (lp=0; lp < nprimd; lp++) {
        prep_libint2(&inteval[n],
                        ilambda, alphaa[ip], A,
                        jlambda, alphab[jp], B,
                        klambda, alphac[kp], C,
                        llambda, alphad[lp], D);
        inteval[n].contrdepth = 1;
        LIBINT2_PREFIXED_NAME (libint2_build_eri)
            [ilambda][jlambda][klambda][llambda](&inteval[n]);
        n++;
    }

    int i, j, k, l;
    int ijkl;
    double integ;
    ijkl=0;
    for (i=ishell; i < ishell + (ilambda+1)*(ilambda+2)/2; i++)
    for (j=jshell; j < jshell + (jlambda+1)*(jlambda+2)/2; j++)
    for (k=kshell; k < kshell + (klambda+1)*(klambda+2)/2; k++)
    for (l=lshell; l < lshell + (llambda+1)*(llambda+2)/2; l++) {
        if (!can_compute_permut(i, j, k, l, permut)) {
            ijkl++;
            continue;
        }
        n = 0;
        integ=0;
        for (ip=0; ip < nprima; ip++)
        for (jp=0; jp < nprimb; jp++)
        for (kp=0; kp < nprimc; kp++)
        for (lp=0; lp < nprimd; lp++) {
            integ += inteval[n].targets[0][ijkl]
                * coefa[ip+nprima*(i-ishell)]
                * coefb[jp+nprimb*(j-jshell)]
                * coefc[kp+nprimc*(k-kshell)]
                * coefd[lp+nprimd*(l-lshell)];
            n++;
        }
        r[ijkl2intindex(i-1,j-1,k-1,l-1)] = integ;
        ijkl++;
    }
    for (n=0; n<MAXprim; n++) {
        LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[n]);
    }
}

} // extern "C"
