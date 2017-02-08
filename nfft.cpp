#include <nfft.h>
#include <nfft3mp.h>

static void simple_test_nfft_2d(void)
{
   int N[5], n[5], M;

   NFFT(plan) p;

   N[0] = 8; n[0] = 8;
   N[1] = 8; n[1] = 8;
   N[2] = 8; n[2] = 8;
   N[3] = 8; n[3] = 8;
   N[4] = 8; n[4] = 8;    

   M = 1000;

   /** init a two dimensional plan */
   NFFT(init_guru)(&p, 5, N, M, n, 6, PRE_PHI_HUT| PRE_FULL_PSI| MALLOC_F_HAT| MALLOC_X| MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE, FFTW_ESTIMATE| FFTW_DESTROY_INPUT);

   /** init pseudo random nodes */
   NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total);

   /** precompute psi, the entries of the matrix B */
   if(p.flags & PRE_ONE_PSI)
      NFFT(precompute_one_psi)(&p);

   /** init pseudo random Fourier coefficients and show them */
   NFFT(vrand_unit_complex)(p.f_hat, p.N_total);

   NFFT(trafo)(&p);
   NFFT(adjoint)(&p);

   /** finalise the two dimensional plan */
   NFFT(finalize)(&p);

   return;
}