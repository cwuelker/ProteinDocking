#ifndef fsglft_h
#define fsglft_h

#include <string.h>
#include <fstream>

void fsglft (int B, double* f_real, double* f_imag, double*** rcoeffs, double*** icoeffs, std::string path);

double R (int n, int l, double r);

double R_l_plus_1_l (int l, double r);

double R_l_plus_2_l (int l, double r);

double alpha (int n, int l, double r_squared);

double beta (int n, int l);

void DRT_Clenshaw (int B, int l, double* r, double* rdata, double* idata, double* data_trans_real, double* data_trans_imag);

#endif
