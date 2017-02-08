#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string.h>
#include <sstream>

#include <fsft.h>
#include <fsglft.h>

void fsglft (int B, double* f_real, double* f_imag, double*** rcoeffs, double*** icoeffs, std::string path)
{
   double*** f_fsft_real = new double**[2 * B];
   double*** f_fsft_imag = new double**[2 * B];
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      f_fsft_real[i] = new double*[B];
      f_fsft_imag[i] = new double*[B];
   
      for ( int l = 0; l < B; ++l )
      {
         f_fsft_real[i][l] = new double[2 * B + 1]; 
         f_fsft_imag[i][l] = new double[2 * B + 1];
      }
   }
   
   double* rdata = new double[4 * B * B];
   double* idata = new double[4 * B * B];
   
   int* perm = new int[2 * B];
   
   for ( int index = 0; index <= B; ++index )
      perm[index] = B - index;
   
   for ( int index = 0; index < B - 1; ++index )
      perm[index + B + 1] = 2 * B - 1 - index;
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      for ( int j = 0; j < 2 * B; ++j )
      {
         for ( int k = 0; k < 2 * B; ++k )
         {
            rdata[2 * B * j + perm[k]] = f_real[2 * B * (2 * B * i + j) + k];
            idata[2 * B * j + perm[k]] = f_imag[2 * B * (2 * B * i + j) + k];
         }
      }
   
      fsft (B, rdata, idata, f_fsft_real[i], f_fsft_imag[i]);
   }
   
   delete[] perm;

   char* filename_r = new char[256];   
   char* filename_a_tilde = new char[256];
   
   sprintf (filename_r, "%s/precomp/half_range_Hermite_quadrature/%d_100_digits_abcissae.txt", path.c_str (), 2 * B);
   sprintf (filename_a_tilde, "%s/precomp/half_range_Hermite_quadrature/%d_100_digits_adapted_weights.txt", path.c_str (), 2 * B);
   
   double* r = new double[2 * B];
   double* a_tilde = new double[2 * B];
   
   std::ifstream in_r (filename_r);
   std::ifstream in_a_tilde (filename_a_tilde);
   
   std::string line_r;
   std::string line_a_tilde;
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      std::getline (in_r, line_r);
      std::getline (in_a_tilde, line_a_tilde);
      
      std::stringstream str_stream_r (line_r);
      std::stringstream str_stream_a_tilde (line_a_tilde);
         
      str_stream_r >> r[i];
      str_stream_a_tilde >> a_tilde[i];
   }
   
   delete[] filename_r;
   delete[] filename_a_tilde;
   
   delete[] rdata;
   delete[] idata;
   
   rdata = new double[2 * B];
   idata = new double[2 * B];
   
   double* data_trans_real;
   double* data_trans_imag;
   
   for ( int l = 0; l < B; ++l )
   {
      data_trans_real = new double[B - l];
      data_trans_imag = new double[B - l];
   
      for ( int m = - l; m <= l; ++m )
      {
         for ( int i = 0; i < 2 * B; ++i )
         {
            rdata[i] = a_tilde[i] * f_fsft_real[i][l][m + B - 1];
            idata[i] = a_tilde[i] * f_fsft_imag[i][l][m + B - 1];
         }
      
         DRT_Clenshaw (B, l, r, rdata, idata, data_trans_real, data_trans_imag);
      
         for ( int n = l + 1; n <= B; ++n )
         {
            rcoeffs[n - 1][l][m + l] = data_trans_real[n - l - 1];
            icoeffs[n - 1][l][m + l] = data_trans_imag[n - l - 1];
         }
      }
   
      delete[] data_trans_real;
      delete[] data_trans_imag;
   }
   
   for ( int i = 0; i < 2 * B; ++i )
   {      
      for ( int l = 0; l < B; ++l )
      {
         delete[] f_fsft_real[i][l]; 
         delete[] f_fsft_imag[i][l];
      }
   }
   
   delete[] rdata;
   delete[] idata;
   delete[] r;
   delete[] a_tilde;
   
   return;
}

double R (int n, int l, double r)
{
   int n_l_difference = n - l;
   
   double r_squared = r * r;
   
   double R;
   
   double temp[n_l_difference];
   
   temp[0] = 1.0;
   
   if ( n_l_difference > 1 )
   {
      temp[1] = 1.5 + l - r_squared;
      
      for ( int index = 2; index < n_l_difference; ++index )
      {
         temp[index] = 0.0;
      }
      
      for ( int index = 2; index < n_l_difference; ++index )
      {
         temp[index] = (((2.0 * index - 0.5 + l - r_squared) * temp[index - 1]) + ((0.5 - index - l) * temp[index - 2])) / index;
      }
   }
   
   R = temp[n_l_difference - 1];
   
   R *= pow (r, l);
   
   double factor = pow (2.0, n + 1.0) / 1.77245385090551588191;
   
   for ( int index = 2 * n - 1; index > 1; index -= 2 )
   {
      factor /= index; 
   }
   
   for ( int index = 2; index < n_l_difference; ++index )
   {
      factor *= index; 
   }
   
   R *= sqrt (factor);

   return R;
}

double R_l_plus_1_l (int l, double r)
{
   return R (l + 1, l, r) * exp (- r * r);
}

double R_l_plus_2_l (int l, double r)
{
   return R (l + 2, l, r) * exp (- r * r);
}

double alpha (int n, int l, double r_squared)
{
   return (2 * n - l - 0.5 - r_squared) / sqrt ((n + 0.5) * (n - l));
}

double beta (int n, int l)
{
   return - sqrt ((n - 0.5) * (n - l - 1) / ((n + 0.5) * (n - l)));
}

void DRT_Clenshaw (int B, int l, double* r, double* rdata, double* idata, double* data_trans_real, double* data_trans_imag)
{
   double* r_squared = new double[2 * B];
   
   for ( int i = 0; i < 2 * B; ++i )
      r_squared[i] = r[i] * r[i];
      
   double* data_temp_real = new double[4 * B];
   double* data_temp_imag = new double[4 * B];

   for ( int i = 0; i < 2 * B; ++i )
   {
      data_temp_real[2 * B + i] = R_l_plus_1_l (l, r[i]);
      data_temp_imag[2 * B + i] = data_temp_real[2 * B + i];
   }
      
   data_trans_real[0] = 0.0;
   data_trans_imag[0] = 0.0;
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      data_trans_real[0] += data_temp_real[2 * B + i] * rdata[i];
      data_trans_imag[0] += data_temp_real[2 * B + i] * idata[i];
   }
      
   if ( B - l == 1 )
      return;

   for ( int i = 0; i < 2 * B; ++i )
   {
      data_temp_real[i] = R_l_plus_2_l (l, r[i]);
      data_temp_imag[i] = data_temp_real[i];
   }
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      data_temp_real[i] = rdata[i] * data_temp_real[i];
      data_temp_imag[i] = idata[i] * data_temp_imag[i];
   }
  
   data_trans_real[1] = 0.0;
   data_trans_imag[1] = 0.0;   
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      data_trans_real[1] += data_temp_real[i];
      data_trans_imag[1] += data_temp_imag[i];
   }
      
   if ( B - l == 2 )
      return;

   for ( int i = 0; i < 2 * B; ++i )
   {
      data_temp_real[2 * B + i] = beta (l + 2, l) * rdata[i] * data_temp_real[2 * B + i];
      data_temp_imag[2 * B + i] = beta (l + 2, l) * idata[i] * data_temp_imag[2 * B + i];
   }
   
   double* temp_real = new double[2 * B];
   double* temp_imag = new double[2 * B];
   
   for ( int index = 1; index <= B - l - 2; ++index )
   {
      if ( index < B - l - 2 )   
      {
         for ( int i = 0; i < 2 * B; ++i )
         {
            temp_real[i] = data_temp_real[i];
            temp_imag[i] = data_temp_imag[i];
         }
      }

      for ( int i = 0; i < 2 * B; ++i )
      {
         data_temp_real[i] = alpha (l + 1 + index, l, r_squared[i]) * data_temp_real[i] + data_temp_real[2 * B + i];
         data_temp_imag[i] = alpha (l + 1 + index, l, r_squared[i]) * data_temp_imag[i] + data_temp_imag[2 * B + i];
      }
   
      data_trans_real[2 + index - 1] = 0.0;
      data_trans_imag[2 + index - 1] = 0.0;
      
      for ( int i = 0; i < 2 * B; ++i )
      {
         data_trans_real[2 + index - 1] += data_temp_real[i];
         data_trans_imag[2 + index - 1] += data_temp_imag[i];
      }
      
      if ( index < B - l - 2 )
      {
         for ( int i = 0; i < 2 * B; ++i )
         {
            data_temp_real[2 * B + i] = beta (l + 2 + index, l) * temp_real[i];
            data_temp_imag[2 * B + i] = beta (l + 2 + index, l) * temp_imag[i];
         }
      }
   }
   
   delete[] data_temp_real;
   delete[] data_temp_imag;
   delete[] temp_real;
   delete[] temp_imag;
   delete[] r_squared;
   
   return;
}
