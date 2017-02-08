const double pi = 3.141592653589793238462643383279502884197169399;

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <pdb.h> // interface for pdb files

int main ()
{
   std::string path = "/Users/wuelker/Documents/ProteinDocking";

   const int B = 32;
   
   const double lambda = 20.0;

   clock_t begin = clock ();
   
   printf ("\n>> Starting FRM for Benzamidine and Trypsin with B = %d, lambda = %.2f... \n\n", B, lambda);

   ligand Benzamidine ("/Users/wuelker/Documents/ProteinDocking/pdb/ligands/Benzamidine.pdb");
   
   Benzamidine.create_artificial_skin (100);
   Benzamidine.evaluate_affinity_function (B, lambda, path);
   Benzamidine.compute_SGL_Fourier_coefficients (B, path);
      
   protein Trypsin ("/Users/wuelker/Documents/ProteinDocking/pdb/proteins/Trypsin.pdb");
   
   Trypsin.create_skin (10);
   Trypsin.evaluate_affinity_function (B, lambda, path);
   Trypsin.compute_SGL_Fourier_coefficients (B, path);
   
   // HERE COMES THE DOCKING
   //
   // Use the SGL Fourier coefficients like this (here real part of SGL Fourier coefficient of Trypsin): 
   //
   //    Trypsin.get_SGL_Fourier_coefficient_real (n, l, m);
   //
   // If you want to translate the ligand along the z axis and get the SGL Fourier coefficients, you can do it like this:
   //
   //    Benzamidine.translate_along_z_axis (distance);
   //    Benzamitine.evaluate_affinity_function (B, lambda, path);
   //    Benzamidine.compute_SGL_Fourier_coefficients (B, path);
   
   clock_t end = clock ();

   double elapsed_msecs = 1000.0 * (double)(end - begin) / CLOCKS_PER_SEC;
   
   printf (">> computation time: %f ms\n\n", elapsed_msecs);

   return 0;
}