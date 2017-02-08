#ifndef ligand_h
#define ligand_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>
#include <ctype.h>

/*
 * ligand class
 */
class ligand
{
public:

   ligand (std::string filename);
   
   virtual ~ligand () {}
   
   void create_artificial_skin (int N);
   void evaluate_affinity_function (int B, double lambda, std::string path);
   void compute_SGL_Fourier_coefficients (int B, std::string path);
   
   double get_SGL_Fourier_coefficient_real (int n, int l, int m) { return SGL_Fourier_coefficients_real[n - 1][l][m + l]; };
   double get_SGL_Fourier_coefficient_imag (int n, int l, int m) { return SGL_Fourier_coefficients_imag[n - 1][l][m + l]; };
   
private:
   
   int n_atoms, n_H_atoms, n_N_atoms, n_C_atoms, n_O_atoms, n_S_atoms, n_unknown_atoms, n_skin_atoms;   
   
   double* coordinates;
   char* type;
   
   double* H_coordinates, * N_coordinates, * C_coordinates, * O_coordinates, * S_coordinates, * unknown_atom_coordinates, * HS_coordinates;
   
   double* affinity_function_real, * affinity_function_imag;
   
   double*** SGL_Fourier_coefficients_real, *** SGL_Fourier_coefficients_imag;
};

/*
 * constructor
 */
ligand::ligand (std::string filename)
{
   printf (">> loading ligand from pdb file \"%s\"... ", filename.c_str());
   fflush (stdout);

   std::ifstream in (filename.c_str ());

   std::string line;

   for ( int index = 0; index < 5; ++index )
      std::getline (in, line);
   
   n_atoms = 0;
   
   std::string buffer;
   
   do
   {
      std::getline (in, line);
   
      std::stringstream str_stream (line);
      
      str_stream >> buffer;
   
      if ( buffer == "CONECT" )
         break;
   
      ++n_atoms;
   
   } while ( true );

   in.close ();
   in.open (filename.c_str ());
   
   for ( int index = 0; index < 5; ++index )
      std::getline (in, line);   
      
   n_H_atoms = 0;
   n_N_atoms = 0;
   n_C_atoms = 0;
   n_O_atoms = 0;
   n_S_atoms = 0;
   n_unknown_atoms = 0;

   coordinates = new double[3 * n_atoms];
   type = new char[n_atoms];

   for ( int index = 0; index < n_atoms; ++index )
   {
      std::getline (in, line);

      double position_x, position_y, position_z;

      char atom_type;

      std::stringstream str_stream (line);

      for ( int skip = 0; skip < 6; ++skip )
         str_stream >> buffer;
   
      if ( !isdigit (buffer.c_str()[0]) )
         str_stream >> buffer;
   
      str_stream >> position_x;
      str_stream >> position_y;
      str_stream >> position_z;
   
      coordinates[3 * index + 0] = position_x;
      coordinates[3 * index + 1] = position_y;
      coordinates[3 * index + 2] = position_z;
   
      for ( int skip = 0; skip < 2; ++skip )
         str_stream >> buffer;
   
      str_stream >> atom_type;
	  
	   switch ( atom_type )
	   {
         case 'H' : ++n_H_atoms;
		      break;
		   case 'N' : ++n_N_atoms;
		      break;
		   case 'C' : ++n_C_atoms;
		      break;
		   case 'O' : ++n_O_atoms;
		      break;
         case 'S' : ++n_S_atoms;
		      break;
         default : ++n_unknown_atoms;
		      break;
      }

      type[index] = atom_type;
   }
   
   H_coordinates = new double[3 * n_H_atoms];
   N_coordinates = new double[3 * n_N_atoms];
   C_coordinates = new double[3 * n_C_atoms];
   O_coordinates = new double[3 * n_O_atoms];
   S_coordinates = new double[3 * n_S_atoms];
   unknown_atom_coordinates = new double[3 * n_unknown_atoms];
   
   int H_counter = 0;
   int N_counter = 0;
   int C_counter = 0;
   int O_counter = 0;
   int S_counter = 0;
   int unknown_atom_counter = 0;
   
   for ( int index = 0; index < n_atoms; ++index )
   {
      switch ( type[index] )
      {
	      case 'H' : 
		      H_coordinates[3 * H_counter + 0] = coordinates[3 * index + 0];
			   H_coordinates[3 * H_counter + 1] = coordinates[3 * index + 1];
			   H_coordinates[3 * H_counter + 2] = coordinates[3 * index + 2];
			   ++H_counter;
		      break;
	      case 'N' : 
		      N_coordinates[3 * N_counter + 0] = coordinates[3 * index + 0];
			   N_coordinates[3 * N_counter + 1] = coordinates[3 * index + 1];
			   N_coordinates[3 * N_counter + 2] = coordinates[3 * index + 2];
			   ++N_counter;
		      break;
         case 'C' : 
		      C_coordinates[3 * C_counter + 0] = coordinates[3 * index + 0];
			   C_coordinates[3 * C_counter + 1] = coordinates[3 * index + 1];
			   C_coordinates[3 * C_counter + 2] = coordinates[3 * index + 2];
			   ++C_counter;
		      break;
		   case 'O' : 
		      O_coordinates[3 * O_counter + 0] = coordinates[3 * index + 0];
			   O_coordinates[3 * O_counter + 1] = coordinates[3 * index + 1];
			   O_coordinates[3 * O_counter + 2] = coordinates[3 * index + 2];
			   ++O_counter;
		      break;
		   case 'S' : 
		      S_coordinates[3 * S_counter + 0] = coordinates[3 * index + 0];
			   S_coordinates[3 * S_counter + 1] = coordinates[3 * index + 1];
			   S_coordinates[3 * S_counter + 2] = coordinates[3 * index + 2];
			   ++S_counter;
		      break;
		   default : 
		      unknown_atom_coordinates[3 * unknown_atom_counter + 0] = coordinates[3 * index + 0];
			   unknown_atom_coordinates[3 * unknown_atom_counter + 1] = coordinates[3 * index + 1];
			   unknown_atom_coordinates[3 * unknown_atom_counter + 2] = coordinates[3 * index + 2];
			   ++unknown_atom_counter;		 
		      break;
      }
   }

   in.close ();

   printf ("done.\n");
   printf ("   -> number of atoms: %d (%d H, %d N, %d C, %d O, %d S, %d unknown)\n\n", n_atoms, n_H_atoms, n_N_atoms, n_C_atoms, n_O_atoms, n_S_atoms, n_unknown_atoms);
}

void ligand::create_artificial_skin (int N)
{
   printf (">> creating artificial skin for ligand... ");
   fflush (stdout);
   
   std::vector<double> skin_x;
   std::vector<double> skin_y;
   std::vector<double> skin_z;
   
   for ( int atom_index = 0; atom_index < n_atoms; ++atom_index )
   {
      double r;
   
      switch ( type[atom_index] )
      {
         case 'O' : r = 1.52;
            break;
         case 'H' : r = 1.1;
            break;
         case 'C' : r = 1.7;
            break;
         case 'N' : r = 1.55;
            break;
         case 'S' : r = 1.8;
            break;
         default : r = 1.0;
            break;
      }
   
      double a = 4 * pi * r * r / N;
      double d = sqrt (a);
   
      int M_theta = static_cast<int> (pi / d + 0.5);
   
      double d_theta = pi / M_theta;
      double d_phi = a / d_theta;
   
      for ( int m = 0; m < M_theta; ++m )
      {
         double theta = pi * (m + 0.5) / M_theta;
      
         int M_phi = static_cast<int> (2 * pi * sin (theta) / d_phi + 0.5);
         
         for ( int n = 0; n < M_phi; ++n )
         {
            double phi = 2 * pi * n / M_phi;
         
            double x = r * sin (theta) * cos (phi) + coordinates[3 * atom_index + 0];
            double y = r * sin (theta) * sin (phi) + coordinates[3 * atom_index + 1];
            double z = r * cos (theta) + coordinates[3 * atom_index + 2];
         
            bool add = true;
         
            for ( int atom_index_2 = 0; atom_index_2 < n_atoms; ++atom_index_2 )
            {
               double r_neighboor;
            
               switch ( type[atom_index_2] )
               {
                  case 'O' : r_neighboor = 1.52;
                     break;
                  case 'H' : r_neighboor = 1.1;
                     break;
                  case 'C' : r_neighboor = 1.7;
                     break;
                  case 'N' : r_neighboor = 1.55;
                     break;
                  case 'S' : r_neighboor = 1.8;
                     break;
                  default : r_neighboor = 1.0;
                     break;
               }
            
               double distance_squared = 0.0; 
            
               distance_squared += (x - coordinates[3 * atom_index_2 + 0]) * (x - coordinates[3 * atom_index_2 + 0]);
               distance_squared += (y - coordinates[3 * atom_index_2 + 1]) * (y - coordinates[3 * atom_index_2 + 1]);
               distance_squared += (z - coordinates[3 * atom_index_2 + 2]) * (z - coordinates[3 * atom_index_2 + 2]);
            
               if ( distance_squared < r_neighboor * r_neighboor - 1e-10 )
               {
                  add = false;
               
                  break;
               }
            }
         
            if ( add == true )
            {
               skin_x.push_back (x);
               skin_y.push_back (y);
               skin_z.push_back (z);
            }
         }
      
      }
   }
   
   n_skin_atoms = (int)skin_x.size ();
   
   HS_coordinates = new double[3 * n_skin_atoms];
   
   for ( int skin_atom_index = 0; skin_atom_index < n_skin_atoms; ++skin_atom_index )
   {
      HS_coordinates[3 * skin_atom_index + 0] = skin_x.at (skin_atom_index);
      HS_coordinates[3 * skin_atom_index + 1] = skin_y.at (skin_atom_index);
      HS_coordinates[3 * skin_atom_index + 2] = skin_z.at (skin_atom_index);
   }
   
   skin_x.clear ();
   skin_y.clear ();
   skin_z.clear ();
   
   printf ("done.\n");
   printf ("   -> number of artificial skin atoms: %d\n\n", n_skin_atoms);
   
   return;
}

void ligand::evaluate_affinity_function (int B, double lambda, std::string path)
{
   printf (">> evaluating affinity function for ligand using FGT... ");
   fflush (stdout);
   
   double sqrt_lambda = sqrt (lambda);
   double lambda_pow = pow (lambda, 0.75);
   
   double special_sqrt = sqrt (2.31);
   double special_exp = exp (2.31);
   
   affinity_function_real = new double[8 * B * B * B];
   affinity_function_imag = new double[8 * B * B * B];
   
   for ( int index = 0; index < 8 * B * B * B; ++index )
   {
      affinity_function_real[index] = 0.0;   
      affinity_function_imag[index] = 0.0;
   }
   
   double* affinity_function_buffer = new double[8 * B * B * B];
   
   double* grid = new double[3 * 8 * B * B * B];
   
   char* filename = new char[256];
   
   sprintf (filename, "%s/precomp/half_range_Hermite_quadrature/%d_100_digits_abcissae.txt", path.c_str (), 2 * B);
   
   double* r = new double[2 * B];
   
   std::ifstream in (filename);
   
   std::string line;
   
   for ( int i = 0; i < 2 * B; ++i )
   {
      std::getline (in, line);
   
      std::stringstream str_stream (line);
   
      str_stream >> r[i];
   }
   
   double* exp_array = new double[2 * B];
   
   for ( int i = 0; i < 2 * B; ++i )
      exp_array[i] = exp (0.5 * r[i] * r[i]);
   
   int grid_point_index = 0;
   
   for ( int i = 0; i < 2 * B; ++i) 
   {
      for ( int j = 0; j < 2 * B; ++j)
      {
         for ( int k = 0; k < 2 * B; ++k )
         {
            grid[3 * grid_point_index + 0] = sqrt_lambda * r[i] * sin ((2 * j + 1) * pi / (4 * B)) * cos (k * pi / B);
            grid[3 * grid_point_index + 1] = sqrt_lambda * r[i] * sin ((2 * j + 1) * pi / (4 * B)) * sin (k * pi / B);
            grid[3 * grid_point_index + 2] = sqrt_lambda * r[i] * cos ((2 * j + 1) * pi / (4 * B));
         
            ++grid_point_index;
         }
      }
   }
   
   // ---------------
   
   if ( n_H_atoms > 0 )
   {
      double* H_gamma = new double[n_H_atoms];
   
      for ( int index = 0; index < n_H_atoms; ++index )
         H_gamma[index] = 16.0 * lambda_pow * special_exp;
   
      fgt (n_H_atoms, 8 * B * B * B, H_coordinates, grid, affinity_function_buffer, 1.1 / special_sqrt, H_gamma, 1e-15);
   
      grid_point_index = 0;
   
      for ( int i = 0; i < 2 * B; ++i ) 
      {
         for ( int j = 0; j < 2 * B; ++j )
         {
            for ( int k = 0; k < 2 * B; ++k )
            {
               affinity_function_imag[grid_point_index] += affinity_function_buffer[grid_point_index] * exp_array[i];
         
               ++grid_point_index;
            }
         }
      }
   
      delete[] H_gamma;
   }
   
   // ---------------
   
   if ( n_N_atoms > 0 )
   {
      double* N_gamma = new double[n_N_atoms];
   
      for ( int index = 0; index < n_N_atoms; ++index )
         N_gamma[index] = 16.0 * lambda_pow * special_exp;
   
      fgt (n_N_atoms, 8 * B * B * B, N_coordinates, grid, affinity_function_buffer, 1.55 / special_sqrt, N_gamma, 1e-15);
   
      grid_point_index = 0;
   
      for ( int i = 0; i < 2 * B; ++i ) 
      {
         for ( int j = 0; j < 2 * B; ++j )
         {
            for ( int k = 0; k < 2 * B; ++k )
            {
                affinity_function_imag[grid_point_index] += affinity_function_buffer[grid_point_index] * exp_array[i];
            
                ++grid_point_index;
            }
         }
      }
   
      delete[] N_gamma;
   }
   
   // ---------------
   
   if ( n_O_atoms > 0 )
   {
      double* O_gamma = new double[n_O_atoms];
   
      for ( int index = 0; index < n_O_atoms; ++index )
         O_gamma[index] = 16.0 * lambda_pow * special_exp;
   
      fgt (n_O_atoms, 8 * B * B * B, O_coordinates, grid, affinity_function_buffer, 1.52 / special_sqrt, O_gamma, 1e-15);
   
      grid_point_index = 0;
   
      for ( int i = 0; i < 2 * B; ++i ) 
      {
         for ( int j = 0; j < 2 * B; ++j )
         {
            for ( int k = 0; k < 2 * B; ++k )
            {
               affinity_function_imag[grid_point_index] += affinity_function_buffer[grid_point_index] * exp_array[i];
            
               ++grid_point_index;
            }
         }
      }
   
      delete[] O_gamma;   
   }
   
   // ---------------
   
   if ( n_C_atoms > 0 )
   {   
      double* C_gamma = new double[n_C_atoms];
   
      for ( int index = 0; index < n_C_atoms; ++index )
         C_gamma[index] = 16.0 * lambda_pow * special_exp;
   
      fgt (n_C_atoms, 8 * B * B * B, C_coordinates, grid, affinity_function_buffer, 1.7 / special_sqrt, C_gamma, 1e-15);
   
      grid_point_index = 0;
   
      for ( int i = 0; i < 2 * B; ++i ) 
      {
         for ( int j = 0; j < 2 * B; ++j )
         {
            for ( int k = 0; k < 2 * B; ++k )
            {
               affinity_function_imag[grid_point_index] += affinity_function_buffer[grid_point_index] * exp_array[i];
            
               ++grid_point_index;
            }
         }
      }
   
      delete[] C_gamma;     
   }
   
   // ---------------
   
   if ( n_S_atoms > 0 )
   {
      double* S_gamma = new double[n_S_atoms];
   
      for ( int index = 0; index < n_S_atoms; ++index )
         S_gamma[index] = 16.0 * lambda_pow * special_exp;
   
      fgt (n_S_atoms, 8 * B * B * B, S_coordinates, grid, affinity_function_buffer, 1.8 / special_sqrt, S_gamma, 1e-15);
   
      grid_point_index = 0;
   
      for ( int i = 0; i < 2 * B; ++i ) 
      {
         for ( int j = 0; j < 2 * B; ++j )
         {
            for ( int k = 0; k < 2 * B; ++k )
            {
               affinity_function_imag[grid_point_index] += affinity_function_buffer[grid_point_index] * exp_array[i];
            
               ++grid_point_index;
            }
         }
      }
   
      delete[] S_gamma;      
   }
   
   // ---------------

   if ( n_unknown_atoms > 0 )
   {
      double* U_gamma = new double[n_unknown_atoms];
   
      for ( int index = 0; index < n_unknown_atoms; ++index )
         U_gamma[index] = 16.0 * lambda_pow * special_exp;
   
      fgt (n_unknown_atoms, 8 * B * B * B, unknown_atom_coordinates, grid, affinity_function_buffer, 1 / special_sqrt, U_gamma, 1e-15);
   
      grid_point_index = 0;
   
      for ( int i = 0; i < 2 * B; ++i ) 
      {
         for ( int j = 0; j < 2 * B; ++j )
         {
            for ( int k = 0; k < 2 * B; ++k )
            {
               affinity_function_imag[grid_point_index] += affinity_function_buffer[grid_point_index] * exp_array[i];
            
               ++grid_point_index;
            }
         }
      }
   
      delete[] U_gamma;
   }
   
   // ---------------
   
   if ( n_skin_atoms > 0 )
   {
      double* HS_gamma = new double[n_skin_atoms];
   
      for ( int index = 0; index < n_skin_atoms; ++index )
         HS_gamma[index] = lambda_pow * special_exp;
   
      fgt (n_skin_atoms, 8 * B * B * B, HS_coordinates, grid, affinity_function_buffer, 1.1 / special_sqrt, HS_gamma, 1e-15);
   
      grid_point_index = 0;
   
      for ( int i = 0; i < 2 * B; ++i ) 
      {
         for ( int j = 0; j < 2 * B; ++j )
         {
            for ( int k = 0; k < 2 * B; ++k )
            {
               affinity_function_real[grid_point_index] += affinity_function_buffer[grid_point_index] * exp_array[i];
            
               ++grid_point_index;
            }
         }
      }
   
      delete[] HS_gamma; 
   }
   
   delete[] r;
   delete[] grid;
   delete[] affinity_function_buffer;
   delete[] exp_array;
   
   printf ("done.\n\n");
   
   return;  
}

void ligand::compute_SGL_Fourier_coefficients (int B, std::string path)
{
   printf (">> computing SGL Fourier coefficients of ligand using FSGLFT... ");
   fflush (stdout);
   
   SGL_Fourier_coefficients_real = new double**[B];
   SGL_Fourier_coefficients_imag = new double**[B];
   
   for ( int n = 1; n <= B; ++n )
   {
      SGL_Fourier_coefficients_real[n - 1] = new double*[n];
      SGL_Fourier_coefficients_imag[n - 1] = new double*[n];
   
      for ( int l = 0; l < n; ++l )
      {
         SGL_Fourier_coefficients_real[n - 1][l] = new double[2 * l + 1];
         SGL_Fourier_coefficients_imag[n - 1][l] = new double[2 * l + 1];
      }
   }
   
   fsglft (B, affinity_function_real, affinity_function_imag, SGL_Fourier_coefficients_real, SGL_Fourier_coefficients_imag, path);
   
   delete[] affinity_function_real;
   delete[] affinity_function_imag;
   
   printf ("done.\n\n");   
   
   return;
}

#endif