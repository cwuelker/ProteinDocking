#ifndef protein_h
#define protein_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

/*
 * protein class
 */
class protein
{
public:

   protein (std::string filename);
   
   virtual ~protein () {}
   
   void create_skin (int N);
   void evaluate_affinity_function (int B, double lambda, std::string path);
   void compute_SGL_Fourier_coefficients (int B, std::string path);
   
   double get_SGL_Fourier_coefficient_real (int n, int l, int m) { return SGL_Fourier_coefficients_real[n - 1][l][m + l]; };
   double get_SGL_Fourier_coefficient_imag (int n, int l, int m) { return SGL_Fourier_coefficients_imag[n - 1][l][m + l]; };

private:

   double* coordinates;
   char* type;
   
   double* H_coordinates, * N_coordinates, * C_coordinates, * O_coordinates, * P_coordinates, * F_coordinates, * S_coordinates, * unknown_atom_coordinates;
   
   double* H_gamma, * N_gamma, * C_gamma, * O_gamma, * P_gamma, * F_gamma, * S_gamma, * unknown_atom_gamma;
   
   double* HS_coordinates, * NS_coordinates, * CS_coordinates, * OS_coordinates, * PS_coordinates, * FS_coordinates, * SS_coordinates, * US_coordinates;

   int n_atoms, n_skin_atoms, n_H_atoms, n_N_atoms, n_C_atoms, n_O_atoms, n_P_atoms, n_F_atoms, n_S_atoms, n_unknown_atoms, n_HS_atoms, n_NS_atoms, n_CS_atoms, n_OS_atoms, n_PS_atoms, n_FS_atoms, n_SS_atoms, n_US_atoms;
   
   double* affinity_function_real, * affinity_function_imag;
   
   double*** SGL_Fourier_coefficients_real, *** SGL_Fourier_coefficients_imag;
};

/*
 * constructor
 */
protein::protein (std::string filename)
{
   printf (">> loading protein from pdb file \"%s\"... ", filename.c_str());
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
   in.open (filename.c_str());
   
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
         default: ++n_unknown_atoms;
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
   
   H_gamma = new double[n_H_atoms];
   N_gamma = new double[n_N_atoms];
   C_gamma = new double[n_C_atoms];
   O_gamma = new double[n_O_atoms];
   S_gamma = new double[n_S_atoms];
   unknown_atom_gamma = new double[n_unknown_atoms];
   
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
            H_gamma[H_counter] = 1.0;
			   ++H_counter;
		      break;
	      case 'N' : 
            N_coordinates[3 * N_counter + 0] = coordinates[3 * index + 0];
			   N_coordinates[3 * N_counter + 1] = coordinates[3 * index + 1];
            N_coordinates[3 * N_counter + 2] = coordinates[3 * index + 2];
            N_gamma[N_counter] = 1.0;
			   ++N_counter;
		      break;
         case 'C' : 
            C_coordinates[3 * C_counter + 0] = coordinates[3 * index + 0];
			   C_coordinates[3 * C_counter + 1] = coordinates[3 * index + 1];
			   C_coordinates[3 * C_counter + 2] = coordinates[3 * index + 2];
            C_gamma[C_counter] = 1.0;
            ++C_counter;
		      break;
		   case 'O' : 
		      O_coordinates[3 * O_counter + 0] = coordinates[3 * index + 0];
			   O_coordinates[3 * O_counter + 1] = coordinates[3 * index + 1];
			   O_coordinates[3 * O_counter + 2] = coordinates[3 * index + 2];
            O_gamma[O_counter] = 1.0;
            ++O_counter;
		      break;
		   case 'S' : 
		      S_coordinates[3 * S_counter + 0] = coordinates[3 * index + 0];
			   S_coordinates[3 * S_counter + 1] = coordinates[3 * index + 1];
			   S_coordinates[3 * S_counter + 2] = coordinates[3 * index + 2];
            S_gamma[S_counter] = 1.0;
            ++S_counter;
		      break;
		   default : 
		      unknown_atom_coordinates[3 * unknown_atom_counter + 0] = coordinates[3 * index + 0];
			   unknown_atom_coordinates[3 * unknown_atom_counter + 1] = coordinates[3 * index + 1];
			   unknown_atom_coordinates[3 * unknown_atom_counter + 2] = coordinates[3 * index + 2];
            unknown_atom_gamma[unknown_atom_counter] = 1.0;
            ++unknown_atom_counter;		 
		      break;
	  }
   }

   in.close ();

   printf ("done.\n");
   printf ("   -> number of atoms: %d (%d H, %d N, %d C, %d O, %d S, %d unknown)\n\n", n_atoms, n_H_atoms, n_N_atoms, n_C_atoms, n_O_atoms, n_S_atoms, n_unknown_atoms);
}

void protein::create_skin (int N)
{
   printf (">> detecting skin of protein... ");
   fflush (stdout);
   
   std::vector<double> skin_x;
   std::vector<double> skin_y;
   std::vector<double> skin_z;
   
   std::vector<char> skin_type;
   
   int H_index = -1;
   int N_index = -1;
   int C_index = -1;
   int O_index = -1;
   int S_index = -1;
   int unknown_atom_index = -1;   
   
   for ( int atom_index = 0; atom_index < n_atoms; ++atom_index )
   {
      double r;
   
      bool already_added = false;
      
      switch ( type[atom_index] )
      {
         case 'O' : r = 1.52; ++O_index;
            break;
         case 'H' : r = 1.1; ++H_index;
            break;
         case 'C' : r = 1.7; ++C_index;
            break;
         case 'N' : r = 1.55; ++N_index;
            break;
         case 'S' : r = 1.8; ++S_index;
            break;
         default : r = 1.0; ++unknown_atom_index;
            break;
      }
      
      double a = 4 * pi * r * r / N;
      
      double d = sqrt (a);
      
      int M_theta = static_cast<int> (pi / d + 0.5);
      
      double d_theta = pi / M_theta;
      double d_phi = a / d_theta;
      
      for ( int m = 0; m < M_theta && !already_added; ++m )
      {
         double theta = pi * (m + 0.5) / M_theta;
         
         double M_phi = static_cast<int> (2 * pi * sin (theta) / d_phi + 0.5);
         
         for ( int n = 0; n < M_phi && !already_added; ++n )
         {
            double phi = 2 * pi  * n / M_phi;
            
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
               
               double distance_squared; 
               
               distance_squared = (x - coordinates[3 * atom_index_2 + 0]) * (x - coordinates[3 * atom_index_2 + 0]);
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
               skin_x.push_back (coordinates[3 * atom_index + 0]);
               skin_y.push_back (coordinates[3 * atom_index + 1]);
               skin_z.push_back (coordinates[3 * atom_index + 2]);
            
               skin_type.push_back (type[atom_index]);
            
               switch ( type[atom_index] )
               {
                  case 'O' : ++n_OS_atoms; O_gamma[O_index] = 0.0;
                     break;
                  case 'H' : ++n_HS_atoms; H_gamma[H_index] = 0.0;
                     break;
                  case 'C' : ++n_CS_atoms; C_gamma[C_index] = 0.0;
                     break;
                  case 'N' : ++n_NS_atoms; N_gamma[N_index] = 0.0;
                     break;
                  case 'S' : ++n_SS_atoms; S_gamma[S_index] = 0.0;
                     break;
                  default : ++n_US_atoms; unknown_atom_gamma[unknown_atom_index] = 0.0;
                     break;
               }
            
               already_added = true;
            }
         }
      }
   }
   
   n_skin_atoms = (int)skin_type.size ();
   
   HS_coordinates = new double[3 * n_HS_atoms];
   OS_coordinates = new double[3 * n_OS_atoms];
   CS_coordinates = new double[3 * n_CS_atoms];
   NS_coordinates = new double[3 * n_NS_atoms];
   SS_coordinates = new double[3 * n_SS_atoms];
   US_coordinates = new double[3 * n_US_atoms];
   
   int HS_counter = 0;
   int NS_counter = 0;
   int CS_counter = 0;
   int OS_counter = 0;
   int SS_counter = 0;
   int unknown_skin_atom_counter = 0;   
   
   for ( int skin_atom_index = 0; skin_atom_index < n_skin_atoms; ++skin_atom_index )
   {
      switch ( skin_type.at (skin_atom_index) )
      {
         case 'H' : 
            HS_coordinates[3 * HS_counter + 0] = skin_x.at (skin_atom_index);
            HS_coordinates[3 * HS_counter + 1] = skin_y.at (skin_atom_index);
            HS_coordinates[3 * HS_counter + 2] = skin_z.at (skin_atom_index);
            ++HS_counter;
            break;
         case 'N' : 
            NS_coordinates[3 * NS_counter + 0] = skin_x.at (skin_atom_index);
            NS_coordinates[3 * NS_counter + 1] = skin_y.at (skin_atom_index);
            NS_coordinates[3 * NS_counter + 2] = skin_z.at (skin_atom_index);
            ++NS_counter;
            break;
         case 'C' : 
            CS_coordinates[3 * CS_counter + 0] = skin_x.at (skin_atom_index);
            CS_coordinates[3 * CS_counter + 1] = skin_y.at (skin_atom_index);
            CS_coordinates[3 * CS_counter + 2] = skin_z.at (skin_atom_index);
            ++CS_counter;
            break;
         case 'O' : 
            OS_coordinates[3 * OS_counter + 0] = skin_x.at (skin_atom_index);
            OS_coordinates[3 * OS_counter + 1] = skin_y.at (skin_atom_index);
            OS_coordinates[3 * OS_counter + 2] = skin_z.at (skin_atom_index);
            ++OS_counter;
            break;
         case 'S' : 
            SS_coordinates[3 * SS_counter + 0] = skin_x.at (skin_atom_index);
            SS_coordinates[3 * SS_counter + 1] = skin_y.at (skin_atom_index);
            SS_coordinates[3 * SS_counter + 2] = skin_z.at (skin_atom_index);
            ++SS_counter;
            break;
         default : 
            US_coordinates[3 * unknown_skin_atom_counter + 0] = skin_x.at (skin_atom_index);
            US_coordinates[3 * unknown_skin_atom_counter + 1] = skin_y.at (skin_atom_index);
            US_coordinates[3 * unknown_skin_atom_counter + 2] = skin_z.at (skin_atom_index);
            ++unknown_skin_atom_counter;		 
            break;
      }
   }
   
   skin_x.clear ();
   skin_y.clear ();
   skin_z.clear ();
   
   skin_type.clear ();
   
   printf ("done.\n");
   printf ("   -> number of skin atoms: %d (%d H, %d N, %d C, %d O, %d S, %d unknown)\n\n", n_skin_atoms, n_HS_atoms, n_NS_atoms, n_CS_atoms, n_OS_atoms, n_SS_atoms, n_US_atoms);
   
   return;
}

void protein::evaluate_affinity_function (int B, double lambda, std::string path)
{
   printf (">> evaluating affinity function for protein using FGT... ");
   fflush (stdout);
   
   double sqrt_lambda = sqrt (lambda);
   double lambda_pow = pow (lambda, 0.75);
   
   affinity_function_real = new double[8 * B * B * B];
   affinity_function_imag = new double[8 * B * B * B];
   
   for ( int index = 0; index < 8 * B * B * B; ++index )
   {
      affinity_function_real[index] = 0.0;   
      affinity_function_imag[index] = 0.0;
   }
   
   double* affinity_function_buffer = new double[8 * B * B * B];
   
   double* grid = new double[3 * 8 * B * B * B];
   
   double special_sqrt = sqrt (2.31);
   double special_exp = exp (2.31);
   
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
   
   delete[] filename;
   
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
   
   if ( n_H_atoms > 0 )
   {      
      for ( int index = 0; index < n_H_atoms; ++index )
         H_gamma[index] *= 16.0 * lambda_pow * special_exp;
      
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
   
   if ( n_N_atoms )
   {      
      for ( int index = 0; index < n_N_atoms; ++index )
         N_gamma[index] *= 16.0 * lambda_pow * special_exp;
      
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
      for ( int index = 0; index < n_O_atoms; ++index )
         O_gamma[index] *= 16.0 * lambda_pow * special_exp;
      
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
   
   if ( n_C_atoms )
   {            
      for ( int index = 0; index < n_C_atoms; ++index )
         C_gamma[index] *= 16.0 * lambda_pow * special_exp;
         
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
      for ( int index = 0; index < n_S_atoms; ++index )
         S_gamma[index] *= 16.0 * lambda_pow * special_exp;
      
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
      for ( int index = 0; index < n_unknown_atoms; ++index )
         unknown_atom_gamma[index] *= 16.0 * lambda_pow * special_exp;
      
      fgt (n_unknown_atoms, 8 * B * B * B, unknown_atom_coordinates, grid, affinity_function_buffer, 1 / special_sqrt, unknown_atom_gamma, 1e-15);
      
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
      
      delete[] unknown_atom_gamma;
   }
   
   // --------------- SKIN
   
   if ( n_HS_atoms > 0 )
   {      
      double* HS_gamma = new double[n_HS_atoms];
      
      for ( int index = 0; index < n_HS_atoms; ++index )
         HS_gamma[index] = lambda_pow * special_exp;
         
         fgt (n_HS_atoms, 8 * B * B * B, HS_coordinates, grid, affinity_function_buffer, 1.1 / special_sqrt, HS_gamma, 1e-15);
         
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
   
   // ---------------
   
   if ( n_NS_atoms > 0 )
   {      
      double* NS_gamma = new double[n_NS_atoms];
      
      for ( int index = 0; index < n_NS_atoms; ++index )
         NS_gamma[index] = lambda_pow * special_exp;
         
      fgt (n_NS_atoms, 8 * B * B * B, NS_coordinates, grid, affinity_function_buffer, 1.55 / special_sqrt, NS_gamma, 1e-15);
         
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
         
      delete[] NS_gamma;
   }
   
   // ---------------
   
   if ( n_OS_atoms > 0 )
   {      
      double* OS_gamma = new double[n_OS_atoms];
      
      for ( int index = 0; index < n_OS_atoms; ++index )
         OS_gamma[index] = lambda_pow * special_exp;
         
      fgt (n_OS_atoms, 8 * B * B * B, OS_coordinates, grid, affinity_function_buffer, 1.52 / special_sqrt, OS_gamma, 1e-15);
         
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
         
      delete[] OS_gamma;   
   }
   
   // ---------------
   
   if ( n_CS_atoms > 0 )
   {            
      double* CS_gamma = new double[n_CS_atoms];
      
      for ( int index = 0; index < n_CS_atoms; ++index )
         CS_gamma[index] = lambda_pow * special_exp;
         
      fgt (n_CS_atoms, 8 * B * B * B, CS_coordinates, grid, affinity_function_buffer, 1.7 / special_sqrt, CS_gamma, 1e-15);
         
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
         
      delete[] CS_gamma;     
   }
   
   // ---------------
   
   if ( n_SS_atoms > 0 )
   {
      double* SS_gamma = new double[n_SS_atoms];
   
      for ( int index = 0; index < n_SS_atoms; ++index )
         SS_gamma[index] = lambda_pow * special_exp;
      
      fgt (n_SS_atoms, 8 * B * B * B, SS_coordinates, grid, affinity_function_buffer, 1.8 / special_sqrt, SS_gamma, 1e-15);
      
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
      
      delete[] SS_gamma;      
   }
   
   // ---------------
   
   if ( n_US_atoms > 0 )
   {      
      double* US_gamma = new double[n_US_atoms];
      
      for ( int index = 0; index < n_US_atoms; ++index )
         US_gamma[index] = lambda_pow * special_exp;
         
      fgt (n_US_atoms, 8 * B * B * B, US_coordinates, grid, affinity_function_buffer, 1 / special_sqrt, US_gamma, 1e-15);
         
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
         
      delete[] US_gamma;
   }
   
   delete[] r;
   delete[] grid;
   delete[] affinity_function_buffer;
   delete[] exp_array;
   
   printf ("done.\n\n");
   
   return;  
}

void protein::compute_SGL_Fourier_coefficients (int B, std::string path)
{
   printf (">> computing SGL Fourier coefficients of protein using FSGLFT... ");
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