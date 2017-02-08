#include <fgt.h>
#include <figtree-master/include/figtree.h>

void fgt (int n_source, int n_target, double* source, double* target, double* f, double h, double* q, double epsilon)
{
   figtree (3, n_source, n_target, 1, source, h, q, target, epsilon, f);

   return;
}

// use it like this:
//
//   for ( int k = 0; k < n_atoms; ++k )
//   {
//	     source[3 * k + 0] = atom_coordinates[k][0];
//	     source[3 * k + 1] = atom_coordinates[k][1];
//	     source[3 * k + 2] = atom_coordinates[k][2];
//   }
//   
//   for ( int x_index = 0; x_index < n_points; ++x_index )
//   {   
//	     target[3 * x_index + 0] = x[x_index][0];
//	     target[3 * x_index + 1] = x[x_index][1];
//	     target[3 * x_index + 2] = x[x_index][2];
//   }