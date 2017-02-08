#include <stdio.h>

#include <fsft.h>

extern "C"
{
   #include <SpharmonicKit27/FST_semi_memo.h>
   #include <SpharmonicKit27/cospmls.h>
}
 
void fsft (int B, double* rdata, double* idata, double** rcoeffs, double** icoeffs)
{   
   const int spharmonic_tablesize = Reduced_Naive_TableSize (B, 0) + Reduced_SpharmonicTableSize (B, 0);
   
   double* resultspace = new double[spharmonic_tablesize];
   double* workspace = new double[16 * B];
   double* workspace_2 = new double[(8 * B * B) + (29 * B)];
   
   double* rcoeffs_compressed = new double[B * B];
   double* icoeffs_compressed = new double[B * B];
   
   double** seminaive_pml_table = SemiNaive_Naive_Pml_Table (B, B, resultspace, workspace); 
   
   FST_semi_memo (rdata, idata, rcoeffs_compressed, icoeffs_compressed, 2 * B, seminaive_pml_table, workspace_2, 0, B);

   for ( int l = 0; l < B; ++l )
   {       
      for ( int m = - l; m <= l; ++m )
      {        
         if ( m < 0 )
         {
            int offset = 0;
         
            for ( unsigned index = 0; index < - m; ++index )
               offset += B - index;

            rcoeffs[l][m + B - 1] = rcoeffs_compressed[l + m + offset];
            icoeffs[l][m + B - 1] = icoeffs_compressed[l + m + offset];
         }
         else if ( m == 0 )
         {
            rcoeffs[l][m + B - 1] = rcoeffs_compressed[l];
            icoeffs[l][m + B - 1] = icoeffs_compressed[l];
         }
         else  
         {
            int offset = 0;
            
            for ( unsigned index = 1; index < m; ++index )
               offset += B - index;
            
            rcoeffs[l][m + B - 1] = rcoeffs_compressed[B * (B - 1) + l - offset];
            icoeffs[l][m + B - 1] = icoeffs_compressed[B * (B - 1) + l - offset];
         }
      }
   }
   
   delete[] rcoeffs_compressed;
   delete[] icoeffs_compressed;
   
   delete[] resultspace;
   delete[] workspace;
   delete[] workspace_2;
   
   // TODO: delete seminaive_pml_table

   return;
}