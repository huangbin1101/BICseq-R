/* Jianhua Zhang. All rights reserved

*/

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void getratios(int *tum_reads, int *norm_reads, int *length, double *adjust, double *ratios){
	int i;
	for (i = 0; i < length[0]; i++) {
	     if(tum_reads[i] == 0 && norm_reads[i] == 0)
			 ratios[i] = 0;
	     else if(norm_reads[i] == 0){
			 if(tum_reads[i] == 1)
			     ratios[i] = (0.465465 - adjust[0]);
			 else
			     ratios[i] = (log2(tum_reads[i]) - adjust[0]);
		 }else if(tum_reads[i] == 0){
			 if(norm_reads[i] == 1)
			     ratios[i] = (-0.465465 - adjust[0]);
			 else
			     ratios[i] = (-log2(norm_reads[i]) - adjust[0]);
		 }else
		     ratios[i] = (log2(tum_reads[i]) - log2(norm_reads[i]) - adjust[0]);
	 }

	 return;
 }


