/* Written by Ruibin Xi March 25 2010

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

int cmp_int(const void *a, const void *b)
        { int tmp1,tmp2;
          tmp1 = *((const int *) a);
          tmp2 = *((const int *) b);
          if(tmp1<tmp2) return -1;
          else if(tmp1>tmp2) return 1;
          else return 0;
        }



void binning(int *tum_reads, int *norm_reads, int *n_tum, int *n_norm,int *bn_size,char **output_file)
	{ FILE *output;
	  int rd_in_bin_tum, rd_in_bin_norm,sum=0, n_tum_reads, n_norm_reads, bin_size;
	  int i,j,k, ind_tum, ind_norm; /*ind_tum and ind_norm records which reads is under consideration*/
	  int flag = 0; /*flag whether or not all reads have been processed*/

	  n_tum_reads = *n_tum;
	  n_norm_reads = *n_norm;
	  bin_size = *bn_size;

	  if(bin_size<=0)
		{ printf("Error, the bin size must be positive interger.\n");
		  return;
		}

	  output = fopen(output_file[0],"w"); /*open for output*/
	  if(output==NULL)
		{ printf("cannot open the file %s, reason \'No such file or directory\'\n",output_file[0]);
		  return;
		}

	/*sort the reads*/
	qsort(tum_reads,n_tum_reads,sizeof(int),cmp_int);
	printf("sorted %d tumor reads\n",n_tum_reads);
	qsort(norm_reads,n_norm_reads,sizeof(int),cmp_int);
	printf("sorted %d normal reads\n",n_norm_reads);



	  k = 0;
	  ind_tum = 0;
	  ind_norm  = 0;
	  while(flag==0)
		{ sum = 0;
		  rd_in_bin_tum = 0;
		  rd_in_bin_norm = 0;
		 
		  flag = 1;
		  
		  /*count the tumor reads in the kth bin*/
		  j = ind_tum;
		  while(tum_reads[j]<=(k+1)*bin_size&&j<n_tum_reads)
			{  rd_in_bin_tum ++;
			   sum++;
			   j++;
			}
		  ind_tum = j;
		  flag = flag*(j>=n_tum_reads);

		  j = ind_norm;
		  while(norm_reads[j]<=(k+1)*bin_size&&j<n_norm_reads)
			{ rd_in_bin_norm ++;
			  sum ++;
			  j++;
			}	
		  ind_norm = j;
		  flag = flag*(j>=n_norm_reads);

		  if(sum>0) /*at least one read in the kth bin, print the bin*/
                        { /*print tumor and its matched normal*/ 
			  fprintf(output,"%d\t%d\t",rd_in_bin_tum,sum);
			  fprintf(output,"%.7f\t",((double) rd_in_bin_tum)/((double) sum));
			  fprintf(output,"%d\t%d",k*bin_size+1,(k+1)*bin_size);
                          fprintf(output,"\n");
                        }
		  k ++; /*next bin*/
		}
	 fclose(output);
	 return;
     	}
