#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "pos_cnt_lst.h"



static int bin_col = 3; /*number of columns in a bin file*/

static int cmp_integer(const void *a, const void *b);

static int *aggregate(int *reads, int n_reads, int *num_pos);
/* aggregate the reads located at the same position together.
 * Returned value is a two column matrix looking like
 * 
 * <position>	<number of read at this position>
 *
 * The argument num_pos will record the size of this matrix (number of rows)
 *
 * reads: the read positions which should be ordered nondecreasingly
 * n_reads: total number of reads (length of the vector reads)
 * */
static void binning_new(int *tumor_1bp_bin, int n_tumor, int *normal_1bp_bin, int n_normal,int bin_size, FILE *output);
/* tumor_1bp_bin: output of the function aggregate for tumor reads
 * n_tumor: nrows of tumor_1bp_bin
 * similar for normal_1bp_bin and n_normal
 * bin_size: the bin size
 * nbins: total number of bins created. 
 *
 * */

void sort_rms_binning(int *tumor,int *n_tmor,int *normal, int *n_nml,int *bin_size, int *w, double *quantile, double *multple, char **output_file);
/* sort the reads, remove the singular reads and bin the data.
   arguments:
	nbins: the total number of bins obtained.
	w: the window size used to determine the singular genomic positions.
	quantile: quantile (e.g. 0.95) used to determine the singular genomic positions.
	multple: if a genomic position has more than multple * quantile, it will be identified as an outlier.
 */



static int cmp_integer(const void *a, const void *b)
	{ double tmp1,tmp2;
	  tmp1 = *((const int *) a);
	  tmp2 = *((const int *) b);
	  if(tmp1<tmp2) return -1;
	  else if(tmp1>tmp2) return 1;
	  else return 0;
	}




static int *aggregate(int *reads, int n_reads, int *num_pos)
	{ /*The first column of read_dist is the position of the read*/
	  /*The second column of read_list is the number of reads at this position.*/
	  /*If certain position does not have any reads, this position will be ignored*/
	  int *reads_dist, pos_min, pos_max,size;
	  int i,j;
	
	  pos_min = reads[0];
	  pos_max = reads[n_reads-1];

	  size = (n_reads<pos_max-pos_min+1)? n_reads : pos_max-pos_min+1;
	  reads_dist = (int *)  malloc(sizeof(int)*(size*2+10));

	  j = 0;
	  i =0;
	  reads_dist[2*j] = reads[i]; /*The first column of read_dist is the position of the read*/
	  reads_dist[2*j+1] = 1; /*The second column of read_list is the number of reads at this position.*/

	  for(i=1;i<n_reads;i++)
		{ if(reads[i]==reads_dist[2*j]) 
			{ reads_dist[2*j+1]++;} /*one more read at position reads_list[2*j]*/
		  else /*get to a new position*/
			{ j++;
			  reads_dist[2*j] = reads[i];
			  reads_dist[2*j+1] = 1;
			}
		}

	 *num_pos = j+1;
	
	 return reads_dist;

	}


static void binning_new(int *tumor_1bp_bin, int n_tumor, int *normal_1bp_bin, int n_normal,int bin_size, FILE *output)
	{ int i,nrow,max_pos,min_pos, *bin;
	  int start,end ,i_tum, i_norm;
	  int total=0;
	  double prob;

	  max_pos = (tumor_1bp_bin[2*(n_tumor-1)]>normal_1bp_bin[2*(n_normal-1)])? tumor_1bp_bin[2*(n_tumor-1)]: normal_1bp_bin[2*(n_normal-1)];
	  min_pos = (tumor_1bp_bin[0]<normal_1bp_bin[0])?tumor_1bp_bin[0]:normal_1bp_bin[0];
	  nrow = (max_pos-min_pos+1)/bin_size+10;
	  bin = (int *) malloc(sizeof(int)*10);

	  i_tum = 0;
	  i_norm = 0;
	  start = min_pos - (min_pos-1)%bin_size;/*the start position of the left most bin*/
	  end = start + bin_size -1;
	  for(i=0;i<nrow;i++)
		{ bin[0] = start;
		  bin[1] = 0;
		  bin[2] = 0;
		  while(tumor_1bp_bin[2*i_tum]<=end&&i_tum<n_tumor)
			{ bin[1] += tumor_1bp_bin[2*i_tum+1];
			  i_tum++;
			}
		  while(normal_1bp_bin[2*i_norm]<=end&&i_norm<n_normal)
			{ bin[2] += normal_1bp_bin[2*i_norm+1];
			  i_norm++;
			}
		  
		  total = bin[1]+bin[2];
		  end = start+bin_size -1;
		  if(total >0){
	  		  prob = ((double) bin[1])/((double) total);
	 		  fprintf(output,"%d\t",bin[1]);
			  fprintf(output,"%d\t",total);
			  fprintf(output,"%0.7f\t",prob);
			  fprintf(output,"%d\t",start);
			  fprintf(output,"%d\n",end);
			}

		  start += bin_size;
		  end = start + bin_size -1;
		}
	  free(bin);
	  return;
	}


void sort_rms_binning(int *tumor, int *n_tmor,int *normal, int *n_nml,int *bin_size, int *w, double *quantile, double *multple,char **output_file)
	{ int *tumor_1bpbin,num_tum_1bp, *normal_1bpbin,num_norm_1bp;
	  int *bin;
	  FILE *output;

	  output = fopen(output_file[0],"w"); /*open for output*/
          if(output==NULL)
                { fprintf(stderr,"cannot open the file %s, reason \'No such file or directory\'\n",output_file[0]);
                  return;
                }


	  qsort(tumor,n_tmor[0],sizeof(int),cmp_integer);
          printf("sorted %d tumor reads\n",n_tmor[0]);
       	  qsort(normal,n_nml[0],sizeof(int),cmp_integer);
      	  printf("sorted %d normal reads.\n",n_nml[0]);
 
	  /*aggregate the reads; */
	  tumor_1bpbin = aggregate(tumor,n_tmor[0],&num_tum_1bp);
	  normal_1bpbin = aggregate(normal,n_nml[0],&num_norm_1bp);

	  /*remove the singular points*/
	  singularity_rm(tumor_1bpbin, num_tum_1bp, w[0], quantile[0], multple[0]);
	  singularity_rm(normal_1bpbin, num_norm_1bp, w[0], quantile[0], multple[0]);

	  /*bin the processed reads*/
	  binning_new(tumor_1bpbin,num_tum_1bp,normal_1bpbin, num_norm_1bp,bin_size[0], output);
	  free(tumor_1bpbin);
	  free(normal_1bpbin);
	  fclose(output);

	 return;
	}


