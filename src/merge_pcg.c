/* 
* Copyright (c) 2011-2012, 2014 Genome Research Ltd.
*
* Author: Zemin Ning <zn1@sanger.ac.uk>, Andrew Whitwham <aw7@sanger.ac.uk>
*
* This file is part of iCAS.
*
* iCAS is free software: you can redistribute it and/or modify it under
* the terms of the GNU General Public License as published by the Free Software
* Foundation; either version 3 of the License, or (at your option) any later
* version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public License along with
* this program. If not, see <http://www.gnu.org/licenses/>. 
*/


/*
* merge_pcg.c - Code modified for Illumina Clone Assembly Pipeline
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "array_sort.h"
#include "matrix.h"

static char **S_Name, **R_Name, **T_Name;
static int *hit_sindex, *hit_rcdex, *hit_locus1, *hit_locus2, *readlength, *superlength;
static int *hit_read1, *hit_read2, *hit_length;
static int *hit_qindex;

/* forward declarations */
void Indel_Process(int nSeq);

int main(int argc, char **argv)
{
    FILE *namef;
    int i;
    int n_reads, nseq = 0;
    float identy;
    char line[2000] = {0}, tempc1[60], RC[1];
    char *filename;
    
    filename = argv[1];

    fflush(stdout);

    if (argc < 2) {
      	printf("Usage: %s <cross_genome_ouput file>\n", argv[0]);
    	exit(1);
    }

    if ((namef = fopen(filename, "r")) == NULL) {
      	fprintf(stderr, "ERROR: unable to open file %s\n", filename);
      	exit(1);
    }

    while (!feof(namef)) {
    	fgets(line, 2000, namef);
	
      	if (feof(namef)) break;
	
      	nseq++;
    }
    
    fclose(namef); 

    if((hit_sindex = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      fprintf(stderr, "Error: out of memory: calloc - hit_sindex\n");
      exit(1);
    }
    
    if ((hit_rcdex = (int *)calloc(nseq, sizeof(int))) == NULL) {
      	fprintf(stderr, "Error: out of memory: calloc - hit_rcdex\n");
      	exit(1);
    }
    
    if ((hit_read1 = (int *)calloc(nseq, sizeof(int))) == NULL) {
      	fprintf(stderr, "Error: out of memory: calloc - hit_read1\n");
      	exit(1);
    }
    
    if ((hit_read2 = (int *)calloc(nseq, sizeof(int))) == NULL) {
      	fprintf(stderr, "Error: out of memory: calloc - hit_read2\n");
      	exit(1);
    }
    
    if ((hit_locus1 = (int *)calloc(nseq, sizeof(int))) == NULL) {
      	fprintf(stderr, "Error: out of memoryfmate: calloc - hit_locus1\n");
      	exit(1);
    }
    
    if ((hit_locus2 = (int *)calloc(nseq, sizeof(int))) == NULL) {
      	fprintf(stderr, "Error: out of memory: calloc - hit_locus2\n");
      	exit(1);
    }
    
    if ((hit_length = (int *)calloc(nseq, sizeof(int))) == NULL) {
      	fprintf(stderr, "Error: calloc - hit_length\n");
      	exit(1);
    }
    
    if ((readlength = (int *)calloc(nseq, sizeof(int))) == NULL) {
      	fprintf(stderr, "Error: calloc - readlength\n");
      	exit(1);
    }
    
    if ((superlength = (int *)calloc(nseq, sizeof(int))) == NULL) {
      	fprintf(stderr, "Error: calloc - superlength\n");
      	exit(1);
    }
    
    if ((hit_qindex = (int *)calloc(nseq, sizeof(int))) == NULL) {
      	fprintf(stderr, "Error: calloc - hit_qindex\n");
      	exit(1);
    }

    R_Name = cmatrix(0, nseq + 1, 0, Max_N_NameBase);
    S_Name = cmatrix(0, nseq + 1, 0, Max_N_NameBase);
    T_Name = cmatrix(0, nseq + 1, 0, 6);

    if ((namef = fopen(filename, "r")) == NULL) {
    	fprintf(stderr, "ERROR: unable to open file %s\n", filename);
      	exit(1);
    }

    /*  read the alignment files         */
    i = 0;
    
    while (fscanf(namef, "%s %d %s %s %d %d %d %d %s %d %f %d %d", tempc1, &hit_sindex[i],
    	    R_Name[i], S_Name[i], &hit_read1[i], &hit_read2[i], &hit_locus1[i],
	    &hit_locus2[i], RC, &hit_qindex[i], &identy, &readlength[i],
	    &superlength[i]) != EOF) {
	    
        if (RC[0] == 'F')
            hit_rcdex[i] = 0;
        else
            hit_rcdex[i] = 1;
	
	strncpy(T_Name[i], S_Name[i], 4);
        i++;
    }
    
    fclose(namef);

    n_reads = i;
    printf("reads: %d %s\n", n_reads, filename);
    Indel_Process(n_reads);

    printf("Job finished for %d reads!\n", nseq);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void Indel_Process(int nSeq)
/* =============================== */
{
    int i, j, k, m, n;
    int num_hits, num_hit1, num_hit2;
    int stopflag;
    float rate;
          
    num_hits = 0;
    k = 0;
     
    for (i = 0; i < nSeq; i++) {
    	stopflag = 0;
        j = i + 1;
	
        while ((j < nSeq) && (stopflag == 0)) {
            if (strcmp(R_Name[i], R_Name[j]) == 0) 
            	j++;
            else
            	stopflag = 1;
        }
	
	num_hit1 = j - i;
	
        if ((j - i) >= 2) {
            int stopflag2 = 0;
            num_hits = 0;
	    num_hit2 = 0;
	    {
	    	int num_blocks = 0;
		
        	for (n = i; n < j; n++) {
        	    m = n + 1;
        	    stopflag2 = 0;
		    
        	    while ((m < (j)) && (stopflag2 == 0)) {
                     	if ((strcmp(S_Name[n], S_Name[m]) == 0) && (m < (j)))
                            m++;
                     	else
                        stopflag2 = 1;
        	    }
		    
		    num_hit2 = m - n;
		    
		    if (num_hit2 == num_hit1)
                    	printf("merge-bq: %s %d %d %d %d %d %d %d %d\n", R_Name[i],
			    hit_qindex[i], 1, hit_sindex[i], hit_rcdex[i], hit_locus1[i],
			    hit_locus2[j - 1], readlength[i],superlength[i]); 
        	    else {
                    	int align_len;
		       
		        if (num_hit2 == 1)
		    	    align_len = abs(hit_locus2[m - 1] - hit_locus1[m - 1]);
		        else
			    align_len = abs(hit_locus2[m - 1] - hit_locus1[n]);
			    
		        rate = align_len;
		        rate = rate / superlength[n];
		       
                        if (hit_rcdex[n] == 0) {
			    if (num_blocks == 0) {
                	    	if (hit_read1[n] >= 200)
                	     	    printf("merge-1s: %s %d %d %d %d %d %d %d %d\n", R_Name[n], hit_qindex[n], 0, hit_qindex[n], 0, 1, hit_read1[n], readlength[n], superlength[n]); 
				    
                	    	printf("merge-1q: %s %d %d %d %d %d %d %d %d\n", R_Name[n], hit_qindex[n], 1, hit_sindex[n], hit_rcdex[n], hit_locus1[n], hit_locus2[m - 1], readlength[n], superlength[n]); 
                	    } else {
                	    	printf("merge-2s: %s %d %d %d %d %d %d %d %d\n", R_Name[n], hit_qindex[n], 0, hit_qindex[n], 0, hit_read2[n - 1], hit_read1[n], readlength[n], superlength[n]); 
                	    	printf("merge-2q: %s %d %d %d %d %d %d %d %d\n", R_Name[n], hit_qindex[n], 1, hit_sindex[n], hit_rcdex[n], hit_locus1[n], hit_locus2[m - 1], readlength[n], superlength[n]); 
                	    }
                        } else {
			    if (num_blocks == 0) {
                	    	if (hit_read1[n] >= 200)
                	     	    printf("merge-1s: %s %d %d %d %d %d %d %d %d\n", R_Name[n], hit_qindex[n], 0, hit_qindex[n], 0, 1, hit_read1[n], readlength[n], superlength[n]);
				     
                	    	printf("merge-1q: %s %d %d %d %d %d %d %d %d\n", R_Name[n], hit_qindex[n], 1, hit_sindex[n], hit_rcdex[n], hit_locus1[n], hit_locus2[m - 1], readlength[n], superlength[n]); 
                	    } else {
                	    	printf("merge-2s: %s %d %d %d %d %d %d %d %d\n", R_Name[n], hit_qindex[n], 0, hit_qindex[n], 0, hit_read2[n - 1], hit_read1[n], readlength[n], superlength[n]); 
                	    	printf("merge-2q: %s %d %d %d %d %d %d %d %d\n", R_Name[n], hit_qindex[n],1, hit_sindex[n], hit_rcdex[n], hit_locus1[n], hit_locus2[m - 1], readlength[n], superlength[n]); 
                	    }
                    	}
        	    }
		   
        	    num_blocks++;
        	    n = m - 1;
        	}
	    }
        } else {
            printf("merge-sq: %s %d %d %d %d %d %d %d %d\n", R_Name[i], hit_qindex[i], 1, hit_sindex[i], hit_rcdex[i], hit_locus1[i], hit_locus2[i], readlength[i], superlength[i]); 
        }
	
        i = j - 1;
    }

}
