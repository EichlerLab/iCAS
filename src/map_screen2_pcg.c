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
* map_screen2_pcg.c - Code modified for Illumina Clone Assembly Pipeline
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "array_sort.h"
#include "matrix.h"

/* forward declarations */
void Indel_Process(char **S_Name, char **R_Name, int *hit_mask, int *superlength, int *hit_read2, int nSeq);

int main(int argc, char **argv) {
    FILE *namef;
    int i, nSeq = 0;
    int n_reads, nseq = 0;
    char line[2000]={0}, tempc1[60], RC[1];
    char *filename;
    int *hit_rcdex, *hit_locus1, *hit_locus2, *readlength;
    int *hit_read1, *hit_length;
    char **S_Name, **R_Name, **T_Name;
    float *hit_identy;
    int *hit_mask, *superlength;
    int *hit_read2;
    
    fflush(stdout);

    if (argc < 2) {
    	printf("Usage: %s <cross_genome_ouput file>\n", argv[0]);
      	exit(1);
    }

    filename = argv[1];

    if ((namef = fopen(filename, "r")) == NULL) {
    	fprintf(stderr, "ERROR: can't open file %s\n", filename);
      	exit(1);
    }
    
    while (!feof(namef)) {
    	fgets(line, 2000, namef);
	
      	if(feof(namef)) break;
	
      	nseq++;
    }
    
    fclose(namef); 

    if ((hit_mask = (int *)calloc(nseq, sizeof(int))) == NULL) {
    	fprintf(stderr, "Error: out of memory: calloc - hit_mask\n");
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
      	fprintf(stderr, "Error: out of memory: calloc - hit_locus1\n");
      	exit(1);
    }
    
    if ((hit_locus2 = (int *)calloc(nseq, sizeof(int))) == NULL) {
      	fprintf(stderr, "Error: out of memory: calloc - hit_locus2\n");
      	exit(1);
    }
    
    if ((hit_length = (int *)calloc(nseq, sizeof(int))) == NULL) {
      	fprintf(stderr, "Error: out of memory: calloc - hit_length\n");
     	exit(1);
    }
    
    if ((readlength = (int *)calloc(nseq, sizeof(int))) == NULL) {
      	fprintf(stderr, "Error: out of memory: calloc - readlength\n");
      	exit(1);
    }
    
    if ((superlength = (int *)calloc(nseq, sizeof(int))) == NULL) {
    	fprintf(stderr, "Error: out of memory: calloc - superlength\n");
      	exit(1);
    }
    
    if ((hit_identy = (float *)calloc(nseq, sizeof(float))) == NULL) {
      	fprintf(stderr, "Error: out of memory: calloc - hit_identy\n");
      	exit(1);
    }

    R_Name = cmatrix(0, nseq + 10, 0, Max_N_NameBase);
    S_Name = cmatrix(0, nseq + 10, 0, Max_N_NameBase);
    T_Name = cmatrix(0, nseq + 10, 0, 6);
    
    printf("reads: %d %s\n", nseq, filename);

    if ((namef = fopen(filename, "r")) == NULL) {
    	fprintf(stderr, "ERROR: unable to open file %s\n", filename);
      	exit(1);
    }

    /*  read the alignment files         */
    i = 0;
    
    while (fscanf(namef, "%s %s %s %s %d %d %d %d %s %d %f %d %d", tempc1, tempc1, R_Name[i], 
    	S_Name[i], &hit_read1[i], &hit_read2[i], &hit_locus1[i], &hit_locus2[i], RC,
	&hit_length[i], &hit_identy[i], &readlength[i], &superlength[i]) != EOF) {
	
    	if (RC[0] == 'F')
            hit_rcdex[i] = 0;
        else {
            int hit1 = hit_locus1[i];
	    int hit2 = hit_locus2[i];
	    hit_locus1[i] = superlength[i] - hit1 +1;
	    hit_locus2[i] = superlength[i] - hit2 +1;
            hit_rcdex[i] = 1;
	}
	
	hit_mask[i] = 0;
	strncpy(T_Name[i], S_Name[i], 4);
        i++;
    }
    
    fclose(namef);

    n_reads = i;
    printf("reads: %d %s\n", n_reads, filename);
    Indel_Process(S_Name, R_Name, hit_mask, superlength, hit_read2, n_reads);
    
    nSeq = nseq;

    nseq = 0;
    
    if ((namef = fopen(filename, "r")) == NULL) {
      	fprintf(stderr, "ERROR: unable to open file %s\n", filename);
      	exit(1);
    }
    
    while (!feof(namef)) {
      	fgets(line, 2000, namef);
	
      	if (feof(namef)) break;
      
      	if (hit_mask[nseq] == 0) printf("%s", line);
      
      	nseq++;
    }
    
    fclose(namef); 

    printf("Job finished for %d reads!\n", nSeq);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void Indel_Process(char **S_Name, char **R_Name, int *hit_mask, int *superlength, int *hit_read2, int nSeq)
/* =============================== */
{
     int i, j, k, m, n;
     int num_hits;
     int stopflag, *readIndex, *readIndex2;
     int offset;
     char **DBname;
          
     if ((readIndex = (int *)calloc(nSeq, sizeof(int))) == NULL) {
        fprintf(stderr, "ERROR Memory_Allocate: calloc - readIndex\n");
        exit(1);
     }
     
     if ((readIndex2 = (int *)calloc(nSeq, sizeof(int))) == NULL) {
        fprintf(stderr, "ERROR Memory_Allocate: calloc - readIndex\n");
        exit(1);
     }
     
     DBname = cmatrix(0, nSeq, 0, Max_N_NameBase);
     
     for (i = 0; i < nSeq; i++) {
     	strcpy(DBname[i], S_Name[i]);
        readIndex[i] = i;
     }
     
     num_hits = 0;
     k        = 0;
     offset   = 0;
     
     for (i = 0; i < nSeq; i++) {
    	stopflag = 0;
        j = i + 1;
	
        while ((j < nSeq) && (stopflag == 0)) {
            if (strcmp(R_Name[i], R_Name[j]) == 0) {
            	j++;
            } else {
            	stopflag = 1;
	    }
        }
	
        if ((j - i) >= 2) {
            int stopflag2 = 0;
	    
	    for(n = i; n < j; n++) {
	     	m = n + 1;
	     	stopflag2 = 0;
		
             	while ((m < j) && (stopflag2 == 0)) {
                    if ((strcmp(S_Name[n], S_Name[m]) == 0) && (m < (j))) {
                    	m++;
                    } else {
                    	stopflag2 = 1;
		    }
	     	}
		
	     	if ((m - n) == 1) {
	            hit_mask[n] = 1;
	     	}
		
	     	n = m - 1;
            }
        } else {
            printf("www: %s %d %d\n", R_Name[i], hit_read2[i], superlength[i]);
        }
	
        i = j - 1;
    }
    
    free(readIndex);
    free(readIndex2);
    free_cmatrix(DBname, 0, nSeq, 0, Max_N_NameBase);
}
