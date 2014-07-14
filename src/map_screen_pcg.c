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
* map_screen_pcg.c - Code modified for Illumina Clone Assembly Pipeline
*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "array_sort.h"
#include "matrix.h"

/* forward declarations */
void Indel_Process(char **S_Name, char **R_Name, int *hit_mask, int *superlength, int nSeq);

int main(int argc, char **argv)
{
    FILE *namef;
    int i;
    int nseq, seq_num;

    char **S_Name, **R_Name, **T_Name;
    int *hit_mask, *superlength;
    int *hit_rcdex, *hit_locus1, *hit_locus2, *readlength;
    int *hit_read1, *hit_read2, *hit_length;
    float *hit_identy;
    char *filename;
    char line[2000] = {0}, tempc1[60], RC[1];
    
    fflush(stdout);

    if (argc < 2) {
    	printf("Usage: %s <cross_genome_ouput file>\n", argv[0]);
      	exit(1);
    }
    
    filename = argv[1];
    nseq = 0;

    if ((namef = fopen(filename, "r")) == NULL) {
    	fprintf(stderr, "ERROR: can't open file %s\n", filename);
      	exit(1);
    }

    while (!feof(namef)) {
    	fgets(line, 2000, namef);
      	if (feof(namef)) break;
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
    
    seq_num = nseq + 10;

    R_Name = cmatrix(0, seq_num, 0, Max_N_NameBase);
    S_Name = cmatrix(0, seq_num, 0, Max_N_NameBase);
    T_Name = cmatrix(0, seq_num, 0, 6);
    
    if ((namef = fopen(filename, "r")) == NULL) {
    	fprintf(stderr, "ERROR: can't open file %s\n", filename);
      	exit(1);
    }

    /* read the alignment files */
    printf("www: %d %s\n", nseq, filename);
    i = 0;
    
    while (fscanf(namef, "%s %s %s %s %d %d %d %d %s %d %f %d %d", tempc1, tempc1, R_Name[i], S_Name[i], &hit_read1[i], &hit_read2[i], &hit_locus1[i], &hit_locus2[i], RC,&hit_length[i], &hit_identy[i], &readlength[i], &superlength[i]) != EOF) {
	if (RC[0] == 'F')
	    hit_rcdex[i] = 0;
	else {
	    int hit1 = hit_locus1[i];
	    int hit2 = hit_locus2[i];
	    hit_locus1[i] = superlength[i] - hit1 + 1;
	    hit_locus2[i] = superlength[i] - hit2 + 1;
	    hit_rcdex[i] = 1;
	}
	
	hit_mask[i] = 0;
	strncpy(T_Name[i], S_Name[i], 4);
	i++;
    }
    
    fclose(namef);

    printf("reads: %d %s\n", i, filename);
    Indel_Process(S_Name, R_Name, hit_mask, superlength, i);

    nseq = 0;
    
    if ((namef = fopen(filename, "r")) == NULL) {
    	fprintf(stderr, "ERROR can't open %s\n", filename);
      	exit(1);
    }
    
    while (!feof(namef)) {
    	fgets(line, 2000, namef);
	
      	if (feof(namef)) break;
	
      	if (hit_mask[nseq] == 0) printf("%s", line);
	
      	nseq++;
    }
    
    fclose(namef);
    
    free(hit_mask);
    free(superlength);
    free(hit_rcdex);
    free(hit_locus1);
    free(hit_locus2);
    free(readlength);
    free(hit_read1);
    free(hit_read2);
    free(hit_length);
    free(hit_identy);
    
    free_cmatrix(S_Name, 0, seq_num, 0, Max_N_NameBase);
    free_cmatrix(R_Name, 0, seq_num, 0, Max_N_NameBase); 
    free_cmatrix(T_Name, 0, seq_num, 0, 6); 
    
    printf("Job finished for %d reads!\n", nseq);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void Indel_Process(char **S_Name, char **R_Name, int *hit_mask, int *superlength, int nSeq)
/* =============================== */
{
    int i, j, k, m, n;
    int num_hits;
    int stopflag, *readIndex, *readIndex2;
    int offset, *dex;
    char **DBname, **ray;

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
    
    for (i = 0; i < (nSeq - 1); i++) {
    	stopflag = 0;
        j = i + 1;
       
        while (( j < nSeq) && (stopflag == 0)) {
	    if(strcmp(R_Name[i], R_Name[j]) == 0) {
	    	j++;
	    } else
	    	stopflag = 1;
        }
       
    	if ((j - i) >= 2) {
    	    int kk, stopflag2 = 0;
	    kk = j - i;
	    ray = DBname;
	    dex = readIndex;
	    
	    ArraySort_String(kk, ray + offset, dex + offset);
	    
	    for (n = i; n < j; n++) {
	        m = n + 1;
	        stopflag2 = 0;
		
	        while ((m < j) && (stopflag2 == 0)) {
		    if ((strcmp(DBname[n], DBname[m]) == 0) && (m < (j))) {
		    	m++;
		    } else
		    	stopflag2 = 1;
	        }
		
	        if ((m - n) <= 3) {
		    for (k = n; k < m;k++) {
		    	int idk = readIndex[k];
			
		    	if (superlength[idk] > 5000) hit_mask[idk] = 1;
		    }
	        }
		
	        n = m - 1;
	    }
    	}

        num_hits = j - i;
        offset = offset + num_hits;
        i = j - 1;
    }
    
    free(readIndex);
    free(readIndex2);
    free_cmatrix(DBname, 0, nSeq, 0, Max_N_NameBase);
}
