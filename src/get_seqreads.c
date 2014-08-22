/* 
* Copyright (c) 2012, 2014 Genome Research Ltd.
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
* get_seqreads.c - Code modified for the Illumina Clone Assembly Pipeline
*/

 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fasta.h"
#include "matrix.h"
#include "array_sort.h"

char **Read_Index(char *filename, int *ctg2wgs_index, fasta *seq, int nRead, int nSeq);

int main(int argc, char **argv) {
    FILE *namef, *fpOutfast;
    long totalBases;
    int i, nSeq, args, n_len;
    int n_contig, n_reads;
    fasta *seq, *seqp;
    char *st, line[2000] = {0};
    long Size_q_pdata;
    char *pdata;
    int fastq_flag = 1;
    int name_flag = 0;
    int *ctg2wgs_index;
    char **rdnames;
    
    if (argc < 2) {
      	printf("Usage: %s <read_name_file> <input_fastq_file> <output_fastq_file\n", argv[0]);
      	exit(1);
    }

    args = 1;

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-fastq")) {
            sscanf(argv[++i], "%d", &fastq_flag); 
            args += 2;
	    
        } else if (!strcmp(argv[i], "-name")) {
            sscanf(argv[++i], "%d", &name_flag);
            args += 2;
        }
    }

    if ((nSeq = extract_fastq_from_file(argv[args + 1], &pdata, &Size_q_pdata)) == 0) {
    	fprintf(stderr, "Error: unable to extract fastq data from %s\n", argv[args + 1]);
	exit(1);
    }

    if ((seq = decode_fastq(pdata, Size_q_pdata, &nSeq, &totalBases)) == NULL) {
    	fprintf(stderr, "Error: no query data found.\n");
	exit(1);
    }

    fastaUC(seq, nSeq);

    if ((namef = fopen(argv[args], "r")) == NULL) {
      	fprintf(stderr, "ERROR unable to open reads file %s\n", argv[args]);
      	exit(1);
    }
    
    n_reads = 0;
    
    while(!feof(namef)) {
      	fgets(line, 2000, namef);
	
      	if (feof(namef)) break;
	
      	n_reads++;
    }
    
    fclose(namef); 

    printf("reads: %d %s\n", n_reads, argv[args+1]);
    
    if ((ctg2wgs_index = (int *)malloc(n_reads * sizeof(int))) == NULL) {
        fprintf(stderr, "ERROR Out of memory : malloc - ctg2wgs_index\n");
        exit(1);
    }
    
    memset(ctg2wgs_index, -1, n_reads*sizeof(int));
    
    rdnames = Read_Index(argv[args], ctg2wgs_index, seq, n_reads, nSeq);

    /*  process contigs one by one   */
    printf("file name  %s \n",argv[args + 2]);
    
    if ((fpOutfast = fopen(argv[args + 2], "w")) == NULL) {
        fprintf(stderr, "ERROR: cannot write to %s\n", argv[args + 2]);
        exit(1);
    }
    
    n_contig = 0;
    
    for (i = 0; i < n_reads; i++) {
        if (ctg2wgs_index[i] >= 0) {
            int nline, k, g, rc;
            int idd = ctg2wgs_index[i];
	    
            n_contig++;
            seqp=seq + idd;
            n_len = seqp->length;

            if (fastq_flag == 0) {
            	fprintf(fpOutfast, ">%s\n", seqp->name);
              	nline=n_len / 60;
              	st = seqp->data;
		
              	for (k = 0; k < nline; k++) {
        	    for (g = 0; g < 60; g++, st++) {
                    	fprintf(fpOutfast, "%c", *st);
		    }
		    
        	    fprintf(fpOutfast, "\n");
              	}
		
              	for (g = 0; g < (n_len - (nline * 60)); g++, st++) {
        	    fprintf(fpOutfast, "%c", *st);
		}
		
              	if ((n_len % 60) != 0) fprintf(fpOutfast, "\n");
		
            } else {
            	if (seqp->name2)
        	    fprintf(fpOutfast, "@%s %s\n", seqp->name, seqp->name2);
              	else
        	    fprintf(fpOutfast, "@%s\n", seqp->name);
		    
              	for (rc = 0; rc < seqp->length; rc++) {
        	    fprintf(fpOutfast, "%c", seqp->data[rc]);
		}
		
              	fprintf(fpOutfast, "\n");
		
              	if (name_flag)
        	    fprintf(fpOutfast, "+%s\n", seqp->name);
              	else
        	    fprintf(fpOutfast, "+\n");
		    
              	for (rc = 0; rc < seqp->length; rc++) {
        	    if (seqp->finished) {
                    	if ((seqp->data[rc] == 'N') || (seqp->data[rc] == 'X'))
                     	    putc(2 + 041, fpOutfast);
                    	else
                     	    putc(40 + 041, fpOutfast);
			    
        	    } else {
                    	putc(seqp->qual[rc] + 041, fpOutfast);
		    }
              	}
		
              	fprintf(fpOutfast, "\n");
            }
    	} else {
            printf("missed: %d %s\n", i, rdnames[i]);
	}
    }
    
    fclose(fpOutfast);

    free(seq->name);
    free(seq);
    free(ctg2wgs_index);
    free_cmatrix(rdnames, 0, n_reads + 1, 0, Max_N_NameBase);
    
        
    printf("Job finished for %d contigs!\n", n_contig);
    
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read index    */
/* =============================== */
char **Read_Index(char *filename, int *ctg2wgs_index, fasta *seq, int nRead, int nSeq) {
/* =============================== */

    int i, j, k, n_reads = nSeq + nRead;
    FILE *namef;
    fasta *seqp;
    int i_reads, num_rd_find, stopflag;
    int *readIndex;
    int mapindex = 0, n_dupreads, num_nodup;
    int *dup_reads;
    char **DBname, dupname[Max_N_NameBase];
    char **rdnames;

    if ((readIndex = (int *)calloc(n_reads, sizeof(int))) == NULL) {
    	fprintf(stderr, "Error: out of memory: calloc - readIndex\n");
      	exit(1);
    } 
    
    DBname  = cmatrix(0, n_reads,   0, Max_N_NameBase);
    rdnames = cmatrix(0, nRead + 1, 0, Max_N_NameBase);
    
    
    if (!DBname || !rdnames) {
    	fprintf(stderr, "Error: unable to make cmatrix, out of memory\n");
	exit(1);
    }
    
    if ((namef = fopen(filename, "r")) == NULL) {
      	fprintf(stderr, "ERROR: Unable to open file %s \n", filename);
      	exit(1);
    }

    i_reads = 0;
    
    while (fscanf(namef, "%s", rdnames[i_reads]) != EOF) {
    	i_reads++;
    }
    
    fclose(namef);

    for (j = 0; j < nSeq; j++) {
        seqp = seq + j;
        strcpy(DBname[j], seqp->name);
        readIndex[j] = j;
    }
    
    n_reads = nSeq;
    ArraySort_String(n_reads, DBname, readIndex);

    if ((dup_reads = (int *)calloc(nSeq, sizeof(int))) == NULL) {
    	fprintf(stderr, "ERROR Contig_Hist: calloc - dup_reads\n");
        exit(1);
    }
    
    n_dupreads = 0;
    
    for (i = 0; i <n_reads - 1; i++) {
    	/* search reads with an index < i     */
    	/* search reads with an index > i     */
        stopflag = 0;
        j = i + 1;
       
        while ((j < n_reads) && (stopflag == 0)) {
            if (strcmp(DBname[j], DBname[i]) == 0) 
            	j++;
            else
            	stopflag = 1;
        }
       
        if (( j - i) > 1) {
            for (k = (i + 1); k < j; k++) {
            	dup_reads[readIndex[k]] = 1;
            	n_dupreads++;
            }
        }
	
        i = j - 1;
    }
    
    printf("number of dup reads: %d %d\n", nSeq, n_dupreads);
    
    num_nodup = nSeq;
    
    for (j = 0; j < nSeq; j++) {
        if (dup_reads[j] != 1) {
            seqp = seq + j;
            strcpy(DBname[j], seqp->name);
	    
        } else {
            sprintf(dupname, "%s%09d", "dupname-", j);
            strcpy(DBname[j], dupname);
        }
	
        readIndex[j] = j;
    }
    
    free(dup_reads); dup_reads = NULL; 

    for (j = 0; j < nRead; j++) {
        strcpy(DBname[j + num_nodup], rdnames[j]);
        readIndex[j + num_nodup] = j + num_nodup;
    }
    
    n_reads = num_nodup + nRead;
    
    printf("before sort: %d %d\n", num_nodup, nRead);
    
    ArraySort_String(n_reads, DBname, readIndex);

    num_rd_find = 0;
    mapindex = 0;
    
    for (i = 0; i < n_reads - 1; i++) {
    	if (readIndex[i] >= num_nodup) mapindex = readIndex[i];
	 
    	/* search reads with an index < i     */
    	/* search reads with an index > i     */
        stopflag = 0;
        j = i + 1;
	
	while ((j < n_reads) && (stopflag == 0)) {
            if (strcmp(DBname[j], DBname[i]) == 0) {
            	if (readIndex[j] >= num_nodup) mapindex = readIndex[j];
		
            	j++;
		
            } else {
            	stopflag = 1;
	    }
	}
       
	if ((j - i) >= 2) {
            for (k = i; k < j; k++) {
            	if (readIndex[k] < num_nodup) {
                    ctg2wgs_index[mapindex - num_nodup] = readIndex[k];
                    num_rd_find++;
                    k = j;
             	}
            }
	}
	
	i = j - 1;
    }
    
    free(readIndex);
    free_cmatrix(DBname, 0, n_reads, 0, Max_N_NameBase);
    
    printf("reads found: %d %d\n", nRead, num_rd_find);
    
    return rdnames;
}

