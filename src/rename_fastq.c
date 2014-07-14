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
* rename_fastq.c - Code modified for Illumina Clone Assembly Pipeline
*/
 
#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"

int main(int argc, char **argv) {
    FILE *fpOutfast;
    long totalBases;
    int i,nSeq,args;
    fasta *seq,*seqp;
    char ctgname[30];
    long Size_q_pdata;
    char *pdata;
    int set_qual = 15;
    int set_len = 10;

    fflush(stdout);
    strcpy(ctgname, "Contig");

    if(argc < 2) {
      printf("Usage: %s <-name contigs> <-len 10> <input fastq file> <output fastq file>>\n", argv[0]);
      exit(1);
    }

    args = 1;

    for (i = 1; i < argc; i++) {
	if (!strcmp(argv[i], "-qual")) {
            sscanf(argv[++i], "%d", &set_qual);
            args = args + 2;
	} else if (!strcmp(argv[i], "-len")) {
            sscanf(argv[++i], "%d", &set_len);
            args = args + 2;
	} else if (!strcmp(argv[i], "-name")) {
            memset(ctgname, '\0', 30);
            sscanf(argv[++i], "%s", ctgname);
            args = args + 2;
	}
    }

    printf("name:  %s\n", ctgname);
    

    if ((nSeq = extract_fastq_from_file(argv[args], &pdata, &Size_q_pdata)) == 0) {
    	fprintf(stderr, "Error: unable to extract fastq data from %s\n", argv[1]);
	exit(1);
    }

    if ((seq = decode_fastq(pdata, Size_q_pdata, &nSeq, &totalBases)) == NULL) {
    	fprintf(stderr, "Error: no query data found.\n");
	exit(1);
    }
    
    printf("Number of shotgun reads  %d \n",nSeq);

    /*  process contigs one by one   */
    if ((fpOutfast = fopen(argv[args+1], "w")) == NULL) {
    	fprintf(stderr, "ERROR main:: unable to write to %s\n", argv[args + 1]);
      	exit(1);
    }
    
    for (i = 0; i < nSeq; i++) {
        int seq_st, seq_ed, rc;

        seqp   = seq + i;
        seq_st = 0;
        seq_ed = seqp->length;
       
        if (seqp->length>=set_len) {
            fprintf(fpOutfast, "@%s_%09d\n", ctgname, i);
	    
            for(rc = seq_st; rc < seq_ed; rc++) fprintf(fpOutfast, "%c", seqp->data[rc]);
	    
            fprintf(fpOutfast, "\n");
            fprintf(fpOutfast, "+%s_%09d\n", ctgname, i);
            putc(041, fpOutfast);
	    
            for (rc = (seq_st + 1); rc < seq_ed; rc++) {
                if (seqp->finished) {
		    int idd = (int)(drand48() * 10);
		    
        	    if ((seqp->data[rc] == 'N') || (seqp->data[rc] == 'X'))
                    	putc(2 + 041, fpOutfast);
        	    else
                    	putc(40 + 041 + idd, fpOutfast);
                } else
            	    putc(seqp->qual[rc] + 041, fpOutfast);
            }
	    
            fprintf(fpOutfast, "\n");
	}
    }
    
    fclose(fpOutfast);

    if (seq){
        free(seq->name);
        free(seq);
        seq = NULL;
    }
        
    printf("Job finished for %d contigs!\n",nSeq);
    return EXIT_SUCCESS;
}
/* end of the main */
