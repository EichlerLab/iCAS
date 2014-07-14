/* 
* Copyright (c) 2011, 2014 Genome Research Ltd.
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
* convert_fastq.c - Code modified for Illumina Clone Assembly Pipeline
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fasta.h"


int main(int argc, char **argv) {
    FILE *fpOutfast;
    long totalBases;
    int i, nSeq;
    fasta *seq, *seqp;
    long Size_q_pdata;
    char *pdata;

    fflush(stdout);
    
    if(argc < 2) {
	printf("Usage: %s <-qual 15> <-base 3> <input fastq file> <output fastq file>>\n", argv[0]);
	exit(1);
    }

    if ((nSeq = extract_fastq_from_file(argv[1], &pdata, &Size_q_pdata)) == 0) {
    	fprintf(stderr, "Error: unable to extract fastq data from %s\n", argv[1]);
	exit(1);
    }

    if ((seq = decode_fastq(pdata, Size_q_pdata, &nSeq, &totalBases)) == NULL) {
    	fprintf(stderr, "Error: no query data found.\n");
	exit(1);
    }

    printf("Number of shotgun reads  %d \n", nSeq);

    /*  process contigs one by one   */
    if ((fpOutfast = fopen(argv[2],"w")) == NULL) {
	printf("ERROR main:: reads group file %s unwritable\n", argv[2]);
	exit(1);
    }
    
    for (i = 0; i < nSeq; i++) {
    	int seq_st, seq_ed, rc;

	seqp = seq + i;
	seq_st = 0;
	seq_ed = seqp->length;
       
	fprintf(fpOutfast, "@%s\n", seqp->name);

	for (rc = seq_st; rc < seq_ed; rc++) {
	    fprintf(fpOutfast, "%c", seqp->data[rc]);
	}

	fprintf(fpOutfast, "\n");
	fprintf(fpOutfast, "+%s\n", seqp->name);
	putc(041, fpOutfast);

	for (rc = (seq_st+1); rc < seq_ed; rc++) {
	    if (seqp->finished) {
	    	int idd = (int)(drand48() * 10);

	    	if ((seqp->data[rc] == 'N') || (seqp->data[rc] == 'X'))
	    	    putc(2 + 041, fpOutfast);
	    	else
	    	    putc(40 + 041 + idd, fpOutfast);
	    } else {
	    	putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	}

	fprintf(fpOutfast, "\n");
    }
    
    fclose(fpOutfast);

    if (seq) {
        free(seq->name);
        free(seq);
        seq = NULL;
    }
        
    printf("Job finished for %d contigs!\n", nSeq);
    return EXIT_SUCCESS;

}

/* end of the main */


