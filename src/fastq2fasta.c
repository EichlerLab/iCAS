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
* fastq2fasta.c -  Code modified for Illumina Clone Assembly Pipeline                  *
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fasta.h"

int main(int argc, char **argv) {
    fasta *exp; 
    int ci;
    long Size_q_pdata, sBase;
    int num_seqque;
    char *pdata;
    FILE *outf, *outq;
    char qual[100];

    if (argc < 2) {
    	printf("Usage: %s <input_fastq> <output_fasta> \n",argv[0]);
      	exit(1);
    }
    
    if ((num_seqque = extract_fastq_from_file(argv[1], &pdata, &Size_q_pdata)) == 0) {
    	fprintf(stderr, "Error: unable to extract fastq data from %s\n", argv[1]);
	exit(1);
    }

    if ((exp = decode_fastq(pdata, Size_q_pdata, &num_seqque, &sBase)) == NULL) {
    	fprintf(stderr, "Error: no query data found.\n");
	exit(1);
    }

    fastaUC(exp, num_seqque);
    
    if ((outf = fopen(argv[2],"w")) == NULL) {
    	fprintf(stderr, "Error: can't open %s for writing.\n", argv[2]);
	exit(1);
    }
    
    strcpy(qual, argv[2]);
    strcat(qual, ".qual");
    
    if ((outq = fopen(qual,"w")) == NULL) {
    	fprintf(stderr, "Error: can't open %s for writing.\n", argv[2]);
	exit(1);
    }
    
    for (ci = 0; ci < num_seqque; ci++) {
    	int i;
    
	putc('>', outf);
	fprintf(outf, "%s", exp[ci].name);

	if (exp[ci].name2)
	    fprintf(outf, " %s", exp[ci].name2);

	putc('\n',outf);

	for (i = 0; i < exp[ci].length; i++) {
	    putc(exp[ci].data[i], outf);

	    if (i%60 == 59 && i + 1 != exp[ci].length)
	    	putc('\n',outf);
	}

	putc('\n', outf);
	putc('>', outq);
	fprintf(outq, "%s", exp[ci].name);

	if (exp[ci].name2)
	    fprintf(outq," %s", exp[ci].name2);

	putc('\n',outq);

	for (i = 0; i < exp[ci].length; i++) {
	    fprintf(outq, " %d", exp[ci].qual[i]);

	    if (i%30 == 29 && i + 1 != exp[ci].length)
	    	putc('\n', outq);
	}

	putc('\n', outq);
    }
    
    fclose(outf);
    fclose(outq);
    
    exit(0);
}
