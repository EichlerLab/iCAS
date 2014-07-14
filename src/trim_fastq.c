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
*    trim_fastq.c - trim the length of the reads in a fastq file
*/    

#include <stdio.h>
#include <stdlib.h>

#include "fasta.h"



int main (int argc, char **argv) {
    fasta *seq, *seqp;
    char *pdata;
    char *infile, *outfile;
    int num_seqs;
    B64_long total_bases, data_size;
    FILE *fp;
    int i;
    char *endptr;
    long trim_len;
    
    if (argc != 4) {
    	fprintf(stderr, "Usage: %s trim_length infile outfile\n", argv[0]);
	exit(1);
    }
    
    trim_len = strtol(argv[1], &endptr, 10);    
    
    if (*endptr || trim_len < 1) {
    	fprintf(stderr, "Error: %s is not a valid number.\n", argv[1]);
	exit(1);
    } 
    
    infile  = argv[2];
    outfile = argv[3];
    
    num_seqs = extract_fastq_from_file(infile, &pdata, &data_size);
    
    if (num_seqs == 0) {
    	fprintf(stderr, "Error, no reads found in %s.\n", infile);
	exit(1);
    }
    
    if ((seq = decode_fastq(pdata, data_size, &num_seqs, &total_bases)) == NULL) {
    	fprintf(stderr, "Error: no query data found.\n");
	exit(1);
    }
    
    fastaUC(seq, num_seqs);
    
    if ((fp = fopen(outfile, "w")) == NULL) {
      	printf("ERROR main:: can't write to %s\n", outfile);
      	exit(1);
    }
    
    for (i = 0; i < num_seqs; i++) {
	int rc;
	int seq_st, seq_ed;
	
	seqp   = seq + i;

	if (trim_len > (seqp->length / 2)) {
	    trim_len = ((seqp->length - 2) / 2);
	}
	
	seq_st = 0;
	seq_ed = trim_len;

	fprintf(fp, "@%s\n", seqp->name);

	for (rc = seq_st; rc < seq_ed; rc++)
	    fprintf(fp, "%c", seqp->data[rc]);

	for (rc = 0;rc < 2; rc++)
	    fprintf(fp, "%c", 'N');

	seq_st = (seqp->length / 2) + 1;
	seq_ed = (seqp->length / 2) + trim_len + 1;

	for (rc = seq_st; rc < seq_ed; rc++)
	    fprintf(fp, "%c", seqp->data[rc]);

	fprintf(fp, "\n");
	fprintf(fp, "+\n");
	seq_st = 0;
	seq_ed = trim_len;

	for (rc = seq_st; rc < seq_ed; rc++)
	    putc(seqp->qual[rc] + 041, fp);

	for (rc = 0; rc < 2; rc++)
	    putc(0 + 041, fp);

	seq_st = (seqp->length / 2) + 1;
	seq_ed = (seqp->length / 2) + trim_len + 1;

	for (rc = seq_st; rc < seq_ed; rc++)
	    putc(seqp->qual[rc] + 041, fp);

	fprintf(fp, "\n");
    }
    
    fclose(fp);
    
    exit(0);
}
    
 
