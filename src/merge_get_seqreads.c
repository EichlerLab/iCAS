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
* merge_get_seqreads.c - Code modified for Illumina Clone Assembly Pipeline
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fasta.h"


int main(int argc, char **argv) {
    FILE *namef, *fpOutfast;
    long totalBases;
    int i, j, k, nSeq;
    int n_reads, n_base;
    fasta *seqp, *sub;
    char temp[100];
    char line[2000] = {0};
    int *hit_indexs, *hit_indexq, *hit_contigs, *hit_index2, *ctg_lengthq, *ctg_lengths;
    int *hit_locus1, *hit_locus2, *hit_rcdex1, *hit_rcdex2;
    long Size_q_pdata;
    fasta *seq;
    char *pdata, *qdata;

    fflush(stdout);

    if (argc < 2) {
      	printf("Usage: %s <merge_link_file> <input_ref_fastq_file> <input_merge_fastq_file> <output_fastq_file> \n",argv[0]);
      	exit(1);
    }

    if ((nSeq = extract_fastq_from_file(argv[2], &pdata, &Size_q_pdata)) == 0) {
    	fprintf(stderr, "Error: unable to extract fastq data from %s\n", argv[1]);
	exit(1);
    }

    if ((seq = decode_fastq(pdata, Size_q_pdata, &nSeq, &totalBases)) == NULL) {
    	fprintf(stderr, "Error: no query data found.\n");
	exit(1);
    }
    
    fastaUC(seq, nSeq);

    if ((nSeq = extract_fastq_from_file(argv[3], &qdata, &Size_q_pdata)) == 0) {
    	fprintf(stderr, "Error: unable to extract fastq data from %s\n", argv[1]);
	exit(1);
    }

    if ((sub = decode_fastq(qdata, Size_q_pdata, &nSeq, &totalBases)) == NULL) {
    	fprintf(stderr, "Error: no query data found.\n");
	exit(1);
    }

    fastaUC(sub, nSeq);

    if ((namef = fopen(argv[1],"r")) == NULL) {
    	fprintf(stderr, "ERROR main:: unable to open %s \n", argv[1]);
    	exit(1);
    }
    
    n_reads = 0;
    
    while (!feof(namef)) {
	fgets(line, 2000, namef);
	if (feof(namef)) break;
	n_reads++;
    }
    
    rewind(namef); 

    if ((hit_contigs = (int *)calloc(n_reads, sizeof(int))) == NULL) {
    	fprintf(stderr, "ERROR Contig_Hist: calloc - hit_contigs\n");
        exit(1);
    }
    
    if ((hit_index2 = (int *)calloc(n_reads, sizeof(int))) == NULL) {
    	fprintf(stderr, "ERROR Contig_Hist: calloc - hit_index2\n");
        exit(1);
    }
    
    if ((hit_indexs = (int *)calloc(n_reads, sizeof(int))) == NULL) {
        fprintf(stderr, "ERROR Contig_Hist: calloc - hit_indexs\n");
        exit(1);
    }
    
    if ((hit_indexq = (int *)calloc(n_reads, sizeof(int))) == NULL) {
        fprintf(stderr, "ERROR Contig_Hist: calloc - hit_indexq\n");
        exit(1);
    }
    
    if ((hit_rcdex1 = (int *)calloc(n_reads, sizeof(int))) == NULL) {
        fprintf(stderr, "ERROR Contig_Hist: calloc - hit_rcdex1\n");
        exit(1);
    }
    
    if ((hit_rcdex2 = (int *)calloc(n_reads, sizeof(int))) == NULL) {
        fprintf(stderr, "ERROR Contig_Hist: calloc - hit_rcdex2\n");
        exit(1);
    }
    
    if ((ctg_lengthq = (int *)calloc(n_reads, sizeof(int))) == NULL) {
        fprintf(stderr, "ERROR Contig_Hist: calloc - ctg_lengthq\n");
        exit(1);
    }
    
    if ((ctg_lengths = (int *)calloc(n_reads, sizeof(int))) == NULL) {
        fprintf(stderr, "ERROR Contig_Hist: calloc - ctg_lengths\n");
        exit(1);
    }
    
    if ((hit_locus1 = (int *)calloc(n_reads, sizeof(int))) == NULL) {
        fprintf(stderr, "ERROR Contig_Hist: calloc - hit_locus1\n");
        exit(1);
    }
    
    if ((hit_locus2 = (int *)calloc(n_reads, sizeof(int))) == NULL) {
        fprintf(stderr, "ERROR Contig_Hist: calloc - hit_locus2\n");
        exit(1);
    }

    i = 0;
    printf("reads: %d %d %s %s\n",n_reads, nSeq, argv[1], argv[2]);

    while (fscanf(namef, "%s %s %d %d %d %d %d %d %d %d", temp, temp, &hit_indexq[i], &hit_contigs[i], &hit_indexs[i], &hit_rcdex1[i], &hit_locus1[i], &hit_locus2[i], &ctg_lengthq[i], &ctg_lengths[i]) != EOF) {
    	i++;
    }
    
    fclose(namef);

    printf("reads2: %d %d %s %s\n", n_reads, nSeq, argv[1], argv[2]);
    
    if ((fpOutfast = fopen(argv[4], "w")) == NULL) {
      	fprintf(stderr, "ERROR main:: unable to open %s \n", argv[4]);
      	exit(1);
    }

    for (i = 0; i < n_reads; i++) {
        int stopflag;
        int rc;
        char *st;

        stopflag = 0;
        j = i + 1;
	
	while ((j < n_reads) && (stopflag == 0)) {
            if (hit_indexq[i] == hit_indexq[j]) {
            	j++;
            } else
            	stopflag = 1;
	}

        printf("contig_%06d\n", i);
        fprintf(fpOutfast, "@contig_%06d\n", i);
        n_base = 0;
	
	for (k = i; k < j; k++) {
            int rc, seq_st, seq_ed, r_len;

            if ((j - i) == 1) {
            	seqp  = (sub + hit_indexs[k]);
              	r_len = seqp->length;
	      	st    = seqp->data;
              	hit_rcdex1[k] = 0;
            } else {
            	if (hit_contigs[k] == 1) {
        	    seqp = (sub + hit_indexs[k]);
		    
        	    if ((k == i)||(k == (j - 1))) {
        	    	r_len = seqp->length;
		    	st = seqp->data;
        	    } else {
		    	r_len = abs(hit_locus2[k] - hit_locus1[k]);
		    	st = seqp->data + hit_locus1[k];
        	    }
              	} else {
        	    seqp = (seq+hit_indexq[k]);
        	    hit_rcdex1[k] = 0; 
		    r_len = abs(hit_locus2[k] - hit_locus1[k]);
		    st = seqp->data + hit_locus1[k];
              	}
            }
	    
	    if (hit_rcdex1[k] == 0) {
              	if (hit_contigs[k] == 1) {
		    if (k==i) {
        	    	if ((j - i) == 1) {
                    	    r_len = seqp->length;
                    	    st = seqp->data;
        	    	} else {
	            	    r_len = hit_locus2[k];
	            	    st = seqp->data;
        	    	}
		    } else if (k == (j - 1)) {
		    	r_len = (seqp->length) - hit_locus1[k];
		    	st = seqp->data + hit_locus1[k];
		    } else {
		    	r_len = abs(hit_locus2[k] - hit_locus1[k]);
		    	st = seqp->data + hit_locus1[k];
		    }
              	} else {
		    r_len = abs(hit_locus2[k] - hit_locus1[k]);
		    st = seqp->data + hit_locus1[k];
              	}
		
              	for (rc = 0; rc < r_len; rc++, st++) fprintf(fpOutfast, "%c", *st);
	    } else {
        	if (hit_contigs[k] == 1) {
		    if (k == i) {
        		seq_st = hit_locus2[k];
        		seq_ed = seqp->length;
        		r_len = seq_ed - seq_st;
		    } else if (k == (j - 1)) {
        		seq_st = 0;
        		seq_ed = hit_locus1[k];
        		r_len = hit_locus1[k];
		    } else {
        	    	seq_ed = hit_locus1[k];
        	    	seq_st = hit_locus2[k];
        	    	r_len = seq_ed - seq_st;
		    }
        	} else {
        	    seq_ed = hit_locus1[k];
        	    seq_st = hit_locus2[k];
        	    r_len = seq_ed - seq_st;
        	}
		
        	for (rc = (seq_ed - 1); rc >= seq_st; rc--) {
        	    if (seqp->data[rc] == 'A')
                     	fprintf(fpOutfast, "%c", 'T');
        	    else if (seqp->data[rc] == 'C')
                     	fprintf(fpOutfast, "%c", 'G');
        	    else if (seqp->data[rc] == 'G')
                     	fprintf(fpOutfast, "%c", 'C');
        	    else if (seqp->data[rc] == 'T')
                     	fprintf(fpOutfast, "%c", 'A');
        	    else
                     	fprintf(fpOutfast, "%c", 'N');
        	}
	    }
		
            n_base = n_base + r_len;
    	}

        fprintf(fpOutfast, "\n");
        fprintf(fpOutfast, "+\n");
	
        for (rc = 0; rc < n_base; rc++) putc(40 + 041, fpOutfast);
       
        fprintf(fpOutfast, "\n");
        i = j - 1; 
    }

    fclose(fpOutfast);
    printf("merge: %d\n", n_reads);
    
    return(0);
}
/* end of the main */
