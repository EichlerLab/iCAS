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
* ssaha_reads.c - Code modified for Illumina Clone Assembly Pipeline
*/

#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <sys/wait.h>
#include <sys/signal.h>
#include <errno.h>
#include "fasta.h"

int main(int argc, char **argv) {
    FILE *fpOutfast, *fpOutfast2 = NULL;
    int i, nSeq, args, seq_len;
    char rdname[100], rdname2[100];
    int name_len = 0, name_len2 = 0, seq_st = 0, seq_ed = 0;
    fasta *seq;
    fasta *seqp;
    B64_long Size_q_pdata;
    char *pdata;

    B64_long sBase;

    /* SSAS default parameters   */
    int num_reads = 0;
    int seq_len1 = 76;
    int file_flag = 2;

    if (argc < 2) {
    	// add the other usable options
    	printf("Usage: %s [-file number -len number] <input_fastq_reads_file> <ouput_fastq_file>\n", argv[0]);
      	exit(1);
    }

    nSeq = 0;
    args = 1;
    
    for (i = 1; i < argc; i++) {
    	if (strcmp(argv[i], "-len") == 0) {
            sscanf(argv[++i], "%d", &seq_len1);
            args += 2;
        } else if(strcmp(argv[i], "-file") == 0) {
            sscanf(argv[++i], "%d", &file_flag);
            args += 2;
        }
    }

    if ((nSeq = extract_fastq_from_file(argv[args], &pdata, &Size_q_pdata)) == 0) {
    	fprintf(stderr, "Error: unable to extract fastq data from %s\n", argv[args]);
	exit(1);
    }

    if ((seq = decode_fastq(pdata, Size_q_pdata, &nSeq, &sBase)) == NULL) {
    	fprintf(stderr, "Error: no query data found.\n");
	exit(1);
    }
    
    fastaUC(seq, nSeq);

    /*  output fastq file   */
    if ((fpOutfast = fopen(argv[args + 1], "w")) == NULL) {
    	printf("ERROR main:: reads group file: %s \n", argv[args + 1]);
      	exit(1);
    }
    
    /*  output fastq2   */
    if ((file_flag == 2) || (file_flag == 3) || (file_flag == 22) || (file_flag == 23) ||
    	(file_flag == 33)|| (file_flag == 66) || (file_flag == 8) || (file_flag == 67) ||
	(file_flag == 6)) {
	
    	if ((fpOutfast2 = fopen(argv[args + 2],"w")) == NULL) {
            printf("ERROR main:: reads group file: %s \n", argv[args + 2]);
            exit(1);
      	}
    }

    if (file_flag == 23) {
    	int j;
	int *read_mask;
	
	if ((read_mask = (int *)calloc(10000, sizeof(int))) == NULL) {
    	    printf("ERROR Memory_Allocate: Align_Process - read_mask\n");
      	    exit(1);
	}
	
    	if (nSeq < 5000)
            num_reads = nSeq;
      	else
            num_reads = 5000;
	
      	for (i = 0; i < num_reads; i++) {
	    seqp = seq + i;
	    seq_len = seqp->length;

            for (j = 0; j < seq_len; j++) {
            	if (seqp->data[j] == 'N')
        	    read_mask[j]++;
            }

            for (j = 0; j < seq_len; j++) {
                if (read_mask[j] == num_reads) {
        	    seq_len1 = j;
        	    j = seq_len;
                }
            }
    	}
	
	free(read_mask);
    }
    
    for (i = 0; i < nSeq; i++) {
	int rc, n_half;

	seqp = seq + i;
	seq_len = seqp->length;

	if (file_flag == 0) {
            seq_st = 0;
	    seq_ed = seqp->length;
	    fprintf(fpOutfast, "@%s\n", seqp->name);

	    for (rc = seq_st; rc < seq_ed; rc++) {
	        if(seqp->data[rc] == '.') {
		    fprintf(fpOutfast, "N");
	        } else {
		    fprintf(fpOutfast, "%c", seqp->data[rc]);
	    	}
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    
	} else if (file_flag == 2) {
            seq_st = 0;
	    seq_ed = seqp->length / 2;
	    
	    fprintf(fpOutfast, "@%s/1\n", seqp->name);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");

            seq_st = seqp->length / 2;
            seq_ed = seqp->length;
	    
            fprintf(fpOutfast2, "@%s/2\n", seqp->name);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        fprintf(fpOutfast2, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast2, "\n");
	    fprintf(fpOutfast2, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast2);
	    }
	    
	    fprintf(fpOutfast2, "\n");
	    
	} else if (file_flag == 22) {
            seq_st = 0;
	    seq_ed = (seqp->length / 2) - 1;
            name_len = strlen(seqp->name);
            memset(rdname, '\0', 100);
            strncpy(rdname, seqp->name, name_len);
	    fprintf(fpOutfast, "@%s/1\n", rdname);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");

            seq_st = (seqp->length / 2) + 1;
            seq_ed = seqp->length;
	    fprintf(fpOutfast2, "@%s/2\n", rdname);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	fprintf(fpOutfast2, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast2, "\n");
	    fprintf(fpOutfast2, "+\n");
	    
	    for (rc= seq_st; rc < seq_ed; rc++) {
	    	putc(seqp->qual[rc] + 041, fpOutfast2);
	    }
	    
	    fprintf(fpOutfast2, "\n");
	    
	} else if (file_flag == 23) {
            seq_st = 0;
	    seq_ed = seq_len1;
            name_len = strlen(seqp->name);
            memset(rdname, '\0', 100);
            strncpy(rdname, seqp->name, name_len);
	    fprintf(fpOutfast, "@%s/1\n", rdname);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");

            seq_st = seq_len1 + 2;
            seq_ed = seqp->length;
	    fprintf(fpOutfast2,"@%s/2\n", rdname);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	fprintf(fpOutfast2, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast2, "\n");
	    fprintf(fpOutfast2, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	putc(seqp->qual[rc] + 041, fpOutfast2);
	    }
	    
	    fprintf(fpOutfast2,"\n");
	    
	} else if (file_flag == 3) {
            n_half = nSeq / 2;
            seq_st = 0;
	    seq_ed = seqp->length;

	    if ((i%2) == 0) {
		fprintf(fpOutfast, "@%s\n", seqp->name);

		for (rc = seq_st; rc < seq_ed; rc++) {
		    fprintf(fpOutfast, "%c", seqp->data[rc]);
		}
		
		fprintf(fpOutfast, "\n");
		fprintf(fpOutfast, "+\n");
		
		for (rc = seq_st; rc < seq_ed; rc++) {
		    putc(seqp->qual[rc] + 041, fpOutfast);
		}
		
		fprintf(fpOutfast, "\n");
		
            } else {
        	fprintf(fpOutfast2, "@%s\n", seqp->name);
		
		for (rc = seq_st; rc < seq_ed; rc++) {
		    fprintf(fpOutfast2, "%c", seqp->data[rc]);
		}
		
		fprintf(fpOutfast2, "\n");
		fprintf(fpOutfast2, "+\n");
		
		for (rc = seq_st; rc < seq_ed; rc++) {
		    putc(seqp->qual[rc] + 041, fpOutfast2);
		}
		
		fprintf(fpOutfast2, "\n");
	    }
	    
	} else if (file_flag == 33) {
            seq_st = 0;
	    seq_ed = (seqp->length / 2);
            name_len = strlen(seqp->name);
            memset(rdname, '\0', 100);
            strncpy(rdname, seqp->name, name_len);
	    fprintf(fpOutfast, "@%s/1\n", rdname);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for(rc = seq_st; rc < seq_ed; rc++) {
	       putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");

            seq_st = (seqp->length / 2) + 1;
            seq_ed = seqp->length;
	    fprintf(fpOutfast2, "@%s/2\n", rdname);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	       fprintf(fpOutfast2, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast2, "\n");
	    fprintf(fpOutfast2, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	putc(seqp->qual[rc] + 041, fpOutfast2);
	    }
	    
	    fprintf(fpOutfast2, "\n");
	    
	} else if (file_flag == 4) {
            seq_st = 0;
	    seq_ed = seqp->length;
	    fprintf(fpOutfast, "@%s\n", seqp->name);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
            seqp = seq + i + 1;
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
            seqp = seq + i;
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
            seqp = seq + i + 1;
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");
            i++;
	    
	} else if (file_flag == 44) {
            seq_st = 0;
	    seq_ed = seqp->length;
            name_len = strlen(seqp->name);
            memset(rdname, '\0', 100);
            strncpy(rdname, seqp->name, name_len - 2);
	    fprintf(fpOutfast, "@%s\n", rdname);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    for (rc = 0;rc < 2; rc++) {
	    	fprintf(fpOutfast, "N");
	    }
	    
            seqp = seq + i + 1;
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
            seqp= seq + i;
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    for (rc = 0; rc < 2; rc++) {
	        putc(0 + 041, fpOutfast);
	    }
	    
            seqp= seq + i + 1;
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    
            i++;
	    
	} else if (file_flag == 5) {
            seq_st = 0;
	    seq_ed = (seqp->length / 2) - 1;
	    fprintf(fpOutfast, "@%s.F\n", seqp->name);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");
            seq_st = (seqp->length / 2) + 1;;
	    seq_ed = seqp->length;
	    fprintf(fpOutfast, "@%s.R\n", seqp->name);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    
	} else if (file_flag == 6) {
            seq_st = 0;
	    seq_ed = (seqp->length / 2) - 1;
	    fprintf(fpOutfast, "@%s/1\n", seqp->name);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast,"\n");

            seq_st = (seqp->length / 2) + 1;
            seq_ed = seqp->length;
	    fprintf(fpOutfast2, "@%s/2\n", seqp->name);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        fprintf(fpOutfast2, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast2, "\n");
	    fprintf(fpOutfast2, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast2);
	    }
	    
	    fprintf(fpOutfast2, "\n");
	    
	} else if (file_flag == 66) {
            seq_st = 0;
            seq_ed = seqp->length;
            name_len = strlen(seqp->name);
            memset(rdname, '\0', 100);
            strncpy(rdname,seqp->name, name_len - 2);
	    fprintf(fpOutfast, "@%s/1\n", rdname);
	    
            for (rc = seq_st; rc < seq_ed; rc++) {
                fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
            fprintf(fpOutfast, "\n");
            fprintf(fpOutfast, "+\n");
	    
            for (rc = seq_st; rc < seq_ed; rc++) {
                putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
            fprintf(fpOutfast, "\n");

            i++;
	    seqp = seq + i;
	    seq_st = 0;
	    seq_ed = seqp->length;
            name_len = strlen(seqp->name);
            memset(rdname, '\0', 100);
            strncpy(rdname, seqp->name, name_len - 2);
	    fprintf(fpOutfast2, "@%s/2\n", rdname);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        fprintf(fpOutfast2, "%c", seqp->data[rc]);
	    }
	    
            fprintf(fpOutfast2, "\n");
	    fprintf(fpOutfast2, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast2);
	    }
	    
	    fprintf(fpOutfast2, "\n");
	    
	} else if (file_flag == 67) {
            seq_st = 0;
            seq_ed = seqp->length;
            name_len = strlen(seqp->name);
            memset(rdname, '\0', 100);
            strncpy(rdname, seqp->name, name_len);
	    fprintf(fpOutfast, "@%s/1\n", rdname);
	    
            for (rc = seq_st; rc < seq_ed; rc++) {
                fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
            fprintf(fpOutfast, "\n");
            fprintf(fpOutfast, "+\n");
	    
            for (rc = seq_st; rc < seq_ed; rc++) {
                putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
            fprintf(fpOutfast, "\n");

            i++;
	    seqp = seq + i;
	    seq_st = 0;
	    seq_ed = seqp->length;
	    fprintf(fpOutfast2, "@%s/2\n", rdname);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        fprintf(fpOutfast2, "%c", seqp->data[rc]);
	    }
	    
            fprintf(fpOutfast2, "\n");
	    fprintf(fpOutfast2, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast2);
	    }
	    
	    fprintf(fpOutfast2, "\n");
	    
	} else if (file_flag == 7) {
	    fprintf(fpOutfast, "@SQ\tSN:%s\tLN:%d\n", seqp->name, seqp->length);
	    
	} else if (file_flag == 8) {
	    if (i <(nSeq - 1)) {
            	seq_st = 0;
            	seq_ed = seqp->length;
        	name_len = strlen(seqp->name);
        	name_len2 = strlen((seqp + 1)->name);
        	memset(rdname, '\0', 100);
        	memset(rdname2, '\0', 100);
        	strncpy(rdname, seqp->name, name_len - 2);
        	strncpy(rdname2, (seqp + 1)->name, name_len - 2);
	   }
	 
           if (name_len && (strcmp(rdname, rdname2) == 0) && (seqp->name[name_len - 1] == '1') && (i < (nSeq - 1))) {
                fprintf(fpOutfast, "@%s/1\n", rdname);
	       
                for (rc = seq_st; rc < seq_ed; rc++) {
        	    fprintf(fpOutfast,"%c",seqp->data[rc]);
		}
		  
        	fprintf(fpOutfast, "\n");
        	fprintf(fpOutfast, "+\n");
		
        	for (rc = seq_st; rc < seq_ed; rc++) {
        	    putc(seqp->qual[rc] + 041, fpOutfast);
		}
		
        	fprintf(fpOutfast, "\n");

        	i++;
		seqp = seq + i;
		seq_st = 0;
		seq_ed = seqp->length;
        	name_len = strlen(seqp->name);
        	memset(rdname, '\0', 100);
        	strncpy(rdname, seqp->name, name_len - 2);
		fprintf(fpOutfast2, "@%s/2\n", rdname);
		
		for (rc = seq_st; rc < seq_ed; rc++) {
		    fprintf(fpOutfast2, "%c", seqp->data[rc]);
		}
		
        	fprintf(fpOutfast2, "\n");
		fprintf(fpOutfast2, "+\n");
		
		for (rc = seq_st; rc < seq_ed; rc++) {
		    putc(seqp->qual[rc] + 041, fpOutfast2);
		}
		
		fprintf(fpOutfast2, "\n");
            }
	    
	} else {
            seq_st = 0;
	    seq_ed = seqp->length / 2;
	    
	    fprintf(fpOutfast, "@%s.p1k\n", seqp->name);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");

            seq_st = seqp->length / 2;
            seq_ed = seqp->length;
            fprintf(fpOutfast, "@%s.q1k\n", seqp->name);
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	    	fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for (rc = seq_st; rc < seq_ed; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");
	}
    }
    
    if (file_flag == 7) fprintf(fpOutfast, "@SQ\tSN:*\tLN:0\n");
    
    fclose(fpOutfast);
    
    if (fpOutfast2) {
    	fclose(fpOutfast2);
    }
    
    free(pdata);
    free(seq);
    
    return EXIT_SUCCESS;
}
