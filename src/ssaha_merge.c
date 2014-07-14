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
* ssaha_merge.c - Code modified for Illumina Clone Assembly Pipeline
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


int main(int argc, char **argv)
{
    FILE *fpOutfast;
    int i, nSeq1, nSeq2, args, out1_len, out2_len;
    char line[7] = {0},rdname[100];
    fasta *seq1, *seq2; 
    fasta *seqp, *seqg;
    B64_long Size_q_pdata;
    char *pdata1, *pdata2;
    B64_long sBase1, sBase2;

    /* SSAS default parameters   */
    int set_len=120;
    int trim_len=1200;
    int file_flag = 0;

    if (argc < 2) {
	printf("Usage: %s -file x <input_fastq_file_1> <input_fastq_file_2> <ouput_fastq_file>\n", argv[0]);
	printf("-file 0: reads in two files are rearranged into one file with .F and .R\n");
	printf("-file 1: reads in two files with different length are rearranged into one file with read length: R1+R2+R1-R2, where R1-R2 = XXX.\n");
	printf("-file 2: reads in two files are rearranged into one file with .F and .R - insert with NN\n");
	printf("-file 3: reads in two files are rearranged into one file with 2*length\n");
	exit(1);
    }

    args = 1;

    for (i = 1; i < argc; i++) {
    	if (!strcmp(argv[i], "-trim")) {
            sscanf(argv[++i], "%d", &trim_len);
            args += 2;
	    
        } else if (!strcmp(argv[i], "-len")) {
            sscanf(argv[++i], "%d", &set_len);
            args +=2;
	    
        } else if (!strcmp(argv[i], "-file")) {
            sscanf(argv[++i], "%d", &file_flag);
            args +=2;
        }
    }

    if ((nSeq1 = extract_fastq_from_file(argv[args], &pdata1, &Size_q_pdata)) == 0) {
    	fprintf(stderr, "Error: unable to extract fastq data from %s\n", argv[args]);
	exit(1);
    }

    if ((seq1 = decode_fastq(pdata1, Size_q_pdata, &nSeq1, &sBase1)) == NULL) {
    	fprintf(stderr, "Error: no query data found.\n");
	exit(1);
    }

    fastaUC(seq1, nSeq1);

    if ((nSeq2 = extract_fastq_from_file(argv[args + 1], &pdata2, &Size_q_pdata)) == 0) {
    	fprintf(stderr, "Error: unable to extract fastq data from %s\n", argv[args + 1]);
	exit(1);
    }

    if ((seq2 = decode_fastq(pdata2, Size_q_pdata, &nSeq2, &sBase2)) == NULL) {
    	fprintf(stderr, "Error: no query data found.\n");
	exit(1);
    }

    fastaUC(seq2, nSeq2);

    /*  output fastq file   */
    if ((fpOutfast = fopen(argv[args + 2], "w")) == NULL) {
      	fprintf(stderr, "ERROR main:: reads group file: %s \n", argv[args + 2]);
      	exit(1);
    }

    if (nSeq1 != nSeq2) {
      	fprintf(stderr, "Error: read number in %s and %s not the same: %d %d\n", argv[args], argv[args + 1], nSeq1, nSeq2);
      	fprintf(stderr, "program stopped!\n");
      	exit(1);
    }
    
    for (i = 0; i < nSeq1; i++) {
    	int rc, seq1_len, seq2_len, name_len;

        seqp = seq1 + i;
        seqg = seq2 + i;
        seq1_len = seqp->length;
        seq2_len = seqg->length;
       
	if (file_flag == 0) {
            name_len = strlen(seqp->name);
	    memset(rdname, '\0', 100);
	    strncpy(rdname, seqp->name, name_len - 2);
	    fprintf(fpOutfast, "@%s.F\n", rdname);
	    
  	    for (rc = 0; rc < seq1_len; rc++) {
	    	fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    putc(0 + 041, fpOutfast);
	    
	    for (rc = 1; rc < seq1_len; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "@%s.R\n", rdname);
	    
  	    for (rc = 0;rc < seq2_len; rc++) {
	        fprintf(fpOutfast, "%c", seqg->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    putc(0 + 041, fpOutfast);
	    
	    for (rc = 1; rc < seq2_len; rc++) {
	        putc(seqg->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    
	} else if (file_flag == 1) {
            if (seq1_len >= seq2_len) {
		fprintf(fpOutfast, "@%s\n", seqp->name);
		
  		for (rc = 0; rc < seq1_len; rc++) {
		    fprintf(fpOutfast, "%c", seqp->data[rc]);
		}
		
  		for (rc = 0; rc < seq2_len; rc++) {
		    fprintf(fpOutfast, "%c", seqg->data[rc]);
		}
		
		for (rc = 0; rc < (seq1_len - seq2_len); rc++) {
		    fprintf(fpOutfast, "X");
		}
		
		fprintf(fpOutfast, "\n");
		fprintf(fpOutfast, "+\n");
		
		for (rc = 0; rc < seq1_len; rc++) {
		    putc(seqp->qual[rc] + 041, fpOutfast);
		}
		
		for (rc = 0; rc < seq2_len; rc++) {
		    putc(seqg->qual[rc] + 041, fpOutfast);
		}
		
		for (rc = 0; rc < (seq1_len - seq2_len); rc++) {
		    putc(0 + 041, fpOutfast);
		}
		
		fprintf(fpOutfast, "\n");
		
            } else {
		fprintf(fpOutfast, "@%s\n", seqp->name);
		
		for (rc = 0; rc < seq1_len; rc++) {
		    fprintf(fpOutfast, "%c", seqp->data[rc]);
		}
		
		for (rc = 0; rc < (seq2_len - seq1_len); rc++) {
		    fprintf(fpOutfast ,"X");
		}
		
		for (rc = 0; rc < seq2_len; rc++) {
		    fprintf(fpOutfast, "%c", seqg->data[rc]);
		}
		
		fprintf(fpOutfast, "\n");
		fprintf(fpOutfast, "+\n");
		
		for (rc = 0; rc < seq1_len; rc++) {
		    putc(seqp->qual[rc] + 041, fpOutfast);
		}
		
		for (rc = 0; rc < (seq2_len - seq1_len); rc++) {
		    putc(0 + 041, fpOutfast);
		}
		
		for (rc = 0; rc < seq2_len; rc++) {
		    putc(seqg->qual[rc] + 041, fpOutfast);
		}
		
		fprintf(fpOutfast, "\n");
            }
	    
	} else if (file_flag == 2) {
            name_len = strlen(seqp->name);
	    memset(rdname, '\0', 100);
	    strncpy(rdname, seqp->name, name_len - 2);
	    fprintf(fpOutfast, "@%s\n", rdname);
	    
	    if (seq1_len < trim_len)
	      	out1_len = seq1_len;
	    else
	      	out1_len = trim_len;
	      
	    if (seq2_len<trim_len)
	      	out2_len = seq2_len;
	    else
	      	out2_len = trim_len;
		
  	    for (rc = 0; rc < out1_len; rc++) {
	        fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    for (rc = 0; rc < 2; rc++) {
	        fprintf(fpOutfast, "N");
	    }
	    
  	    for (rc = 0; rc < out2_len; rc++) {
	        fprintf(fpOutfast, "%c", seqg->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for (rc = 0; rc < out1_len; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    for (rc = 0; rc < 2; rc++) {
	    	putc(0 + 041, fpOutfast);
	    }
	    
	    for (rc = 0; rc < out2_len; rc++) {
	        putc(seqg->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    
	} else if(file_flag == 20) {
            name_len = strlen(seqp->name);
	    memset(rdname, '\0', 100);
	    strncpy(rdname, seqp->name, name_len - 2);
	    fprintf(fpOutfast, "@%s\n", rdname);

  	    for (rc = 0; rc < set_len; rc++) {
	        fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
	    for (rc = 0; rc < 2; rc++) {
	        fprintf(fpOutfast, "N");
	    }
	    
  	    for (rc = 0; rc < set_len; rc++) {
	        fprintf(fpOutfast, "%c", seqg->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for (rc = 0; rc < set_len; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    for (rc = 0; rc < 2; rc++) {
	        putc(0 + 041, fpOutfast);
	    }
	    
	    for (rc = 0; rc < set_len; rc++) {
	        putc(seqg->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    
	} else if (file_flag == 22) {
            if (seq1_len >= seq2_len) {
		char *ed;
		ed = strchr(seqp->name, ':');
		
		if (i == 0) {
		    for (rc = 4; rc < 10; rc++) {
		     	line[rc - 4] = seqp->name[rc];
		    }
		}
		
		fprintf(fpOutfast, "@%s%s\n", line,ed);
		
  		for (rc = 0; rc < seq1_len; rc++) {
		    fprintf(fpOutfast, "%c", seqp->data[rc]);
		}
		
		for (rc = 0; rc < 2; rc++) {
		    fprintf(fpOutfast, "N");
		}
		
  		for (rc = 0; rc < seq2_len; rc++) {
		    fprintf(fpOutfast, "%c", seqg->data[rc]);
		}
		
		fprintf(fpOutfast, "\n");
		fprintf(fpOutfast, "+\n");
		putc(0 + 041, fpOutfast);
		
		for (rc = 1; rc < seq1_len; rc++) {
		    putc(seqp->qual[rc] + 041, fpOutfast);
		}
		
		for (rc = 0; rc < 2; rc++) {
		    putc(0 + 041, fpOutfast);
		}
		
		for (rc = 0; rc < seq2_len; rc++) {
		    putc(seqg->qual[rc] + 041, fpOutfast);
		}
		
		fprintf(fpOutfast, "\n");
            }
	    
	} else if (file_flag == 3) {
	    fprintf(fpOutfast, "@%s\n", seqp->name);
	    
  	    for (rc = 0; rc < seq1_len; rc++) {
	        fprintf(fpOutfast, "%c", seqp->data[rc]);
	    }
	    
  	    for (rc = 0; rc < seq2_len; rc++) {
	        fprintf(fpOutfast, "%c", seqg->data[rc]);
	    }
	    
	    fprintf(fpOutfast, "\n");
	    fprintf(fpOutfast, "+\n");
	    
	    for (rc = 0; rc < seq1_len; rc++) {
	        putc(seqp->qual[rc] + 041, fpOutfast);
	    }
	    
	    for (rc = 0; rc < seq2_len; rc++) {
	        putc(seqg->qual[rc] + 041, fpOutfast);
	    }
	    
	    fprintf(fpOutfast, "\n");
	}
    }
    
    fclose(fpOutfast);
    free(pdata1);
    free(seq1);
    free(pdata2);
    free(seq2);
    
    return EXIT_SUCCESS;

}
/* end of the main */
