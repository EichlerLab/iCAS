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
* phusion2_main.c - Code modified for Illumina Clone Assembly Pipeline
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
#include "ssaha.h"

static int KLEN=17;

int main(int argc, char **argv)
{
	int i;
    int SEG_LEN=KLEN;
//	__delayed_free = 0;//so we can free memory right away

	printf("SsahaHist Version 1.02\n");
	if(argc < 2)
	{
		printf("Usage: %s command [options] \n",argv[0]);
		printf("\n");
		printf("Commands:\n");
		printf("          hash     generate hash files from sequencing reads\n");
		printf("          Options: phusion2 hash <-kmer 31>\n");
		printf("\n");
		printf("          sort     sort hash files files\n");
		printf("          Options: phusion2 sort <-kmer 31> reads_0*.fastq (all the read files)\n");
		printf("\n");
		printf("          edge     generate relation matrix file (sparse)\n");
		printf("          Options: phusion2 edge <-kmer 31> <-depth 20> <-edge 50> reads_0*.fastq\n");
		printf("\n");
		printf("          matrix   generate relation matrix file (compact)\n");
		printf("          Options: phusion2 matrix <-kmer 31> <-edge 50> reads_0*.fastq\n");
		printf("\n");
		printf("          clust    cluster reads into groups:\n");
		printf("          Options: phusion2 clust <-kmer 31> <-depth 20> <-match 4> <-set 60000> > clust-d20.out \n");
		printf("\n");
		printf("          remap    remap the reads which were not used initially\n");
		printf("          Options: phusion2 remap <-kmer 31> <-depth 20> <-match 2> <-set 60> clust-d20.dat clust_remap.dat \n");
		printf("\n");
		printf("          merge    merge clusters using remapped reads:\n");
		printf("          Options: phusion2 merge <-kmer 31> <-depth 20> <-match 2> <-set 60> clust_remap.dat clust_merge.dat \n");
		printf("\n");
		printf("          fmate    output read files in patches:\n");
		printf("          Options: phusion2 fmate <-kmer 31> <-edge 100> <-split 1> clust_merge.dat reads_0*.fastq\n");
		printf("\n");
		printf("Note: for all default values, the use of any parameters is optional.\n");
		exit(1);
	}

        for(i=1;i<argc;i++)
        {
           if(!strcmp(argv[i],"-kmer"))
           {
             sscanf(argv[++i],"%d",&SEG_LEN);
           }
        } 
	ssaha_init(argv, argc, SEG_LEN);
        return EXIT_SUCCESS;
}
