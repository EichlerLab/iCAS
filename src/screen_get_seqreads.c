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
* screen_get_seqreads.c - Code modified for Illumina Clone Assembly Pipeline
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
#include "array_sort.h"
#include "matrix.h"

/* SSAS default parameters   */
static int fastq_flag =1;
static int name_flag =0;

int main(int argc, char **argv)
{
    FILE *fp,*namef,*fpOutfast;
    long totalBases;
    int i,j,k,nSeq,args,nRead,n_hits;
    int n_reads,nseq,n_input;
    fasta *sub;
    int *ctg_head,*ctg_list,*ctg_index,*hit_list,*hit_head;
    char line[2000]={0};
    int indexq,indexs,num_ctgq,num_ctgs;
    int *hit_indexs,*hit_indexq,*hit_contigs,*hit_index2,*ctg_lengthq,*ctg_lengths,*hit_locus1,*hit_locus2,*hit_rcdex2,*hit_mask;
    fasta *segg;
    long Size_q_pdata;
    int num_seqque;
    char *pdata;
    char temp[100];
    int num_temp;
    fasta *seq;
    int max_index = 0, this_index;

    seq=NULL;
    fflush(stdout);
    system("ps aux | grep sreads; date");
    if(argc < 2)
    {
      printf("Usage: %s <merge_link_file> <input_ref_fastq_file> <input_merge_fastq_file> <output_fastq_file> \n",argv[0]);
      exit(1);
    }

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-fastq"))
       {
         sscanf(argv[++i],"%d",&fastq_flag); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-name"))
       {
         sscanf(argv[++i],"%d",&name_flag);
         args=args+2;
       }
    }

    n_input=0;
    nseq=0;
    
    if ((fp = fopen(argv[args + 1], "rb")) == NULL) {
    	 printf("Error: Cannot open file %s\n", argv[args + 1]);
	 exit(1);
    }
    
    fseek(fp, 0, SEEK_END);
    Size_q_pdata = ftell(fp) + 1;
    fclose(fp);
    
    if ((pdata = (char*)calloc(Size_q_pdata, sizeof(char))) == NULL) {
        printf("Error: calloc pdata\n");
	exit(1);
    }
    
    num_seqque = extractFastq(argv[args+1],pdata,Size_q_pdata);
    
    if ((segg = (fasta*)calloc((num_seqque+1), sizeof(fasta))) == NULL) {
      	printf("Error: calloc segg\n");
	exit(1);
    }
    
    if ((seq = decodeFastq(argv[args + 1], &num_seqque, &totalBases, pdata, Size_q_pdata,segg)) == NULL) {
    	printf("Error: no query data found.\n");
	exit(1);
    }
    
    nSeq = num_seqque;
    fastaUC(seq,nSeq);
    num_ctgq = nSeq;

    if ((fp = fopen(argv[args + 2],"rb")) == NULL) {
    	printf("Error: Cannot open file %s\n",argv[args + 2]);
	exit(1);
    }
    
    fseek(fp, 0, SEEK_END);
    Size_q_pdata = ftell(fp) + 1;
    fclose(fp);
    
    if ((pdata = (char*)calloc(Size_q_pdata, sizeof(char))) == NULL) {
    	printf("Error: calloc pdata\n");
	exit(1);
    }
    
    num_seqque = extractFastq(argv[args+2],pdata,Size_q_pdata);

    if((segg = (fasta*)calloc((num_seqque + 2),sizeof(fasta))) == NULL) {
    	printf("Error: calloc segg\n");
	exit(1);
    }

    if ((sub = decodeFastq(argv[args + 2], &num_seqque, &totalBases, pdata, Size_q_pdata, segg)) == NULL) {
    	printf("Error: no query data found.\n");
	exit(1);
    }
    
    nSeq = num_seqque;
    fastaUC(sub,nSeq);
    nRead=nSeq;
    num_ctgs = nSeq;

    if ((namef = fopen(argv[args],"r")) == NULL) {
    	fprintf(stderr, "ERROR main:: reads group file %s\n", argv[args]);
      	exit(1);
    }
    
    n_reads = 0;
    
    while (!feof(namef)) {
    	fscanf(namef, "%*s %*s %d %*d %*d %*d %*d %*d %*d %*d", &this_index);
	
	if (this_index > max_index) max_index = this_index;
 	
      	if(feof(namef)) break;
	
      	n_reads++;
    }
    
    fclose(namef); 

    if((hit_contigs= (int *)calloc(n_reads,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - hit_contigs\n");
       exit(1);
    }
    if((hit_index2= (int *)calloc(n_reads,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - hit_index2\n");
       exit(1);
    }
    if((hit_indexs= (int *)calloc(n_reads,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - hit_indexs\n");
       exit(1);
    }
    if((hit_indexq= (int *)calloc(n_reads,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - hit_indexq\n");
       exit(1);
    }
    if((hit_rcdex2= (int *)calloc(n_reads,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - hit_rcdex2\n");
       exit(1);
    }
    if((ctg_lengthq= (int *)calloc(n_reads,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - hit_contigs\n");
       exit(1);
    }
    if((ctg_lengths= (int *)calloc(n_reads,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - ctg_lengths\n");
       exit(1);
    }
    if((hit_mask= (int *)calloc(n_reads,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - ctg_mask\n");
       exit(1);
    }
    
    if ((hit_list = (int *)calloc(max_index + 1, sizeof(int))) == NULL) {
    	fprintf(stderr, "ERROR Out of memory: calloc - hit_list\n");
        exit(1);
    }
    
    if((hit_head= (int *)calloc(n_reads,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - ctg_head\n");
       exit(1);
    }
    if((ctg_list= (int *)calloc(num_ctgs,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - ctg_list\n");
       exit(1);
    }
    if((ctg_head= (int *)calloc(num_ctgs,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - ctg_head\n");
       exit(1);
    }
    if((hit_locus1= (int *)calloc(n_reads,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - hit_contigs\n");
       exit(1);
    }
    if((hit_locus2= (int *)calloc(n_reads,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - hit_locus2\n");
       exit(1);
    }

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR Memory_Allocate:: reads group file \n");
      exit(1);
    }

    i = 0;
    printf("reads: %d %d %d %s %s\n",n_reads,nRead,nSeq,argv[args],argv[args+1]);
    
    while (fscanf(namef, "%s %s %d %d %d %d %d %d %d %d", temp, temp, &hit_indexq[i], &hit_contigs[i], &hit_indexs[i], &num_temp,
     	&hit_locus1[i], &hit_locus2[i], &ctg_lengthq[i], &ctg_lengths[i]) != EOF) {
	
    	hit_list[hit_indexq[i]]++;
	
      	if (hit_contigs[i] == 1) {
            ctg_list[hit_indexs[i]]++;
      	}
	
	i++;
    }
    
    fclose(namef);

    n_hits = 0;
    for(i=0;i<n_reads;i++)
    {
       int stopflag;

       stopflag=0;
       j=i+1;
       while((j<n_reads)&&(stopflag==0))
       {
          if(hit_indexq[i]==hit_indexq[j])
          {
            j++;
          }
          else
            stopflag=1;
       }
       for(k=i;k<j;k++)
       {
          if(k==i)
            hit_index2[k] = 1;
          if(k==(j-1))
            hit_index2[k] = 1;
          hit_head[k] = j-i;
    	printf("merge: %d %d %d\n",hit_head[k],k,hit_indexq[k]);
       }
       hit_list[n_hits] = j-i;
       i = j-1;
       n_hits++;
    }

    ctg_head[0] = 0;
    max_index = 0;
    
    for(i = 1; i < num_ctgs; i++) {
    	ctg_head[i] = ctg_head[i-1] + ctg_list[i-1];
	
	if (ctg_head[i] > max_index) max_index = ctg_head[i];
    }
    
    

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR Memory_Allocate:: reads group file \n");
      exit(1);
    }

    for(i=0;i<num_ctgs;i++)
       ctg_list[i] = 0;

    if ((ctg_index = (int *)calloc(max_index + n_reads + 1, sizeof(int))) == NULL) {
    	fprintf(stderr, "ERROR out of memory: calloc - ctg_index\n");
        exit(1);
    }
    
    i = 0;
    while(fscanf(namef,"%s %s %d %d %d %s %s %s %s %s",temp,temp,&indexq,&k,&indexs,temp,temp,temp,temp,temp)!=EOF)
    {
      if(k==1)
      {
        ctg_index[ctg_head[indexs]+ctg_list[indexs]] = i;
        ctg_list[indexs]++;
      }
      hit_mask[i] = 1;
      i++;
    }
    fclose(namef);


    for(i=0;i<num_ctgs;i++)
    {
       if(ctg_list[i]>1)
       {
         int max_ctg =0;
         int max_ctgid = 0;

         for(j=0;j<ctg_list[i];j++)
         {
            int idd = ctg_head[i]+j;
            int idt = ctg_index[idd];
            if(hit_head[idt] >= max_ctg)
            {
              max_ctg = hit_head[idt];
              max_ctgid = idt;
              printf("www: %d %d %d %d %d %d %d\n",hit_index2[idt],ctg_list[i],hit_head[idt],hit_locus1[idt],hit_locus2[idt],ctg_lengthq[idt],ctg_lengths[idt]);
            } 
         }
         hit_mask[max_ctgid] = 1;
         for(j=0;j<ctg_list[i];j++)
         {
            int idd = ctg_head[i]+j;
            int idt = ctg_index[idd];
            if(idt!=max_ctgid)
            {
              hit_mask[idt] = 0;
              printf("kickout: %d %d %d %d %d %d %d\n",hit_index2[idt],ctg_list[i],hit_head[idt],hit_locus1[idt],hit_locus2[idt],ctg_lengthq[idt],ctg_lengths[idt]);
            }
         }
       }
    }

    if((fpOutfast = fopen(argv[args+3],"w")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    n_reads=0;
    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      if(hit_mask[n_reads] == 1)
       fprintf(fpOutfast,"%s",line);
      n_reads++;
    }
    fclose(namef); 

    fclose(fpOutfast);
    return(0);

}
/* end of the main */

