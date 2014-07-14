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
* parse_cross_genome.c - Code modified for Illumina Clone Assembly Pipeline
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

static char **S_Name,**R_Name,**T_Name;
static int *hit_score,*hit_rcdex,*hit_locus1,*hit_locus2,*readlength,*superlength;
static int *hit_read1,*hit_read2,*hit_length;
static float *hit_identy;

/* forward declarations */
void Clean_Process(char **argv,int args,int nSeq);

int main(int argc, char **argv)
{
    FILE *namef;
    int i,nSeq,args;
    int n_contig,n_reads,n_readsMaxctg,nseq;
    fasta *seq;
    char line[2000]={0},tempc1[60],RC[1];

    seq=NULL;
    args=1;
      printf("here: %d\n",argc);
      printf("www: %s\n",argv[args]);
    if(argc < 3)
    {
      printf("Usage: %s [-edge 2000] <cross_genome_ouput file>\n",argv[0]);

      exit(1);
    }

    nSeq=0;

    nseq=0;
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: args \n");
      exit(1);
    }
    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      nseq++;
    }
    fclose(namef); 

    if((hit_score = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_score\n");
      exit(1);
    }
    if((hit_rcdex = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_rcdex\n");
      exit(1);
    }
    if((hit_read1 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_read1\n");
      exit(1);
    }
    if((hit_read2 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_read2\n");
      exit(1);
    }
    if((hit_locus1 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus1\n");
      exit(1);
    }
    if((hit_locus2 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus2\n");
      exit(1);
    }
    if((hit_length = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_length\n");
      exit(1);
    }
    if((readlength = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - readlength\n");
      exit(1);
    }
    if((superlength = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - superlength\n");
      exit(1);
    }
    if((hit_identy = (float *)calloc(nseq,sizeof(float))) == NULL)
    {
      printf("fmate: calloc - hit_identy\n");
      exit(1);
    }

    nSeq=nseq;
    R_Name=cmatrix(0,nseq+10,0,Max_N_NameBase);
    S_Name=cmatrix(0,nseq+10,0,Max_N_NameBase);
    T_Name=cmatrix(0,nseq+10,0,6);
    n_readsMaxctg=0;
    n_contig=0;
    n_reads=0;

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("name: %s\n",argv[args]);
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    i=0;
    while(fscanf(namef,"%s %d %s %s %d %d %d %d %s %d %f %d %d",tempc1,&hit_score[i],R_Name[i],S_Name[i],&hit_read1[i],&hit_read2[i],&hit_locus1[i],&hit_locus2[i],RC,&hit_length[i],&hit_identy[i],&readlength[i],&superlength[i])!=EOF)
    {
        if(RC[0] == 'F')
          hit_rcdex[i]=0;
        else
	{
          int hit1 = hit_locus1[i];
	  int hit2 = hit_locus2[i];
	  hit_locus1[i] = superlength[i] - hit1 +1;
	  hit_locus2[i] = superlength[i] - hit2 +1;
          hit_rcdex[i]=1;
	}
	strncpy(T_Name[i],S_Name[i],4);
        i++;
    }
    fclose(namef);

    n_reads=i;
    printf("reads: %d %s\n",n_reads,argv[args]);
    Clean_Process(argv,args,n_reads);

    printf("Job finished for %d reads!\n",nSeq);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void Clean_Process(char **argv,int args,int nSeq)
/* =============================== */
{
     int i,j,k,n;
     FILE *namef;
     int stopflag;
     void ArraySort_Mix(int n, long *arr, int *brr);

     if((namef = fopen(argv[args+1],"w")) == NULL)
     {
       printf("ERROR main:: reads group file: %s \n",argv[args]);
       exit(1);
     }

     k = 0;
     for(i=0;i<nSeq;i++)
     {
        stopflag=0;
        j=i+1;
        while((j<nSeq)&&(stopflag==0))
        {
          if(strcmp(R_Name[i],R_Name[j])==0)
          {
            j++;
          }
          else
            stopflag=1;
        }
	k = j-1;
        if(((j-i)>=2)&&(readlength[i]<superlength[i])) 
        {
	  if(superlength[i]==superlength[k])
	  {
	    int gap = readlength[i]-(hit_read2[k]-hit_read1[i]);
	    if((hit_read1[i] < 4000)&&(gap <= 5000))
              fprintf(namef,"%s\n",R_Name[i]);
//              fprintf(namef,"contig_name: %s %d %d %d\n",R_Name[i],readlength[i],hit_read2[k]-hit_read1[i],superlength[i]);
	  }
          for(n=i;n<j;n++)
          {
          }
        }
        else
        {
        }
        i=j-1;
     }

}

