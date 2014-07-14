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
* phusion2.c - Code modified for Illumina Clone Assembly Pipeline                  *
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

static B64_long *hist;
static int N_SET=60000;
static int lsize=4000;
static int mhist=16;
static int mhist2=2;
static int r_len=220;
static int i_start=0;
static int i_end=500;
static int min_size = 10;
static int num_ukmer = 10;
static int num_cut = 10;
static int num_copy = 10000;
static int i_split=0;
static int i_patch=0;
static int i_edge=50;
static int nmatch=6;
static int max_edge=200;
static int num_reads=400000000;
static char *SCGname=NULL;//SingleCopyGenome
static int n_contig;
static int n_group = 30000;
static int n_block;
static char plate_name[10];
static int rshift = 6;

void Hash_Process(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);
void Random_Process(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);
void Sort_Process(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);
void Screen_Process(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);
void Edge_Process(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);
void Matrix_Process(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);
void Clust_Process(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);
void Remap_Process(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);
void Merge_Process(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);
void Fmate_Process(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);
void Extract_Process(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);
void Name_Process(int brr,int len,int crr,int m1,int m2,int step,int pmod, char **argv, int argc, int args);
unsigned int **uimatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
B64_long **limatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
void ArraySort_Mix(int n, B64_long *arr, int *brr);
int **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);


/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void fastaSort(char **argv, int argc, int SEG_LEN)
/* =============================================  */
{
     B64_long i,n_Sbase=SEG_LEN;
     int args,ac;
     int m_st,m_ed,step_len = 2,pmod;
     int nseq = 0;
     double Qerr[100];
     Qerr[0] = 1.0;

     for(i=1;i<100;i++) Qerr[i] = pow((double)10.,(double)-.1*i);

/*   sort all the names of genes or name entries   */
     printf("Input data starts ...\n");
     args=2;
     for(i=1;i<argc;i++)
     {
       if(!strcmp(argv[i],"-kmer"))
       {
	 i++;
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-depth"))
       {
         sscanf(argv[++i],"%d",&mhist);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-step"))
       {
         sscanf(argv[++i],"%d",&step_len);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-set"))
       {
         sscanf(argv[++i],"%d",&N_SET);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-length"))
       {
         sscanf(argv[++i],"%d",&r_len);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-shift"))
       {
         sscanf(argv[++i],"%d",&rshift);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-match"))
       {
         sscanf(argv[++i],"%d",&nmatch);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-plate"))
       {
         sscanf(argv[++i],"%s",plate_name);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-min"))
       {
         sscanf(argv[++i],"%d",&min_size);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-uniq"))
       {
         sscanf(argv[++i],"%d",&num_ukmer);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-copy"))
       {
         sscanf(argv[++i],"%d",&num_copy);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-cut"))
       {
         sscanf(argv[++i],"%d",&num_cut);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-split"))
       {
         sscanf(argv[++i],"%d",&i_split);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-size"))
       {
         sscanf(argv[++i],"%d",&lsize);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-group"))
       {
         sscanf(argv[++i],"%d",&n_group);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-start"))
       {
         sscanf(argv[++i],"%d",&i_start);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-end"))
       {
         sscanf(argv[++i],"%d",&i_end);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-matrix"))
       {
         sscanf(argv[++i],"%d",&max_edge);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-patch"))
       {
         sscanf(argv[++i],"%d",&i_patch);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-edge"))
       {
         sscanf(argv[++i],"%d",&i_edge);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-reads"))
       {
         sscanf(argv[++i],"%d",&num_reads);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-SCG"))
       {
	 SCGname = argv[++i];
         args=args+2;
       } 
     }
     for(ac=0;ac<argc;ac++)
	printf("%s ",argv[ac]);
	
     printf("\n");
     printf("kmer size:         %ld\n",n_Sbase);
     printf("depth:             %d\n",mhist);
     printf("ph2 depth:         %d\n",mhist2);
     printf("matrix size:       %d\n",max_edge);

     if((hist = (B64_long *)calloc(400,sizeof(B64_long))) == NULL)
     {
        printf("ERROR ssaha_init: calloc - hist\n");
        exit(1);
     }

     m_st=nmatch;
     m_ed=0;
     pmod=1;

     nseq = 0;
     if(strcmp(argv[1],"hash")==0)
       Hash_Process(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);
     else if(strcmp(argv[1],"random")==0)
       Random_Process(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);
     else if(strcmp(argv[1],"sort")==0)
       Sort_Process(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);
     else if(strcmp(argv[1],"screen")==0)
       Screen_Process(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);
     else if(strcmp(argv[1],"edge")==0)
       Edge_Process(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);
     else if(strcmp(argv[1],"matrix")==0)
       Matrix_Process(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);
     else if(strcmp(argv[1],"clust")==0)
       Clust_Process(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);
     else if(strcmp(argv[1],"remap")==0)
       Remap_Process(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);
     else if(strcmp(argv[1],"merge")==0)
       Merge_Process(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);
     else if(strcmp(argv[1],"fmate")==0)
       Fmate_Process(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);
     else if(strcmp(argv[1],"name")==0)
       Name_Process(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);
     else if(strcmp(argv[1],"extract")==0)
       Extract_Process(nseq,SEG_LEN,mhist,m_st,m_ed,step_len,pmod,argv,argc,args);
     else
     {
       printf("Invalid input command!\n");
       exit(1);
     }
     printf("All jobs finished\n");
     fflush(stdout);
     system("ps aux |date");
     return;
}


int ssaha_init( char **argv, int argc, int SEG_LEN)
{

    fastaSort(argv,argc,SEG_LEN);
    return(1);   
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Hash_Process(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i,j,k,iseq,n_Sbase=SEG_LEN;
     B64_long ns=0,IntSeg=0,IntBase=0;
     char *b;
     fasta *seqp;
     fasta *seq;
     char nameout_hash[100],nameout_idex[100],nameout_size[100],nameout_list[100];
     int nSeq;
     B64_long totalBases,nsize[10];
     B64_long mhistc = 0;
     B64_long mhistcc = 0;
     int nsorts = 1024;
     int nshift = (n_Sbase<<1)-10;//2^10=nsorts
     B64_long offset,sizesm;
     int ac,max_nhit;
     int seqc,num_refseq;
     unsigned B64_long seqcharc = 0;
     B64_long *patch_array,*patch_head;
     int *patch_index,*patch_list; 
     int **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     unsigned int **uimatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     B64_long **limatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     void ArraySort_Mix(int n, B64_long *arr, int *brr);
     B64_long kmask = (1L<<(n_Sbase<<1))-1;
     double Qerr[100];
     B64_long gtBases=0;
     Qerr[0] = 1.0;
     FILE *fp,*fpHASH,*fpIDEX,*fpLIST,*fpSIZE;
     fasta *segg;
     B64_long Size_q_pdata;
     int num_seqque;
     char *pdata;

     if((fp=fopen(argv[args],"rb"))==NULL) {
     	printf("Error: Cannot open file %s\n", argv[args]);
	exit(1);
     }
     
     fseek(fp, 0, SEEK_END);
     Size_q_pdata = ftell(fp) + 1;
     fclose(fp);
     
     if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL) {
     	printf("Error: calloc pdata\n");
	exit(1);
     }
     
     num_seqque = extractFastq(argv[args],pdata,Size_q_pdata);
     
     if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL) {
        printf("Error: calloc segg\n");
	exit(1);
     }
     
     if((seq=decodeFastq(argv[args],&num_seqque,&totalBases,pdata,Size_q_pdata,segg))==NULL) {
     	printf("Error: no query data found.\n");
	exit(1);
     }
     
     nRead = num_seqque;
     nSeq = num_seqque;
     fastaUC(seq,nRead);
     printf("Number of reads: %d \n",nSeq);
     printf("Number of bases: %ld \n",totalBases);

     mhistc = 0;
     gtBases = 0;
     nsorts = num_reads;
     for(i=0;nsorts > 10;i++)
     {
        nsorts = nsorts>>1;
     } 
     nsorts = 1L<<i;
     nshift = (n_Sbase<<1)-i;//2^i=nsorts
     if((patch_list= (int *)calloc(nsorts,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - patch_list\n");
       exit(1);
     }
     if((patch_head= (B64_long *)calloc(nsorts,sizeof(B64_long))) == NULL)
     {
       printf("phusion: calloc - patch_head\n");
       exit(1);
     }
     num_refseq = 0;
     for(ac=(args);ac<argc;ac++)
     {
        for(iseq=0;iseq<nSeq;iseq++)
        {
           double e=0.0,eh[32];
 	   B64_long IntBaseRC;
           B64_long IntSegRC=0;
	   int pos=0;
           char *q;

           IntSeg=0;
           seqp=seq+iseq;
           b = (char*)seqp->data;
           if(seqp->finished) {
             q = NULL;
           } else {
             q = NULL; //(char*)seqp->qual;
           }

           ns=(seqp->length);
           seqcharc += strlen(seqp->name);
           k=n_Sbase;
	   if(k > ns) k=ns;

           for(;k>0;k--,b++,pos++)
           {
              if     (*b == 'A') IntBase=0;
              else if(*b == 'C') IntBase=1;
              else if(*b == 'G') IntBase=2;
              else if(*b == 'T') IntBase=3;
              else
              {
	        e += eh[pos] = 1.;
	        continue;
              }
              e += eh[pos] = (q != NULL) ? Qerr[q[pos]] : 0;
	      IntBaseRC=(IntBase^3)<<((n_Sbase-k)<<1);
              IntBase=IntBase<<((k-1)<<1);
              IntSeg+=IntBase;
	      IntSegRC+=IntBaseRC;
           }
           for(j=ns-n_Sbase;--j>=0;b++,pos++)
           {
	      int p=pos%n_Sbase;
	      if(e < .99) 
              {
	        if(IntSeg > IntSegRC) 
                {
                  int itt = IntSegRC>>nshift;
                  patch_list[itt]++;
	        } 
                else 
                {
                  int itt = IntSeg>>nshift;
                  patch_list[itt]++;
	        }
	        gtBases++;
	      }
	      e -= eh[p];
	      IntSeg = (IntSeg<<2)&kmask;
	      IntSegRC = (IntSegRC>>2)&kmask;
              if     (*b == 'A') IntBase=0;
              else if(*b == 'C') IntBase=1;
              else if(*b == 'G') IntBase=2;
              else if(*b == 'T') IntBase=3;
              else
              {
	        e += eh[p] = 1.;
	        continue;
              }
	      e += eh[p] = (q != NULL) ? Qerr[q[pos]] : 0; 
	      IntBaseRC=(IntBase^3)<<((n_Sbase-1)<<1);
              IntSeg+=IntBase;
	      IntSegRC+=IntBaseRC;
           }

	   if(e < .99) 
           {
	     if(IntSeg > IntSegRC) 
             {
               int itt = IntSegRC>>nshift;
               patch_list[itt]++;
	     } 
             else 
             {
               int itt = IntSeg>>nshift;
               patch_list[itt]++;
	     }
	     gtBases++;
           }
        }
     }

     printf("data reading finished ...\n");
     fflush(stdout);
     system("ps aux | date");


/*     memset(nameout_list,'\0',100);
     printf("Output list file ... %d\n",nsorts);
     sprintf(nameout_list,"%s.list",argv[args]);
     fpLIST=fopen(nameout_list,"wb");
     offset = 0;
     fseek(fpLIST,offset,0);
     sizesm=nsorts;
     fwrite(patch_list,sizeof(int),sizesm,fpLIST);
     fclose(fpLIST);

     exit(1);      */
     printf("Input data finished one set, nsorts=%d, totalgoodBP=%ld\n",nsorts,gtBases);
     fflush(stdout);
     system("ps aux | grep HASH; date");
     max_nhit=0;
     patch_head[0]=0;
     for(i=0;i<nsorts;i++) 
     {
        if(patch_list[i]>max_nhit)
          max_nhit = patch_list[i];
        if(i>0)
          patch_head[i] = patch_head[i-1] + patch_list[i-1];
//        patch_list[i] = 0;
     }
     memset(patch_list,0,4*nsorts);  
     mhistcc=patch_head[nsorts-1]+patch_list[nsorts-1]+max_nhit+10;

     printf("setting array memory: %ld\n",mhistcc);
     fflush(stdout);
     system("ps aux | grep HASH; date");
     if((patch_array= (B64_long *)calloc(mhistcc,sizeof(B64_long))) == NULL)
     {
       printf("ssaha: calloc - patch_array\n");
       exit(1);
     }
     if((patch_index= (int *)calloc(mhistcc,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }

     fflush(stdout);
     system("ps aux | grep HASH; date");

     seqc = 0;

     mhistcc=0;

     for(iseq=0;iseq<nSeq;iseq++,seqc++)
     {
        char *q;
        double e=0,eh[32];
        B64_long IntBaseRC;
        B64_long IntSegRC=0;
        B64_long IntSeg=0;
        int pos=0;

        seqp=seq+iseq;
        b = (char*)seqp->data; 
        if(seqp->finished) {
          q = NULL;
        } else {
          q = NULL; //(char*)seqp->qual;
        } 
        ns=(seqp->length);

        for(k=n_Sbase;k>0;k--,b++,pos++)
        {
           if     (*b == 'A') IntBase=0;
           else if(*b == 'C') IntBase=1;
           else if(*b == 'G') IntBase=2;
           else if(*b == 'T') IntBase=3;
           else
           {
	     e += eh[pos] = 1.;
	     continue;
           }
           e += eh[pos] = (q != NULL) ? Qerr[q[pos]] : 0;
	   IntBaseRC=(IntBase^3)<<((n_Sbase-k)<<1);
           IntBase=IntBase<<((k-1)<<1);
           IntSeg+=IntBase;
	   IntSegRC+=IntBaseRC;
        }
        for(j=ns-n_Sbase;--j>=0;b++,pos++)
        {
	   int p=pos%n_Sbase;
	   if(e < .99) 
           {
             B64_long patch_pos;
    	     if(IntSeg > IntSegRC) 
             {
               int itt = IntSegRC>>nshift;
               patch_pos = patch_head[itt]+patch_list[itt];
               patch_array[patch_pos] = IntSegRC;
               patch_index[patch_pos]=seqc;
               patch_list[itt]++;
	       mhistcc++;
	     } 
             else 
             {
               int itt = IntSeg>>nshift;
               patch_pos = patch_head[itt]+patch_list[itt];
               patch_array[patch_pos] = IntSeg;
               patch_index[patch_pos]=seqc;
               patch_list[itt]++;
	       mhistcc++;
	     }
	     gtBases++;
	   }
	   e -= eh[p];
	   IntSeg = (IntSeg<<2)&kmask;
	   IntSegRC = (IntSegRC>>2)&kmask;
           if     (*b == 'A') IntBase=0;
           else if(*b == 'C') IntBase=1;
           else if(*b == 'G') IntBase=2;
           else if(*b == 'T') IntBase=3;
           else
           {
	     e += eh[p] = 1.;
             continue;
           }
	   e += eh[p] = 0;
	   IntBaseRC=(IntBase^3)<<((n_Sbase-1)<<1);
           IntSeg+=IntBase;
	   IntSegRC+=IntBaseRC;
        }
        if(e < .99) 
        {
          B64_long patch_pos;
	  if(IntSeg > IntSegRC) 
          {
            int itt = IntSegRC>>nshift;
            patch_pos = patch_head[itt]+patch_list[itt];
            patch_array[patch_pos] = IntSegRC;
            patch_index[patch_pos]=seqc;
            patch_list[itt]++;
	    mhistcc++;
	  } 
          else 
          {
            int itt = IntSeg>>nshift;
            patch_pos = patch_head[itt]+patch_list[itt];
            patch_array[patch_pos] = IntSeg;
            patch_index[patch_pos]=seqc;
            patch_list[itt]++;
	    mhistcc++;
	  }
	  gtBases++;
        }
     }
     memset(nameout_hash,'\0',100); 
     memset(nameout_idex,'\0',100);
     sprintf(nameout_hash,"%s.hash",argv[args]); 
     sprintf(nameout_idex,"%s.idex",argv[args]);
     
     if ((fpHASH=fopen(nameout_hash,"wb")) == NULL) {
     	printf("Error: can't write hash table to %s.\n", nameout_hash);
	exit(1);
     }

     printf("Output hashtable file ... %ld kmers\n",mhistcc);
     offset = 0;
     fseek(fpHASH,offset,0);
     sizesm=mhistcc;
     fwrite(patch_array,sizeof(B64_long),sizesm,fpHASH);
     fclose(fpHASH);

     printf("Output read index file ...\n");
     
     if ((fpIDEX=fopen(nameout_idex,"wb")) == NULL) {
     	printf("Error: can't write index file to %s.\n", nameout_idex);
	exit(1);
     }
     
     offset = 0;
     fseek(fpIDEX,offset,0);
     sizesm=mhistcc;
     fwrite(patch_index,sizeof(int),sizesm,fpIDEX);
     fclose(fpIDEX);

     memset(nameout_list,'\0',100);
     printf("Output list file ... %d\n",nsorts);
     sprintf(nameout_list,"%s.list",argv[args]);
     
     if ((fpLIST=fopen(nameout_list,"wb")) == NULL) {
     	printf("Error: can't write list file to %s\n.", nameout_list);
	exit(1);
     }
      
     offset = 0;
     fseek(fpLIST,offset,0);
     sizesm=nsorts;
     fwrite(patch_list,sizeof(int),sizesm,fpLIST);
     fclose(fpLIST);

     memset(nameout_size,'\0',100);
     printf("Output size file ... %d\n",nSeq);
     sprintf(nameout_size,"%s.size",argv[args]);

     if ((fpSIZE=fopen(nameout_size,"wb")) == NULL) {
     	printf("Error: can't write size file to %s\n.", nameout_size);
	exit(1);
     }
      
     offset = 0;
     fseek(fpSIZE,offset,0);
     sizesm=2;
     nsize[0] = nSeq;
     nsize[1] = mhistcc;
     nsize[2] = mhistcc;
     fwrite(nsize,sizeof(B64_long),sizesm,fpSIZE);
     fclose(fpSIZE);

     free(patch_list);
     free(patch_head);

     free(patch_array);
     free(patch_index);

}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Random_Process(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i,j,iseq;
     int seq_st,seq_ed,rc;
     fasta *seqp;
     fasta *seq;
     int nSeq,num_gaps,*reads_mask,*reads_link,*reads_gaps;
     B64_long totalBases,dsize_4;
     FILE *fp,*fpOutfast;
     fasta *segg;
     B64_long Size_q_pdata;
     int num_seqque;
     char *pdata;

     if((fp=fopen(argv[args],"rb"))==NULL) {
     	printf("Error: Cannot open file %s\n", argv[args]);
	exit(1);
     }
      
     fseek(fp, 0, SEEK_END);
     Size_q_pdata = ftell(fp) + 1;
     fclose(fp);
     
     if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL) {
     	printf("Error: calloc pdata\n");
	exit(1);
     }
     
     num_seqque = extractFastq(argv[args],pdata,Size_q_pdata);
     
     if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL) {
       printf("Error: calloc segg\n");
       exit(1);
     }
     
     if((seq=decodeFastq(argv[args],&num_seqque,&totalBases,pdata,Size_q_pdata,segg))==NULL) {
       printf("Error: no query data found.\n");
       exit(1);
     }
     
     nRead = num_seqque;
     nSeq = num_seqque;
     fastaUC(seq,nRead);
     printf("Number of reads: %d \n",nSeq);
     printf("Number of bases: %ld \n",totalBases);

     if((reads_mask= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - reads_mask\n");
       exit(1);
     }
     if((reads_link= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - reads_link\n");
       exit(1);
     }
     if((reads_gaps= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - reads_gaps\n");
       exit(1);
     }

     dsize_4 = nSeq;
     dsize_4 = 4*dsize_4;
     memset(reads_mask,-1,dsize_4);
     for(iseq=0;iseq<nSeq;iseq++)
     {
        if(iseq%10 == 0)
	  continue;
	i = (int)(drand48()*nRead);
	while(reads_mask[i]>=0)
	{
	  i = (int)(drand48()*nRead);
	}
	reads_mask[i] = iseq;
	reads_link[iseq] = i;
     }
     num_gaps = 0;
     for(j=0;j<nSeq;j++)
     {
        if(reads_mask[j]<0)
	{
	  reads_gaps[num_gaps] = j;
	  num_gaps++;
	}
     }
//     num_gaps = 0;
     for(j=0;j<nSeq;j++)
     {
        if(j%10 == 0)
	{
	  reads_link[j] = reads_gaps[num_gaps-1];
	  num_gaps--;
	}
     }


/*   output fastq file   */
     if((fpOutfast = fopen(argv[args+1],"w")) == NULL)
     {
       printf("ERROR main:: reads group file: %s \n",argv[args]);
       exit(1);
     }
     for(iseq=0;iseq<nRead;iseq++)
     {
        int idk = reads_link[iseq];
        seqp=seq+idk;
        seq_st = 0;
	seq_ed = seqp->length;
        fprintf(fpOutfast,"@%s\n",seqp->name);
        for(rc=seq_st;rc<seq_ed;rc++)
           fprintf(fpOutfast,"%c",seqp->data[rc]);
        fprintf(fpOutfast,"\n+\n");
        for(rc=seq_st;rc<seq_ed;rc++)
           putc(seqp->qual[rc]+041,fpOutfast);
        fprintf(fpOutfast,"\n");
     }
     fclose(fpOutfast);
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Extract_Process(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i,k;
     fasta *seqp;
     fasta *seq;
     char nameout_hash[100],nameout_size[100];
     B64_long tKmer,nsize[10];
     B64_long totalBases;
     B64_long offset,sizesm;
     int ac;
     int tRead,iMat;
     int i_reads,n_reads,c_reads,out_flag; //id_read will never exceed mhist
     int *read_index;
     void ArraySort_Mix(int n, B64_long *arr, int *brr);
     char *ptr;
     unsigned char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     FILE *namef,*fp,*fpHASH = NULL,*fpSIZE;
     fasta *segg;
     B64_long Size_q_pdata;
     int num_seqque;
     char *pdata;

     tKmer = 0;
     tRead = 0;
     memset(nameout_size,'\0',100);
     sprintf(nameout_size,"matrix.size");
     
     if ((fpSIZE = fopen(nameout_size, "rb")) == NULL) {
        printf("Error: can't open file %s.\n", nameout_size);
        exit(1);
     }

     offset = 0;
     fseek(fpSIZE,offset,0);
     sizesm=2;
     fread(nsize,sizeof(B64_long),sizesm,fpSIZE);
     tKmer = nsize[0];
     tRead = nsize[1];
     printf("reads: %ld %d %s\n",tKmer,tRead,argv[args]);

     if((read_index = (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - reads_index\n");
       exit(1);
     }

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }

/*   read the phusion output file         */
     i_reads=0;
     n_reads=0;
     while(!feof(namef))
     {
       int nPair=0;
       char line2[200],line[200],base[200];
      
       fgets(line,200,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       if((strncmp(line,"readnames",9))==0)
       {
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
	 {
	 }
	 i=0;
	 for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
	 {
	    if(i==5)
	    {
	      memset(base,'\0',200);
	      strcat(base,ptr);
	      c_reads = atoi(base);
	    }
	 }
	 n_reads = n_reads+c_reads;
       }
       else
       {      
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==0)
            {
              int idt;
              memset(base,'\0',500);
              strcat(base,ptr);
              idt = atoi(base);
              read_index[idt] = 1;
              i_reads++;
            }
         }
       }
     }

     printf("reads: %d %d %d\n",tRead,n_reads,i_reads);
     nRead = tRead - n_reads;

     n_reads = 0;
     i_reads = 0;
     iMat = 0;
     out_flag=1;
     for(ac=(args+1);ac<argc;ac++)
     {
        int len;
	
        if((fp=fopen(argv[ac],"rb"))==NULL) {
	    printf("Error: Cannot open file %s\n",argv[ac]);
	    exit(1);
	}
        
	fseek(fp, 0, SEEK_END);
        Size_q_pdata = ftell(fp) + 1;
        fclose(fp);
        
	if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL) {
          printf("Error: calloc pdata\n");
	  exit(1);
	}
	
        num_seqque = extractFastq(argv[ac],pdata,Size_q_pdata);
	
        if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL) {
          printf("Error: calloc segg\n");
	  exit(1);
	}
	
        if((seq=decodeFastq(argv[ac],&num_seqque,&totalBases,pdata,Size_q_pdata,segg))==NULL) {
	    printf("Error: no query data found.\n");
	    exit(1);
	}
	
        nRead = num_seqque;
        printf("read file: %s %d %d %d\n",argv[ac],iMat,nRead,N_SET);
        for(i=0;i<nRead;i++)
        {
           if(read_index[n_reads+i] == 0)
	   {
	     if(((i_reads%N_SET)==0)&&(out_flag==1))
	     {
	       memset(nameout_hash,'\0',100);
	       sprintf(nameout_hash,"extract-%03d.fastq",iMat); 
	       
	       if ((fpHASH=fopen(nameout_hash,"w")) == NULL) {
	            printf("Error: can't write to %s.\n", nameout_hash);
		    exit(1);
	       }
	       
               printf("output file: %s %d %d %d\n",nameout_hash,n_reads,nRead,N_SET);
	       out_flag = 0;
	       iMat++;
	     }
	     seqp = seq+i;
	     len = seqp->length;
	     fprintf(fpHASH,"@%s\n",seqp->name);
	     for(k=0;k<len;k++)
	        fprintf(fpHASH,"%c",seqp->data[k]);
	     fprintf(fpHASH,"\n+\n");
	     for(k=0;k<len;k++)
	        putc(seqp->qual[k]+041,fpHASH);
             fprintf(fpHASH,"\n");
	     out_flag=1;
	     i_reads++;

	     if(((i_reads%N_SET)==0)&&(out_flag==1))
	       fclose(fpHASH);
	   }
        }
        n_reads = n_reads + nRead;
	free(pdata);
	free(segg);
     }
     fclose(fpHASH);
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Name_Process(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i;
     fasta *seqp;
     fasta *seq;
     int nSeq,ac,name_len;
     FILE *fp,*fpOutfast;
     fasta *segg;
     B64_long Size_q_pdata,totalBases;
     int num_seqque;
     char *pdata;

     name_len = strlen(plate_name);
     if((fpOutfast = fopen(argv[args],"w")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }
     for(ac=(args+1);ac<argc;ac++)
     {
        if((fp=fopen(argv[ac],"rb"))==NULL) {
	    printf("Error: Cannot open file %s\n", argv[ac]);
	    exit(1);
	}
	
        fseek(fp, 0, SEEK_END);
        Size_q_pdata = ftell(fp) + 1;
        fclose(fp);
	
        if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL) {
          printf("Error: calloc pdata\n");
	  exit(1);
	}
	
        num_seqque = extractFastq(argv[ac],pdata,Size_q_pdata);
        
	if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL) {
          printf("Error: calloc segg\n");
	  exit(1);
	}
	
        if((seq=decodeFastq(argv[ac],&num_seqque,&totalBases,pdata,Size_q_pdata,segg))==NULL) {
          printf("Error: no query data found.\n");
	  exit(1);
	}
	
        nSeq = num_seqque;
        for(i=0;i<nSeq;i++)
        {
    	   int seq_st,seq_ed,rc;

           seqp= seq + i;
           seq_st = 0; 
           seq_ed = seqp->length;
           if(strncmp(seqp->name,plate_name,name_len)==0)
           {
             fprintf(fpOutfast,"@%s\n",seqp->name);
             for(rc=seq_st;rc<seq_ed;rc++)
                fprintf(fpOutfast,"%c",seqp->data[rc]);
             fprintf(fpOutfast,"\n+\n");
             for(rc=seq_st;rc<seq_ed;rc++)
                putc(seqp->qual[rc]+041,fpOutfast);
             fprintf(fpOutfast,"\n");
           }
        }        
	free(pdata);
	free(segg);
     }
     fclose(fpOutfast);
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Fmate_Process(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i,j,k;
     fasta *seqp;
     fasta *seq;
     char nameout_hash[100],nameout_size[100],nameout_list[100];
     B64_long tKmer,nsize[10],dsize_4;
     B64_long totalBases;
     B64_long offset,sizesm;
     int ac,*map_sort,*pat_sort,map_sort2[100];
     int tRead,nMat,iMat,iBulk,stopflag;
     int *ctg_head,*ctg_list;
     unsigned char **rdname,**rdbase,**rdqual; 
     int i_reads  = 0,n_reads,c_reads; //id_read will never exceed mhist
     int *read_index,*read_length;
     void ArraySort_Mix(int n, B64_long *arr, int *brr);
     char *ptr;
     unsigned char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     FILE *namef,*fp,*fpHASH,*fpLIST,*fpSIZE;
     fasta *segg;
     B64_long Size_q_pdata;
     int num_seqque;
     char *pdata;

     tKmer = 0;
     tRead = 0;
     memset(nameout_size,'\0',100);
     sprintf(nameout_size,"matrix.size");
     
     if ((fpSIZE = fopen(nameout_size,"rb")) == NULL) {
     	printf("Error: can't open %s\n", nameout_size);
	exit(1);
     }
     
     offset = 0;
     fseek(fpSIZE,offset,0);
     sizesm=2;
     fread(nsize,sizeof(B64_long),sizesm,fpSIZE);
     tKmer = nsize[0];
     tRead = nsize[1];
     printf("reads: %ld %d\n",tKmer,tRead);

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }

/*   read the phusion output file         */
     n_contig=0;
     n_reads = 0;
     while(!feof(namef))
     {
       int nPair=0;
       char line2[200],line[200],base[200];
      
       fgets(line,200,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       if((strncmp(line,"readnames",9))==0)
       { 
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==5)
            {
              memset(base,'\0',200);
              strcat(base,ptr);
                c_reads = atoi(base);
            }
         }
         n_reads = n_reads+c_reads;
         n_contig++;
       }
     }

     printf("contig: %d %d\n",n_contig,n_reads);
     if((ctg_head = (int *)calloc(n_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_head\n");
       exit(1);
     }
     if((ctg_list = (int *)calloc(n_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_list\n");
       exit(1);
     }
     if((pat_sort = (int *)calloc(n_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - pat_sort\n");
       exit(1);
     }
     if((read_index = (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - reads_index\n");
       exit(1);
     }
     if((map_sort= (int *)calloc(1000,sizeof(int))) == NULL)
     {        
       printf("phusion: calloc - map_sort\n");
       exit(1); 
     }        
     fclose(namef);

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }

     dsize_4 = tRead;
     dsize_4 = 4*dsize_4;
     memset(read_index,-1,dsize_4);
/*   read the phusion output file         */
     n_contig=0;
     ctg_head[n_contig] = 0;
     while(!feof(namef))
     {
       int nPair=0;
       char line2[200],line[200],base[200];
      
       fgets(line,200,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       if((strncmp(line,"readnames",9))==0)
       { 
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==5)
            {
              memset(base,'\0',200);
              strcat(base,ptr);
              c_reads = atoi(base);
            }
         }
         ctg_list[n_contig] = c_reads;
         if(n_contig>0)
          ctg_head[n_contig] = ctg_head[n_contig-1]+ctg_list[n_contig-1];
         n_contig++;
         i_reads = 0;
       }
       else
       {      
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==0)
            {
              int idt;
              memset(base,'\0',200);
              strcat(base,ptr);
              idt = atoi(base);
              read_index[idt] = ctg_head[n_contig-1]+i_reads;
              i_reads++;
            }
         }
       }
     }

     ctg_head[n_contig] = n_reads;
     printf("matrix: %d \n",n_contig);
     tRead = n_reads;
     if(tRead%(i_edge*1000000)==0)
       nMat = tRead/(i_edge*1000000);
     else
       nMat = tRead/(i_edge*1000000) + 1;
     nRead = i_edge;
     nRead = nRead*1000000 + 1000000;

     if(tRead < nRead)
       nRead = tRead;
     rdname = cmatrix(0,nRead,0,100);
     rdbase = cmatrix(0,nRead,0,r_len);
     rdqual = cmatrix(0,nRead,0,r_len);
     if((read_length = (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - reads_length\n");
       exit(1);
     }

     printf("matrix: %d %d %d\n",n_contig,nRead,nMat);
     printf("matrix done!\n");
     if(i_end>nMat)
       i_end = nMat;
     nRead = i_edge;
     nRead = nRead*1000000;

     iMat = 0;
     map_sort[-1] = 0;
     map_sort[0] = 0;

     for(i=0;i<n_contig;i++)
     {
        if(((ctg_head[i]-ctg_head[map_sort[iMat]]) > nRead)&&((ctg_head[i-1]-ctg_head[map_sort[iMat]]) <= nRead))
	{
	  iMat++;
	  map_sort[iMat] = i-1;
          printf("map: %ld %d %d %d %d\n",i,iMat,ctg_head[i],ctg_head[i-1],map_sort[iMat]);
	}
     }

     map_sort[nMat] = n_contig;
     for(i=0;i<=nMat;i++)
        printf("sort: %ld %d %d %d\n",i,map_sort[i],ctg_head[map_sort[i]],i_end);
     iBulk = 0;
     pat_sort[iBulk] = 0;

     for(iMat=0;iMat<i_end;iMat++)
     {
	map_sort2[iMat] = iBulk;
        for(i=map_sort[iMat];i<map_sort[iMat+1];i++)
	{
	   stopflag=0;
	   j = i+1;
	   while((j<map_sort[iMat+1])&&(stopflag==0))
	   {
	     if((ctg_head[j]-ctg_head[i])<=N_SET)
	     {
//          printf("www: %d %d %d\n",j,i,ctg_head[j]-ctg_head[i]);
	       j++;
	     }
	     else
	     {
	       pat_sort[iBulk] = i;
               printf("pat1: %ld %ld %d %d %d %d\n",i,j,iBulk,ctg_head[j],ctg_head[j-1],pat_sort[iBulk]);
	       iBulk++;
	       stopflag=1;
	     }
	   }
	   i=j-1;
        }
     }

     map_sort2[nMat] = iBulk;
     map_sort[nMat] = n_contig;
     pat_sort[iBulk] = n_contig; 
     for(iMat=0;iMat<=i_end;iMat++)
        printf("sort: %d %d %d %d\n",iMat,map_sort2[iMat],map_sort[iMat],ctg_head[map_sort[iMat]]);
     for(iMat=i_start;iMat<i_end;iMat++)
     {
        int mat_low;
        int mat_hig,rlen,nlen;

	mat_low = ctg_head[pat_sort[map_sort2[iMat]]];
	mat_hig = ctg_head[pat_sort[map_sort2[iMat+1]]];
        if(iMat == (nMat-1))
          mat_hig = tRead;
        
        n_reads = 0;
        for(ac=(args+1);ac<argc;ac++)
        {
	   printf("read file: %s %d %d %d\n",argv[ac],nMat,mat_low,mat_hig);

           if((fp=fopen(argv[ac],"rb"))==NULL) {
	     	printf("Error: Cannot open file %s\n",argv[ac]);
		exit(1);
	    }
	    
           fseek(fp, 0, SEEK_END);
           Size_q_pdata = ftell(fp) + 1;
           fclose(fp);
	   
           if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL) {
             	printf("Error: calloc pdata\n");
		exit(1);
	   }
	   
           num_seqque = extractFastq(argv[ac],pdata,Size_q_pdata);
           
	   if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL) {
             	printf("Error: calloc segg\n");
		exit(1);
	   }
	   
           if((seq=decodeFastq(argv[ac],&num_seqque,&totalBases,pdata,Size_q_pdata,segg))==NULL) {
             	printf("Error: no query data found.\n");
		exit(1);
	   }

	   nRead = num_seqque;
	   printf("read file: %s %d %d %d %d\n",argv[ac],iMat,mat_low,mat_hig,nRead);
           for(i=0;i<nRead;i++)
           {
              int idk = read_index[n_reads+i];

              if(read_index[n_reads+i] >= 0)
	      {
                seqp = seq + i;
//                if((idk<mat_low)||(idk>=mat_hig))
//                        continue;
                if((idk>=mat_low)&&(idk<mat_hig))
                {
                  rlen = seqp->length;
                  read_length[idk-mat_low] = seqp->length;
		  memset(rdname[idk-mat_low],'\0',100);
		  nlen = strlen(seqp->name);
                  strncpy(rdname[idk-mat_low],seqp->name,nlen);
//		  printf("name: %s %s %d %d\n",seqp->name,rdname[idk-mat_low],nlen,rlen);
                  strncpy(rdbase[idk-mat_low],seqp->data,rlen);
                  for(j=0;j<rlen;j++)
                  {
                     rdbase[idk-mat_low][j] = seqp->data[j];
		     if(seqp->qual[j]>=99)
		       seqp->qual[j] = 99;
                     rdqual[idk-mat_low][j] = seqp->qual[j];
                  }
                }
	      }
           }
	   printf("read done: %s %d %d %d %d %d\n",argv[ac],iMat,mat_low,mat_hig,nRead,n_reads);
	   n_reads = n_reads + nRead;
	   free(pdata);
	   free(segg);
        }
	mat_low = ctg_head[pat_sort[map_sort2[iMat]]];
        printf("www: %d %d %d %d\n",iMat,map_sort2[iMat],map_sort2[iMat+1],mat_low);
        for(ac=map_sort2[iMat];ac<map_sort2[iMat+1];ac++)
	{
	   int rlen2;
           memset(nameout_hash,'\0',100);
           sprintf(nameout_hash,"%04d.fastq",ac);
	    
           if ((fpHASH=fopen(nameout_hash,"w")) == NULL) {
	    	printf("Error: can't write hash to %s\n.", nameout_hash);
		exit(1);
	   }
	   
           memset(nameout_list,'\0',100);
           sprintf(nameout_list,"readname%04d",ac); 

           if ((fpLIST=fopen(nameout_list,"w")) == NULL) {
	    	printf("Error: can't write list to %s\n", nameout_list);
		exit(1);
	   }
	   
           printf("www0: %d %d %d %d %d %d\n",ac,iMat,pat_sort[ac],pat_sort[ac+1],ctg_head[pat_sort[ac]],mat_low);
	   for(i=pat_sort[ac];i<pat_sort[ac+1];i++)
	   {
	     if(ctg_list[i]>min_size)
	     {
	      if(i_split)
	        fprintf(fpLIST,"readnames for contig %ld %ld %d\n",i,i,2*ctg_list[i]);
	      else
	        fprintf(fpLIST,"readnames for contig %ld %ld %d\n",i,i,ctg_list[i]);
	      for(j=0;j<ctg_list[i];j++)
	      {
	         offset = ctg_head[i]-mat_low+j;
                 rlen = read_length[offset];
	         if(i_split)
                 {
	           fprintf(fpHASH,"@%s.F\n",rdname[offset]);
	           rlen2 = (read_length[offset]/2)-1;
                 }
	         else
                 {
	           fprintf(fpHASH,"@%s\n",rdname[offset]);
                   rlen2  = read_length[offset];
                 }
	         for(k=0;k<rlen2;k++)
	            fprintf(fpHASH,"%c",rdbase[offset][k]);
	         fprintf(fpHASH,"\n");
	         fprintf(fpHASH,"+\n");
	         putc(041,fpHASH);
	         for(k=1;k<rlen2;k++)
	            putc(rdqual[offset][k]+041,fpHASH);
	         fprintf(fpHASH,"\n");

                 if(i_split)
	         {
	           fprintf(fpHASH,"@%s.R\n",rdname[offset]);
	           for(k=(rlen2+2);k<rlen;k++)
	              fprintf(fpHASH,"%c",rdbase[offset][k]);
	           fprintf(fpHASH,"\n+\n");
	           putc(041,fpHASH);
	           for(k=(rlen2+3);k<rlen;k++)
	              putc(rdqual[offset][k]+041,fpHASH);
	           fprintf(fpHASH,"\n");
	           fprintf(fpLIST,"%s.R\n",rdname[offset]);
	         }
	         if(i_split)
	           fprintf(fpLIST,"%s.F\n",rdname[offset]);
                 else
	           fprintf(fpLIST,"%s\n",rdname[offset]);
		 memset(rdname[offset],'\0',100);
		 memset(rdbase[offset],'\0',r_len);
		 memset(rdqual[offset],'\0',r_len);
                 read_length[offset] = 0;
	      }
	     }
	   }
	   fclose(fpHASH);
	   fclose(fpLIST);
	}
     }
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Merge_Process(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i,k;
     B64_long offset,sizesm;
     B64_long tKmer,nsize[10];
     B64_long *n_head;
     char nameout_hash[100],nameout_idex[100],nameout_list[100],nameout_size[100];
     int *read_index,*ctg2_reads,*n_list; 
     int j,i_reads = 0,c_reads,tRead;
     int **r_matrix,**c_matrix,*rarray,*ctg_list,*ctg2_list,*ctg2_head;
     unsigned char *carray;
     char *ptr;
     unsigned char *rd_used,*rd_rept,*rd_group,*rd_mask;
     FILE *namef,*fpHASH,*fpIDEX,*fpLIST,*fpSIZE;
     B64_long r_size,dsize_4;
     int step_flag,NLINK,max_contig;
     int *rd_name,*rd_link;
     int *num_front,*ind_front,*ctg_index,*ctg_links;
     int n_maxblock,split_flag,n_ungroup,n_blockreads = 0,n_reads; 

     tKmer = 0;
     tRead = 0;
     memset(nameout_size,'\0',100);
     sprintf(nameout_size,"matrix.size");
     
     if ((fpSIZE = fopen(nameout_size, "rb")) == NULL) {
        printf("Error: can't open file %s.\n", nameout_size);
        exit(1);
     }

     offset = 0;
     fseek(fpSIZE,offset,0);
     sizesm=2;
     fread(nsize,sizeof(B64_long),sizesm,fpSIZE);
     tKmer = nsize[0];
     tRead = nsize[1];
     printf("size: %ld %d\n",tKmer,tRead);

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }

/*   read the phusion output file         */
     n_contig=0;
     n_reads = 0;
     c_reads = 0;
     while(!feof(namef))
     {
       int nPair=0;
       char line2[200],line[200],base[200];
      
       fgets(line,200,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       if((strncmp(line,"readnames",9))==0)
       { 
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==5)
            {
              memset(base,'\0',200);
              strcat(base,ptr);
                c_reads = atoi(base);
            }
         }
         n_reads = n_reads+c_reads;
         n_contig++;
       }
     }
     fclose(namef);

     printf("contig: %d %d\n",n_contig,n_reads);
     if((ctg_list = (int *)calloc(n_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_list\n");
       exit(1);
     }
     if((ctg2_list = (int *)calloc(n_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_list\n");
       exit(1);
     }
     if((ctg2_head = (int *)calloc(n_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg2_head\n");
       exit(1);
     }
     if((read_index = (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - reads_index\n");
       exit(1);
     }
     if((ctg2_reads = (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - reads_index\n");
       exit(1);
     }
     if((n_list= (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - n_list\n");
       exit(1);
     }
     if((n_head= (B64_long *)calloc(tRead,sizeof(B64_long))) == NULL)
     {
       printf("ssaha: calloc - patch_head\n");
       exit(1);
     }

     dsize_4 = tRead;
     dsize_4 = 4*dsize_4;
     memset(read_index,-1,dsize_4);
     memset(nameout_list,'\0',100); 
     sprintf(nameout_list,"matrix.list");
     
     if ((fpLIST = fopen(nameout_list, "rb")) == NULL) {
        printf("Error: can't open file %s.\n", nameout_list);
        exit(1);
     }

     offset = 0;
     fseek(fpLIST,offset,0);
     sizesm=tKmer;
     fread(n_list,sizeof(int),sizesm,fpLIST);
     fclose(fpLIST);

     n_head[0] = 0;
     for(i=1;i<=tRead;i++)
        n_head[i] = n_head[i-1] + n_list[i-1];
     r_size = tKmer;
     r_size = r_size + 5000;
     if((rarray= (int *)calloc(r_size,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }
     if((carray= (unsigned char *)calloc(r_size,sizeof(unsigned char))) == NULL)
     {
       printf("ssaha: calloc - cmatrix\n");
       exit(1);
     }

     memset(nameout_hash,'\0',100);
     sprintf(nameout_hash,"relation.matrix.link"); 
     
     if ((fpHASH = fopen(nameout_hash, "rb")) == NULL) {
        printf("Error: can't open file %s\n.", nameout_hash);
        exit(1);
     }

     offset = 0;
     fseek(fpHASH,offset,0);
     sizesm=tKmer;
     fread(rarray,sizeof(int),sizesm,fpHASH);
     fclose(fpHASH);

     memset(nameout_idex,'\0',100);
     sprintf(nameout_idex,"relation.matrix.read"); 
     
     if ((fpIDEX = fopen(nameout_idex, "rb")) == NULL) {
        printf("Error: can't open file %s\n.", nameout_idex);
        exit(1);
     }

     offset = 0;
     fseek(fpIDEX,offset,0);
     sizesm=tKmer;
     fread(carray,sizeof(unsigned char),sizesm,fpIDEX);
     fclose(fpIDEX);

     nRead = n_contig; 
     r_matrix = imatrix(0,nRead,0,200);
     c_matrix = imatrix(0,nRead,0,200);

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }

     n_maxblock = 0;
     dsize_4 = tRead;
     dsize_4 = 4*dsize_4;
     memset(read_index,-1,dsize_4);
/*   read the phusion output file         */
     n_contig=0;
     while(!feof(namef))
     {
       int nPair=0;
       char line2[200],line[200],base[200];
      
       fgets(line,200,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       if((strncmp(line,"readnames",9))==0)
       { 
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==5)
            {
              memset(base,'\0',200);
              strcat(base,ptr);
              c_reads = atoi(base);
            }
         }
         ctg2_list[n_contig] = c_reads;
	 if(n_contig>0)
	   ctg2_head[n_contig] = ctg2_head[n_contig-1]+ctg2_list[n_contig-1];
	 if(c_reads>n_maxblock)
	   n_maxblock = c_reads;
         n_contig++;
         i_reads = 0;
       }
       else
       {
	 nPair=0;
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==0)
            {
              int idt;
              memset(base,'\0',200);
              strcat(base,ptr);
              idt = atoi(base);
	      ctg2_reads[ctg2_head[n_contig-1]+i_reads] = idt;
              read_index[idt] = n_contig-1;
              i_reads++;
            }
         }
       }
     }

     if((ctg_index= (int *)calloc(n_maxblock,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - ctg_index\n");
       exit(1);
     }
     if((ctg_links= (int *)calloc(n_maxblock,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - ctg_links\n");
       exit(1);
     }
     printf("matrix: %d \n",n_contig);
     for(i=0;i<n_contig;i++)
        ctg_list[i] = 0;

     for(i=0;i<tRead;i++)
     {
        if(read_index[i] >= 0)
        {
          int i_contig = read_index[i];
          int cidk = ctg_list[i_contig];

	  if(cidk >= (max_edge-2)) continue;
	  if(ctg2_list[i_contig] >= n_group) continue;
//	  if(ctg2_list[i_contig] >= 10000) continue;
          for(j=0;j<n_list[i];j++)
          {
             B64_long k_offset = n_head[i]+j;
             int a_index = rarray[k_offset];
             int c_index = carray[k_offset];
             int j_contig = read_index[a_index];
             if((i_contig!=j_contig)&&(c_index>=10)&&(j_contig)&&(j_contig >= 0))
//             if((i_contig!=j_contig)&&(c_index>=20)&&(j_contig)&&(j_contig >= 0))
             {
               int idi = read_index[a_index];
               int match = 0;
               int km;
               for(km=0;km<cidk;km++)
               {
                  if(r_matrix[i_contig][km]==idi)
                  {
                    c_matrix[i_contig][km]++;
                    match = 1;
                    break;
                  } 
               }
               if(match==0)
               {
                 r_matrix[i_contig][cidk]=idi;
                 c_matrix[i_contig][cidk]=1;
                 cidk++;
                 ctg_list[i_contig]++;
               }
             }
          }   
        }
     }

/*     for(i=0;i<n_contig;i++)
     {
        printf("contigs: %d %d %d ",i,ctg_list[i],ctg2_list[i]);
        for(j=0;j<ctg_list[i];j++)
        {
           printf("| %d %d %d ",r_matrix[i][j],c_matrix[i][j],ctg2_list[r_matrix[i][j]]); 
        }
        printf("\n");
     }   */
     nRead = n_contig;
     if((ind_front= (int *)calloc(nRead,sizeof(int))) == NULL)//worst case all reads are linked together
     {
       printf("ssaha: calloc - ind_front\n");
       exit(1);
     }
     if((num_front= (int *)calloc(nRead,sizeof(int))) == NULL)//worst case all reads are linked together
     {
       printf("ssaha: calloc - num_front\n");
       exit(1);
     }
     if((rd_name= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - rd_name\n");
       exit(1);
     }
     if((rd_link= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - rd_link\n");
       exit(1);
     }
     if((rd_used= (unsigned char *)calloc(nRead,sizeof(unsigned char))) == NULL)
     {
       printf("ssaha: calloc - rd_used\n");
       exit(1);
     }
     if((rd_rept= (unsigned char *)calloc(nRead,sizeof(unsigned char))) == NULL)
     {
       printf("ssaha: calloc - rd_rept\n");
       exit(1);
     }
     if((rd_group= (unsigned char *)calloc(nRead,sizeof(unsigned char))) == NULL)
     {
       printf("phusion: calloc - rd_group\n");
       exit(1);
     }
     if((rd_mask= (unsigned char *)calloc(nRead,sizeof(unsigned char))) == NULL)
     {
       printf("phusion: calloc - rd_mask\n");
       exit(1);
     }

     if((namef = fopen(argv[args+1],"w")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }
     n_maxblock=nRead;
     split_flag=1;
     step_flag=1;
     NLINK=10;
     n_contig=0;
     mat_ed=500;
     min_size = 1;
     printf("links: %d %d\n",NLINK,nRead);
     n_ungroup = 0; 
   while((n_maxblock>0*N_SET)&&split_flag&&step_flag)
   {
     printf("Number of reads: %d %d %d %d %d\n",nRead,n_ungroup,n_blockreads,mat_ed,n_maxblock);
     n_maxblock=0;
     n_blockreads=0;

/*   get the readnames for each contigs   */
     printf("match loop finished ....\n");
     fflush(stdout);
     system("ps aux | grep phusion; date");

     max_contig=0;
     dsize_4 = nRead;
     dsize_4 = 4*nRead;
     memset(rd_used,0,nRead);
     memset(rd_rept,0,nRead);
     memset(rd_name,0,dsize_4);
     memset(ind_front,0,dsize_4);
     memset(num_front,0,dsize_4);
     n_ungroup=0;
     for(i=0;i<nRead;i++)
     {
        if(rd_group[i]==0)
          n_ungroup++;
        rd_used[i]=rd_group[i];
     }
     n_reads=0;
     n_block=0;
     for(i=0;i<nRead;i++)
     {
        int id_pt,n_front,n_midd,n_end,n_over;
        int a_index,a_links,fit2,flag_link=0,i_links=0;

        n_reads=0;
        n_over=0;
        n_front=ctg_list[i];
        if((rd_used[i]==0)&&(n_front>0))
        {
          for(j=0;j<n_front;j++)
          {
             a_index=r_matrix[i][j];
             a_links=c_matrix[i][j];
             ind_front[j]=a_index;
             num_front[j]=a_links;
             if(a_links>=NLINK)
             {
               if(flag_link==0)
               {
                 flag_link=1;
                 i_links=a_links;
               }
               n_over++;
             }
          }
          if(n_over==0) 
            n_front=0;
          else
          {
            rd_used[i]=1;
            rd_name[n_reads]=i;
            n_reads++;
            rd_rept[i]=1;
            rd_link[i]=i_links;
            for(j=0;j<n_front;j++)
            {
               a_index=r_matrix[i][j];
               a_links=c_matrix[i][j];
               if((a_links>=NLINK)&&(rd_used[a_index]==0))
               {
                 rd_name[n_reads]=a_index;
                 rd_rept[a_index]=1;
                 rd_link[a_index]=a_links;
                 n_reads++;
               }
            }
          }
          while(n_front>0)
          {
            n_midd=0;
            for(j=0;j<n_front;j++)
            {
               id_pt=ind_front[j];
               if((num_front[j]>=NLINK)&&(rd_used[id_pt]==0))
               {
                 rd_used[id_pt]=1;
                 for(k=0;k<ctg_list[id_pt];k++)
                 {
                    a_index=r_matrix[id_pt][k];
                    a_links=c_matrix[id_pt][k];
                    if((a_links>=NLINK)&&(rd_rept[a_index]==0)&&(rd_used[a_index]==0))
                    {
                      rd_name[n_reads]=a_index;
                      n_reads++;
                      rd_rept[a_index]=1;
                      rd_link[a_index]=a_links;
                      num_front[n_front+n_midd]=a_links;
                      ind_front[n_front+n_midd]=a_index;
                      n_midd++;
                    }
                 }
               }
            }
            n_end=0;
            for(j=0;j<(n_front+n_midd);j++)
            {
               if((rd_used[ind_front[j]]==0)&&(num_front[j]>=NLINK))
               {
                 num_front[n_end]=num_front[j];
                 ind_front[n_end]=ind_front[j];
                 n_end++; 
               }
            }
            n_front=n_end;
          }
          if(n_reads>=min_size)
          {
//            fit1=(n_reads>N_SET)&&(NLINK!=(mat_ed-s_len));
            fit2=(n_reads>N_SET)&&(NLINK<(mat_ed-s_len));
            if(fit2)
            {
              for(j=0;j<n_reads;j++)
              {
                 rd_group[rd_name[j]]=0;
              }
              n_blockreads=n_blockreads+n_reads;
              n_block++;
              if(n_reads>=n_maxblock)
                n_maxblock=n_reads;
            }
            else
            {
	      int *rd_index,n2_reads =0;
              int num_gp_reads = 0;

	      if((rd_index= (int *)calloc(n_reads,sizeof(int))) == NULL)
	      {
	        printf("ssaha: calloc - rd_index\n");
	        exit(1);
	      }
              for(j=0;j<n_reads;j++)
	      {
	         rd_index[j] = rd_name[j];
		 rd_mask[rd_name[j]]=0;
	      }
              for(j=0;j<n_reads;j++)
              {
	         if(rd_mask[rd_name[j]]==0)
		 {
		   rd_mask[rd_name[j]] = 1;
		   rd_index[n2_reads] = rd_name[j];
		   num_gp_reads = num_gp_reads + ctg2_list[rd_name[j]];
		   n2_reads++;
		 }
              }
	      if(n2_reads>=min_size)
	      {
                printf("readnames for contig %d read_cluster %d %d\n",n_contig,n2_reads,num_gp_reads);
                fprintf(namef,"readnames for contig %d read_cluster %d\n",n_contig,num_gp_reads);
                for(j=0;j<n2_reads;j++)
                {
                   printf("%d read %d %d\n",rd_index[j],rd_link[rd_index[j]],ctg2_list[rd_index[j]]);
	           for(k=0;k<ctg2_list[rd_index[j]];k++)
                      fprintf(namef,"%d read %d\n",ctg2_reads[ctg2_head[rd_index[j]]+k],10);
		   rd_used[rd_index[j]] = 1;
                   rd_group[rd_index[j]]=1;
/*                 printf("%d read %d\n",rd_name[j],rd_link[rd_name[j]]);
		   rd_used[rd_name[j]] = 1;
                   rd_group[rd_name[j]]=1;      */
                }
                n_contig++;
	      }
	      free(rd_index);
            }
            if(n_reads>=max_contig)
              max_contig=n_reads;
          }
        }
     }
     printf("phase %d finished ....\n",pmod);
     if(max_contig<N_SET)
       step_flag=0;
     NLINK=NLINK+s_len;
     printf("links: %d %d\n",NLINK,max_contig); 
     mat_st=nmatch;
     mat_ed=500;
     if(NLINK>=mat_ed)
     {
       step_flag=0;
       NLINK=NLINK-s_len;
     }
   }
   for(i=0;i<nRead;i++)
   {
      if(rd_group[i]==0)
      {
        printf("readnames for contig %d read_cluster %d %d\n",n_contig,1,ctg2_list[i]);
        printf("%ld read %d %d\n",i,rd_link[i],ctg2_list[i]);
        fprintf(namef,"readnames for contig %d read_cluster %d\n",n_contig,ctg2_list[i]);
	for(j=0;j<ctg2_list[i];j++)
           fprintf(namef,"%d read %d\n",ctg2_reads[ctg2_head[i]+j],10);
	n_contig++;
      }
   }
   fclose(namef);
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Remap_Process(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i;
     B64_long offset,sizesm,dsize_4;
     B64_long tKmer,nsize[10];
     B64_long *n_head;
     char nameout_hash[100],nameout_idex[100],nameout_list[100],nameout_size[100];
     int *read_index,*ctg2_reads,*n_list; 
     int j,i_reads = 0,c_reads,tRead;
     int **r_matrix,**c_matrix,*rarray,*ctg_list,*ctg_head,*ctg2_list,*ctg2_head;
     int **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     unsigned char *carray;
     char *ptr;
     FILE *namef,*fp,*fpHASH,*fpIDEX,*fpLIST,*fpSIZE;
     B64_long r_size,n_remap;
     int n_reads,*reads_remap,*ctg_maps,*map_index,*map_link,print_flag; 

     tKmer = 0;
     tRead = 0;
     memset(nameout_size,'\0',100);
     sprintf(nameout_size,"matrix.size");
     
     if ((fpSIZE = fopen(nameout_size, "rb")) == NULL) {
        printf("Error: can't open file %s.\n", nameout_size);
        exit(1);
     }

     offset = 0;
     fseek(fpSIZE,offset,0);
     sizesm=2;
     fread(nsize,sizeof(B64_long),sizesm,fpSIZE);
     tKmer = nsize[0];
     tRead = nsize[1];
     printf("size: %ld %d\n",tKmer,tRead);

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }

/*   read the phusion output file         */
     n_contig=0;
     n_reads = 0;
     while(!feof(namef))
     {
       int nPair=0;
       char line2[200],line[200],base[200];
      
       fgets(line,200,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       if((strncmp(line,"readnames",9))==0)
       { 
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==5)
            {
              memset(base,'\0',200);
              strcat(base,ptr);
                c_reads = atoi(base);
            }
         }
         n_reads = n_reads+c_reads;
         n_contig++;
       }
     }
     fclose(namef);

     printf("contig: %d %d %s\n",n_contig,n_reads,argv[args]);
     if((ctg_list = (int *)calloc(n_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_list\n");
       exit(1);
     }
     if((ctg_head = (int *)calloc(n_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_head\n");
       exit(1);
     }
     if((ctg_maps = (int *)calloc(n_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_maps\n");
       exit(1);
     }
     if((map_index = (int *)calloc(max_edge,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - map_index\n");
       exit(1);
     }
     if((map_link = (int *)calloc(max_edge,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - map_link\n");
       exit(1);
     }
     if((ctg2_list = (int *)calloc(n_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg2_list\n");
       exit(1);
     }
     if((ctg2_head = (int *)calloc(n_contig,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg2_head\n");
       exit(1);
     }
     if((read_index = (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - reads_index\n");
       exit(1);
     }
     if((ctg2_reads = (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - reads_index\n");
       exit(1);
     }
     if((n_list= (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - n_list\n");
       exit(1);
     }
     if((n_head= (B64_long *)calloc(tRead,sizeof(B64_long))) == NULL)
     {
       printf("ssaha: calloc - n_head\n");
       exit(1);
     }

     r_size = tRead;
     r_size = 4*r_size;
     memset(read_index,-1,r_size);
     memset(nameout_list,'\0',100); 
     sprintf(nameout_list,"matrix.list");
     
     if ((fpLIST = fopen(nameout_list, "rb")) == NULL) {
        printf("Error: can't open file %s.\n", nameout_list);
        exit(1);
     }

     offset = 0;
     fseek(fpLIST,offset,0);
     sizesm=tKmer;
     fread(n_list,sizeof(int),sizesm,fpLIST);
     fclose(fpLIST);

     n_head[0] = 0;
     for(i=1;i<=tRead;i++)
        n_head[i] = n_head[i-1] + n_list[i-1];
     r_size = tKmer;
     r_size = r_size + 5000;
     if((rarray= (int *)calloc(r_size,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }
     if((carray= (unsigned char *)calloc(r_size,sizeof(unsigned char))) == NULL)
     {
       printf("ssaha: calloc - cmatrix\n");
       exit(1);
     }

     memset(nameout_hash,'\0',100);
     sprintf(nameout_hash,"relation.matrix.link"); 
     
     if ((fpHASH = fopen(nameout_hash, "rb")) == NULL) {
        printf("Error: can't open file %s\n.", nameout_hash);
        exit(1);
     }

     offset = 0;
     fseek(fpHASH,offset,0);
     sizesm=tKmer;
     fread(rarray,sizeof(int),sizesm,fpHASH);
     fclose(fpHASH);

     memset(nameout_idex,'\0',100);
     sprintf(nameout_idex,"relation.matrix.read"); 
     
     if ((fpIDEX = fopen(nameout_idex, "rb")) == NULL) {
        printf("Error: can't open file %s\n.", nameout_idex);
        exit(1);
     }

     offset = 0;
     fseek(fpIDEX,offset,0);
     sizesm=tKmer;
     fread(carray,sizeof(unsigned char),sizesm,fpIDEX);
     fclose(fpIDEX);

     nRead = n_contig; 
     r_matrix = imatrix(0,nRead,0,200);
     c_matrix = imatrix(0,nRead,0,200);

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }

     dsize_4 = tRead;
     dsize_4 = 4*dsize_4;
     memset(read_index,-1,dsize_4);
/*   read the phusion output file         */
     n_contig=0;
     while(!feof(namef))
     {
       int nPair=0;
       char line2[200],line[200],line3[200],base[200];
      
       fgets(line,200,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       strcpy(line3,line);
       if((strncmp(line,"readnames",9))==0)
       { 
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==5)
            {
              memset(base,'\0',200);
              strcat(base,ptr);
              c_reads = atoi(base);
            }
         }
         ctg2_list[n_contig] = c_reads;
	 if(n_contig>0)
	   ctg2_head[n_contig] = ctg2_head[n_contig-1]+ctg2_list[n_contig-1];
         n_contig++;
         i_reads = 0;
       }
       else
       {      
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==0)
            {
              int idt;
              memset(base,'\0',200);
              strcat(base,ptr);
              idt = atoi(base);
	      ctg2_reads[ctg2_head[n_contig-1]+i_reads] = idt;
              read_index[idt] = n_contig-1;
              i_reads++;
            }
         }
       }
     }
     fclose(namef);

     n_remap = 0;
     for(i=0;i<tRead;i++)
     {
        if(read_index[i] < 0)
        {
          int remap_flag = 0;
          int i_remap = 0;
          int m_remap = 0;
          int m_index = 0;
	  if(n_list[i] >=  max_edge)
	    n_list[i] = max_edge-1;
//         printf("index: %d %d %d ",i,n_list[i],i_reads);
          for(j=0;j<n_list[i];j++)
          {
             B64_long k_offset = n_head[i]+j;
             int a_index = rarray[k_offset];
             if(read_index[a_index]>=0)
             {
               map_link[i_remap] = read_index[a_index];
               ctg_maps[read_index[a_index]]++;
               i_remap++;
               remap_flag = 1;
             }
//           printf("| %d %d %d %d ",j,a_index,c_index,read_index[a_index]);
          }
//        printf("\n");
          for(j=0;j<i_remap;j++)
          {
             if(ctg_maps[map_link[j]]>m_remap)
             {
               m_remap = ctg_maps[map_link[j]];
               m_index = map_link[j];
             }
          }
          if(remap_flag)
          {
//         printf("index-xxx: %d %d %d\n",i,m_index,ctg_list[m_index]);
            ctg_list[m_index]++;
            n_remap++;
          }
          for(j=0;j<i_remap;j++)
             ctg_maps[map_link[j]] = 0;
        }
     }
     printf("remap: %ld\n",n_remap);

     n_remap = n_remap + 10;
     if((reads_remap= (int *)calloc(n_remap,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - n_list\n");
       exit(1);
     }
     printf("remap2: %ld\n",n_remap);

     ctg_head[0] = 0;
     ctg_maps[0] = 0;
     for(i=1;i<n_contig;i++)
     {
        ctg_head[i] = ctg_head[i-1] + ctg_list[i-1];
	ctg_list[i-1] = 0;
	ctg_maps[i] = 0;
     }
     ctg_list[n_contig-1] = 0;
     for(i=0;i<tRead;i++)
     {
        if(read_index[i] < 0)
        {
          int remap_flag = 0;
          int i_remap = 0;
          int m_remap = 0;
          int m_index = 0;
          for(j=0;j<n_list[i];j++)
          {
             B64_long k_offset = n_head[i]+j;
             int a_index = rarray[k_offset];
             if(read_index[a_index]>=0)
             {
               map_link[i_remap] = read_index[a_index];
               ctg_maps[read_index[a_index]]++;
               i_remap++;
               remap_flag = 1;
             }
          }
          for(j=0;j<i_remap;j++)
          {
             if(ctg_maps[map_link[j]]>m_remap)
             {
               m_remap = ctg_maps[map_link[j]];
               m_index = map_link[j];
             }
          }
          if(remap_flag)
          {
            int offset = ctg_head[m_index]+ctg_list[m_index];
            reads_remap[offset] = i; 
//          printf("www: %d %d %d %d %d\n",i,n_contig,m_index,ctg_list[m_index],ctg2_list[m_index]);
            ctg_list[m_index]++;
            n_remap++;
          }
          for(j=0;j<i_remap;j++)
             ctg_maps[map_link[j]] = 0;
        }
     }

/*     
     for(i=0;i<n_contig;i++)
     {
        if(ctg_list[i]>0)
        {
          printf("contigs: %d %d %d ",i,ctg_list[i],ctg2_list[i]);
          for(j=0;j<ctg_list[i];j++)
          {
             int offset = ctg_head[i]+j;
             printf("| %d %d ",j,reads_remap[offset]);
          }
          printf("\n");
        }
     }  */
     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }
     if((fp = fopen(argv[args+1],"w")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }

/*   read the phusion output file         */
     n_contig=0;
     n_reads = 0;
     print_flag = 0;
     while(!feof(namef))
     {
       int nPair=0;
       char line2[200],line[200],base[200];
      
       fgets(line,200,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       if((strncmp(line,"readnames",9))==0)
       { 
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==5)
            {
              memset(base,'\0',200);
              strcat(base,ptr);
                c_reads = atoi(base);
            }
         }
	 if(print_flag)
	 {
	   for(j=0;j<ctg_list[n_contig-1];j++)
	   {
              int offset = ctg_head[n_contig-1]+j;
              fprintf(fp,"%d read 1\n",reads_remap[offset]);
	   }
	 }
	 fprintf(fp,"readnames for contig %d read_cluster %d\n",n_contig,c_reads+ctg_list[n_contig]);
         n_reads = n_reads+c_reads;
         n_contig++;
	 print_flag = 0;
       }
       else
       {      
	 fprintf(fp,"%s",line);
	 if(ctg_list[n_contig-1]>0)
	   print_flag = 1;
       }
     }
     fclose(namef);
     if(print_flag)
     {
       for(j=0;j<ctg_list[n_contig-1];j++)
       {
          int offset = ctg_head[n_contig-1]+j;
          fprintf(fp,"%d read 1\n",reads_remap[offset]);
       }
     }
     fclose(fp);

}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Clust_Process(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i,k;
     B64_long offset,sizesm;
     B64_long tKmer,nsize[10];
     B64_long *n_head;
     char nameout_hash[100],nameout_idex[100],nameout_list[100],nameout_size[100];
     int *n_list; 
     int j,tRead;
     int *rarray;
     unsigned char *carray;
     unsigned char *rd_used,*rd_rept,*rd_group,*rd_mask;
     FILE *fpHASH,*fpIDEX,*fpLIST,*fpSIZE;
     B64_long r_size,dsize_4;
     int step_flag,NLINK,max_contig;
     int *rd_name,*rd_link;
     int *num_front,*ind_front;
     int n_maxblock,split_flag,n_ungroup,n_blockreads = 0,n_reads; 

     tKmer = 0;
     tRead = 0;
     memset(nameout_size,'\0',100);
     sprintf(nameout_size,"matrix.size");
     
     if ((fpSIZE = fopen(nameout_size, "rb")) == NULL) {
        printf("Error: can't open file %s.\n", nameout_size);
        exit(1);
     }

     offset = 0;
     fseek(fpSIZE,offset,0);
     sizesm=2;
     fread(nsize,sizeof(B64_long),sizesm,fpSIZE);
     tKmer = nsize[0];
     tRead = nsize[1];
     printf("size: %ld %d\n",tKmer,tRead);

     if((n_list= (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - n_list\n");
       exit(1);
     }
     if((n_head= (B64_long *)calloc(tRead,sizeof(B64_long))) == NULL)
     {
       printf("ssaha: calloc - patch_head\n");
       exit(1);
     }
     memset(nameout_list,'\0',100); 
     sprintf(nameout_list,"matrix.list");
     
     if ((fpLIST = fopen(nameout_list, "rb")) == NULL) {
        printf("Error: can't open file %s.\n", nameout_list);
        exit(1);
     }

     offset = 0;
     fseek(fpLIST,offset,0);
     sizesm=tKmer;
     fread(n_list,sizeof(int),sizesm,fpLIST);
     fclose(fpLIST);

     n_head[0] = 0;
     for(i=1;i<=tRead;i++)
        n_head[i] = n_head[i-1] + n_list[i-1];
     r_size = tKmer;
     r_size = r_size + 5000;
     if((rarray= (int *)calloc(r_size,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }
     if((carray= (unsigned char *)calloc(r_size,sizeof(unsigned char))) == NULL)
     {
       printf("ssaha: calloc - cmatrix\n");
       exit(1);
     }

     memset(nameout_hash,'\0',100);
     sprintf(nameout_hash,"relation.matrix.link"); 
     
     if ((fpHASH = fopen(nameout_hash, "rb")) == NULL) {
        printf("Error: can't open file %s\n.", nameout_hash);
        exit(1);
     }

     offset = 0;
     fseek(fpHASH,offset,0);
     sizesm=tKmer;
     fread(rarray,sizeof(int),sizesm,fpHASH);
     fclose(fpHASH);

     memset(nameout_idex,'\0',100);
     sprintf(nameout_idex,"relation.matrix.read"); 
     
     if ((fpIDEX = fopen(nameout_idex, "rb")) == NULL) {
        printf("Error: can't open file %s\n.", nameout_idex);
        exit(1);
     }

     offset = 0;
     fseek(fpIDEX,offset,0);
     sizesm=tKmer;
     fread(carray,sizeof(unsigned char),sizesm,fpIDEX);
     fclose(fpIDEX);

     nRead = tRead;
     if((ind_front= (int *)calloc(nRead,sizeof(int))) == NULL)//worst case all reads are linked together
     {
       printf("ssaha: calloc - ind_front\n");
       exit(1);
     }
     if((num_front= (int *)calloc(nRead,sizeof(int))) == NULL)//worst case all reads are linked together
     {
       printf("ssaha: calloc - num_front\n");
       exit(1);
     }
     if((rd_name= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - rd_name\n");
       exit(1);
     }
     if((rd_link= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - rd_link\n");
       exit(1);
     }
     if((rd_used= (unsigned char *)calloc(nRead,sizeof(unsigned char))) == NULL)
     {
       printf("ssaha: calloc - rd_used\n");
       exit(1);
     }
     if((rd_rept= (unsigned char *)calloc(nRead,sizeof(unsigned char))) == NULL)
     {
       printf("ssaha: calloc - rd_rept\n");
       exit(1);
     }
     if((rd_group= (unsigned char *)calloc(nRead,sizeof(unsigned char))) == NULL)
     {
       printf("phusion: calloc - rd_group\n");
       exit(1);
     }
     if((rd_mask= (unsigned char *)calloc(nRead,sizeof(unsigned char))) == NULL)
     {
       printf("phusion: calloc - rd_mask\n");
       exit(1);
     }

     n_maxblock=nRead;
     split_flag=1;
     step_flag=1;
     NLINK=mat_st;
     n_contig=0;
     mat_ed=50;
     printf("links: %d %d\n",NLINK,nRead);
     n_ungroup = 0; 
   while((n_maxblock>0*N_SET)&&split_flag&&step_flag)
   {
     printf("Number of reads: %d %d %d %d %d\n",nRead,n_ungroup,n_blockreads,mat_ed,n_maxblock);
     n_maxblock=0;
     n_blockreads=0;

/*   get the readnames for each contigs   */
     printf("match loop finished ....\n");
     fflush(stdout);
     system("ps aux | grep phusion; date");

     max_contig=0;
     dsize_4 = nRead;
     dsize_4 = 4*dsize_4;
     memset(rd_used,0,nRead);
     memset(rd_rept,0,nRead);
     memset(rd_name,0,dsize_4);
     memset(ind_front,0,dsize_4);
     memset(num_front,0,dsize_4);
     n_ungroup=0;
     for(i=0;i<nRead;i++)
     {
        if(rd_group[i]==0)
          n_ungroup++;
        rd_used[i]=rd_group[i];
     }
     n_reads=0;
     n_block=0;
     for(i=0;i<nRead;i++)
     {
        B64_long k_offset;
        int id_pt,n_front,n_midd,n_end,n_over;
        int a_index,a_links,fit2,flag_link=0,i_links=0;

        n_reads=0;
        n_over=0;
        n_front=n_list[i];
        if((rd_used[i]==0)&&(n_front>0))
        {
	  k_offset = n_head[i];
          for(j=0;j<n_front;j++)
          {
             a_index=rarray[k_offset+j];
             a_links=carray[k_offset+j];
             ind_front[j]=a_index;
             num_front[j]=a_links;
//	printf("www: %d %d %d %d %d %d\n",i,rd_used[a_index],n_front,j,ind_front[j],num_front[j]);
             if(a_links>=NLINK)
             {
               if(flag_link==0)
               {
                 flag_link=1;
                 i_links=a_links;
               }
               n_over++;
             }
          }
          if(n_over==0) 
            n_front=0;
          else
          {
            rd_used[i]=1;
            rd_name[n_reads]=i;
            n_reads++;
            rd_rept[i]=1;
            rd_link[i]=i_links;
	    k_offset = n_head[i];
            for(j=0;j<n_front;j++)
            {
               a_index=rarray[k_offset+j];
               a_links=carray[k_offset+j];
               if((a_links>=NLINK)&&(rd_used[a_index]==0))
               {
                 rd_name[n_reads]=a_index;
                 rd_rept[a_index]=1;
                 rd_link[a_index]=a_links;
                 n_reads++;
               }
            }
          }
          while(n_front>0)
          {
            n_midd=0;
            for(j=0;j<n_front;j++)
            {
               id_pt=ind_front[j];
               if((num_front[j]>=NLINK)&&(rd_used[id_pt]==0))
               {
                 rd_used[id_pt]=1;
		 k_offset = n_head[id_pt];
                 for(k=0;k<n_list[id_pt];k++)
                 {
                    a_index=rarray[k_offset+k];
                    a_links=carray[k_offset+k];
//                printf("www1: %d %d %d %d %d %d\n",id_pt,rd_used[a_index],n_list[id_pt],k,a_index,a_links);	
                    if((a_links>=NLINK)&&(rd_rept[a_index]==0)&&(rd_used[a_index]==0))
                    {
                      rd_name[n_reads]=a_index;
                      n_reads++;
                      rd_rept[a_index]=1;
                      rd_link[a_index]=a_links;
                      num_front[n_front+n_midd]=a_links;
                      ind_front[n_front+n_midd]=a_index;
                      n_midd++;
                    }
                 }
               }
            }
            n_end=0;
            for(j=0;j<(n_front+n_midd);j++)
            {
               if((rd_used[ind_front[j]]==0)&&(num_front[j]>=NLINK))
               {
                 num_front[n_end]=num_front[j];
                 ind_front[n_end]=ind_front[j];
                 n_end++; 
               }
            }
            n_front=n_end;
          }
          if(n_reads>=min_size)
          {
//            fit1=(n_reads>N_SET)&&(NLINK!=(mat_ed-s_len));
            fit2=(n_reads>N_SET)&&(NLINK<(mat_ed-s_len));
            if(fit2)
            {
              for(j=0;j<n_reads;j++)
              {
                 rd_group[rd_name[j]]=0;
              }
              n_blockreads=n_blockreads+n_reads;
//              printf("contig: %d %d %d %d\n",n_contig,n_reads,n_block,n_blockreads);
              n_block++;
              if(n_reads>=n_maxblock)
                n_maxblock=n_reads;
            }
            else
            {
	      int *rd_index,n2_reads =0;

	      if((rd_index= (int *)calloc(n_reads,sizeof(int))) == NULL)
	      {
	        printf("ssaha: calloc - rd_index\n");
	        exit(1);
	      }
              for(j=0;j<n_reads;j++)
	      {
	         rd_index[j] = rd_name[j];
		 rd_mask[rd_name[j]]=0;
	      }
              for(j=0;j<n_reads;j++)
              {
	         if(rd_mask[rd_name[j]]==0)
		 {
		   rd_mask[rd_name[j]] = 1;
		   rd_index[n2_reads] = rd_name[j];
		   n2_reads++;
		 }
              }
	      if(n2_reads>=min_size)
	      {
                printf("readnames for contig %d read_cluster %d\n",n_contig,n2_reads);
                for(j=0;j<n2_reads;j++)
                {
                   printf("%d read %d\n",rd_index[j],rd_link[rd_index[j]]);
		   rd_used[rd_index[j]] = 1;
                   rd_group[rd_index[j]]=1;
/*                 printf("%d read %d\n",rd_name[j],rd_link[rd_name[j]]);
		   rd_used[rd_name[j]] = 1;
                   rd_group[rd_name[j]]=1;      */
                }
                n_contig++;
	      }
	      free(rd_index);
            }
            if(n_reads>=max_contig)
              max_contig=n_reads;
          }
        }
     }
     printf("phase %d finished ....\n",pmod);
     if(max_contig<N_SET)
       step_flag=0;
     NLINK=NLINK+s_len;
     printf("links: %d %d\n",NLINK,max_contig); 
     mat_st=nmatch;
     mat_ed=50;
     if(NLINK>=mat_ed)
     {
       step_flag=0;
       NLINK=NLINK-s_len;
     }
   }
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Matrix_Process(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i;
     B64_long offset,sizesm;
     B64_long tKmer,nsize[10];
     B64_long *n_head;
     int **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     unsigned char **cmatrix(long nrl,long nrh,long ncl,long nch);
     char nameout_hash[100],nameout_idex[100],nameout_list[100],nameout_size[100];
     int *n_list; 
     int j,tRead,iMat,nMat;
     int *r_matrix,*rarray;
     unsigned char *c_matrix,*carray;
     FILE *fpHASH,*fpIDEX,*fpLIST,*fpSIZE;
     B64_long l_size,r_size,rsize,rsize2;

     tKmer = 0;
     tRead = 0;
     memset(nameout_size,'\0',100);
     sprintf(nameout_size,"matrix.size");
     
     if ((fpSIZE = fopen(nameout_size, "rb")) == NULL) {
        printf("Error: can't open file %s.\n", nameout_size);
        exit(1);
     }

     offset = 0;
     fseek(fpSIZE,offset,0);
     sizesm=2;
     fread(nsize,sizeof(B64_long),sizesm,fpSIZE);
     tKmer = nsize[0];
     tRead = nsize[1];
     printf("size: %ld %d\n",tKmer,tRead);

     if((n_list= (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - n_list\n");
       exit(1);
     }
     if((n_head= (B64_long *)calloc((tRead+1),sizeof(B64_long))) == NULL)
     {
       printf("ssaha: calloc - patch_head\n");
       exit(1);
     }
     memset(nameout_list,'\0',100); 
     sprintf(nameout_list,"matrix.list");
     
     if ((fpLIST = fopen(nameout_list, "rb")) == NULL) {
        printf("Error: can't open file %s.\n", nameout_list);
        exit(1);
     }

     offset = 0;
     fseek(fpLIST,offset,0);
     sizesm=tKmer;
     fread(n_list,sizeof(int),sizesm,fpLIST);
     fclose(fpLIST);

     n_head[0] = 0;
     for(i=1;i<=tRead;i++)
        n_head[i] = n_head[i-1] + n_list[i-1];
     l_size = tKmer;
     l_size = l_size + 5000;
     if((rarray= (int *)calloc(l_size,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }
     if(tRead%(i_edge*1000000)==0)
       nMat = tRead/(i_edge*1000000);
     else
       nMat = tRead/(i_edge*1000000) + 1;
     nRead = i_edge;
     nRead = nRead*1000000 + 5000;
     r_size = nRead;
     r_size = r_size*max_edge;
     printf("r_size: %ld %d %d %d\n",r_size,max_edge,i_edge,nRead);
     if((r_matrix= (int *)calloc(r_size,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - r_matrix\n");
       exit(1);
     }
     printf("reads: %d %d %d %ld %ld\n",tRead,nRead,nMat,n_head[tRead],tKmer);
     
     for(iMat=i_patch;iMat<nMat;iMat++)
     {
        int mat_low;
        int mat_hig;
        B64_long nKmer;

        mat_low = iMat*i_edge*1000000;
        mat_hig = (iMat+1)*i_edge*1000000;
        if(iMat == (nMat-1))
          mat_hig = tRead;

        memset(nameout_idex,'\0',100);
        sprintf(nameout_idex,"matrix.read.%03d",iMat); 
        
	if ((fpIDEX = fopen(nameout_idex, "rb")) == NULL) {
           printf("Error: can't open file %s\n.", nameout_idex);
           exit(1);
	}

        offset = 0;
        fseek(fpIDEX,offset,SEEK_END);
	nKmer = ftell(fpIDEX);
        fclose(fpIDEX);

	rsize = nKmer;
        rsize2 = max_edge;	
        rsize = rsize/rsize2;
	nRead = rsize;
        memset(nameout_hash,'\0',100);
        sprintf(nameout_hash,"matrix.link.%03d",iMat); 
        
	if ((fpHASH = fopen(nameout_hash, "rb")) == NULL) {
           printf("Error: can't open file %s\n.", nameout_hash);
           exit(1);
	}

	offset = 0;
        fseek(fpHASH,offset,0);
        sizesm=nKmer;
        fread(r_matrix,sizeof(int),sizesm,fpHASH);
        fclose(fpHASH);
        printf("patch: %d %d %ld\n",iMat,nRead,nKmer);

        for(i=0;i<nRead;i++)
	{
	   for(j=0;j<n_list[mat_low+i];j++)
	   {
	      long k_offset = i*max_edge;
	      rarray[n_head[mat_low+i]+j] = r_matrix[k_offset+j];
	   }
	}
     }

     memset(nameout_hash,'\0',100);
     sprintf(nameout_hash,"relation.matrix.link"); 
     
     if ((fpHASH = fopen(nameout_hash, "wb")) == NULL) {
        printf("Error: can't write to file %s\n.", nameout_hash);
        exit(1);
     }

     offset = 0;
     fseek(fpHASH,offset,0);
     sizesm=tKmer;
     fwrite(rarray,sizeof(int),sizesm,fpHASH);
     fclose(fpHASH);

     free(rarray);
     free(r_matrix);

     if((carray= (unsigned char *)calloc(l_size,sizeof(unsigned char))) == NULL)
     {
       printf("ssaha: calloc - cmatrix\n");
       exit(1);
     }
     if((c_matrix= (unsigned char *)calloc(r_size,sizeof(unsigned char))) == NULL)
     {
       printf("ssaha: calloc - c_matrix\n");
       exit(1);
     }
     printf("reads: %d %d %d %ld %ld\n",tRead,nRead,nMat,n_head[tRead],tKmer);
     
     for(iMat=i_patch;iMat<nMat;iMat++)
     {
        int mat_low;
        int mat_hig;
        B64_long nKmer;

        mat_low = iMat*i_edge*1000000;
        mat_hig = (iMat+1)*i_edge*1000000;
        if(iMat == (nMat-1))
          mat_hig = tRead;

        memset(nameout_idex,'\0',100);
        sprintf(nameout_idex,"matrix.read.%03d",iMat); 
        
	if ((fpIDEX = fopen(nameout_idex, "rb")) == NULL) {
           printf("Error: can't open file %s\n.", nameout_idex);
           exit(1);
	}

        offset = 0;
        fseek(fpIDEX,offset,SEEK_END);
	nKmer = ftell(fpIDEX);
        fclose(fpIDEX);

	rsize = nKmer;
        rsize = rsize/max_edge;
	nRead = rsize;
        printf("patch: %d %d %ld\n",iMat,nRead,nKmer);

        memset(nameout_idex,'\0',100);
        sprintf(nameout_idex,"matrix.read.%03d",iMat); 
        
	if ((fpIDEX = fopen(nameout_idex, "rb")) == NULL) {
           printf("Error: can't open file %s\n.", nameout_idex);
           exit(1);
	}

        offset = 0;
        fseek(fpIDEX,offset,0);
        sizesm=nKmer;
        fread(c_matrix,sizeof(unsigned char),sizesm,fpIDEX);
        fclose(fpIDEX);

        for(i=0;i<nRead;i++)
	{
	   for(j=0;j<n_list[mat_low+i];j++)
	   {
	      long k_offset = i*max_edge;
	      carray[n_head[mat_low+i]+j] = c_matrix[k_offset+j];
	   }
	}
     }

     memset(nameout_idex,'\0',100);
     sprintf(nameout_idex,"relation.matrix.read"); 
     
     if ((fpIDEX = fopen(nameout_idex, "wb")) == NULL) {
        printf("Error: can't write to file %s\n.", nameout_idex);
        exit(1);
     }

     offset = 0;
     fseek(fpIDEX,offset,0);
     sizesm=tKmer;
     fwrite(carray,sizeof(unsigned char),sizesm,fpIDEX);
     fclose(fpIDEX);
}

/*   Subroutine to screen out paired reads with unique kmers   */
/* =============================================  */
void Screen_Process(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i,k,n_Sbase=SEG_LEN;
     B64_long mhistc,num_maxmers,offset,sizesm;
     B64_long mhistcc = 0,tKmer,nPatchs,iPatchs,nsize[10];
     B64_long *patch_array,*patch_head;
     int **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     unsigned char **cmatrix(long nrl,long nrh,long ncl,long nch);
     char nameout_hash[100],nameout_idex[100],nameout_list[100],nameout_size[100];
     int *patch_index,*patch_list,*patch_list2,*read_mask,*read_mask2,*read_mask3; 
     int num_sect,n_reads; //id_read will never exceed mhist
     int nshift = (n_Sbase<<1)-10;//2^10=nsorts
     int j,ac,nsorts,nclip,*map_sort,tRead,iMat,nMat;
     FILE *fp,*fpHASH,*fpIDEX,*fpLIST,*fpSIZE;
     B64_long l_size,all_size;
     void ArraySort_Mix(int n, B64_long *arr, int *brr);
     fasta *segg,*seq,*seqp;
     B64_long Size_q_pdata,totalBases;
     int num_seqque;
     char *pdata;

     mhistc = 0;
     l_size = lsize*1000000L;
     nsorts = num_reads;
     for(i=0;nsorts > 10;i++)
     {     
        nsorts = nsorts>>1; 
     }        
     nsorts = 1L<<i;
     nshift = (n_Sbase<<1)-i;//2^i=nsorts
     if((map_sort= (int *)calloc(1000,sizeof(int))) == NULL)
     {        
       printf("phusion: calloc - map_sort\n");
       exit(1); 
     }        
     if((patch_list= (int *)calloc(nsorts + 1,sizeof(int))) == NULL)
     {        
       printf("phusion: calloc - patch_list\n");
       exit(1); 
     }        
     if((patch_list2= (int *)calloc(nsorts + 1,sizeof(int))) == NULL)
     {        
       printf("phusion: calloc - patch_list\n");
       exit(1); 
     }        
     if((patch_head= (B64_long *)calloc(nsorts + 1,sizeof(B64_long))) == NULL)
     {        
       printf("phusion: calloc - patch_head\n");
       exit(1);
     }
     if((hist = (B64_long *)calloc(400,sizeof(long))) == NULL)
     {
       printf("ERROR ssaha_init: calloc - hist\n");
       exit(1);
     }

     num_maxmers = 0;
     tKmer = 0;
     tRead = 0;
     for(ac=args;ac<argc;ac++)
     {
	memset(nameout_size,'\0',100);
	sprintf(nameout_size,"%s.size",argv[ac]);
	
	if ((fpSIZE = fopen(nameout_size, "rb")) == NULL) {
           printf("Error: can't open file %s.\n", nameout_size);
           exit(1);
	}

	offset = 0;
	fseek(fpSIZE,offset,0);
        sizesm=2;
        fread(nsize,sizeof(B64_long),sizesm,fpSIZE);
        printf("size: %ld %ld %ld\n",nsize[0],nsize[1],l_size);
        mhistc = nsize[1];
        tKmer = tKmer + mhistc;
        if(mhistc > num_maxmers)
          num_maxmers = mhistc;
        fclose(fpSIZE);
        tRead = tRead + nsize[0];
     }
     mhistcc = num_maxmers + 10;

     for(ac=args;ac<argc;ac++)
     {
        memset(nameout_list,'\0',100); 
	sprintf(nameout_list,"%s.list",argv[ac]); 
        
	if ((fpLIST = fopen(nameout_list, "rb")) == NULL) {
           printf("Error: can't open file %s.\n", nameout_list);
           exit(1);
	}

        offset = 0;
	fseek(fpLIST,offset,0);
        sizesm=nsorts;
        fread(patch_list2,sizeof(int),sizesm,fpLIST);
	fclose(fpLIST);

        printf("size2: %ld\n",mhistcc);

        for(i=0;i<nsorts;i++)
        {
           patch_list[i] = patch_list[i] + patch_list2[i]; 
        }
     }
     nPatchs = tKmer;
     nPatchs = nPatchs/l_size;
     nPatchs = nPatchs+1; 
     patch_head[0]=0;
     iPatchs = 0;
     map_sort[iPatchs] = 0;
     for(i=0;i<=nsorts;i++) 
     {
        if(i>0)
          patch_head[i] = patch_head[i-1] + patch_list[i-1];
	if((patch_head[i]-patch_head[map_sort[iPatchs]])>l_size)
	{
	  iPatchs++;
	  map_sort[iPatchs] = i;
	printf("got: %ld %d %ld %d %d %ld\n",i,patch_list[i],patch_head[i],map_sort[iPatchs],nsorts,patch_head[i]-patch_head[map_sort[iPatchs]]);
	}
     }
     map_sort[nPatchs] = nsorts; 

     all_size = 0;
     for(j=0;j<nPatchs;j++)
     {
        int n_st,n_ed;
	B64_long gsize;
	
	n_st = map_sort[j];
	n_ed = map_sort[j+1];
	gsize = patch_head[n_ed-1] - patch_head[n_st];
	if(gsize > all_size)
	  all_size = gsize;
	printf("size3: %ld %ld %d %d %d %ld %ld\n",gsize,all_size,n_st,n_ed,map_sort[j],tKmer,patch_head[nsorts]);
     }

     if((patch_array= (B64_long *)calloc(all_size+5000000,sizeof(B64_long))) == NULL)
     {
       printf("ssaha: calloc - patch_array\n");
       exit(1);
     }
     if((patch_index= (int *)calloc(all_size+5000000,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }

     /*   allocate memery for arrays to calculate the link number   */


     nRead = i_edge*1000000 + 50000;
     printf("reads: %d %d %ld\n",tRead,nRead,nPatchs);
     if((read_mask= (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - read_mask\n");
       exit(1);
     }
     if((read_mask2= (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - read_mask2\n");
       exit(1);
     }
     if((read_mask3= (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - read_mask2\n");
       exit(1);
     }
     if(tRead%(i_edge*1000000)==0)
       nMat = tRead/(i_edge*1000000);
     else
       nMat = tRead/(i_edge*1000000) + 1;
     for(iMat=i_patch;iMat<nMat;iMat++)
     {
        int mat_low;
        int mat_hig;

        mat_low = iMat*i_edge*1000000;
        mat_hig = (iMat+1)*i_edge*1000000;
        if(iMat == (nMat-1))
          mat_hig = tRead;
        nclip = 0;
        printf("mat: %d %d %d %d\n",iMat,nMat,mat_low,mat_hig);
        for(iPatchs=0;iPatchs<nPatchs;iPatchs++)
        {
           int n_st,n_ed;
	   n_st = map_sort[iPatchs];
	   n_ed = map_sort[iPatchs+1];
           memset(nameout_hash,'\0',100); 
           memset(nameout_idex,'\0',100);
           memset(nameout_size,'\0',100);
           sprintf(nameout_hash,"sorted.hash.%03ld",iPatchs); 
           sprintf(nameout_idex,"sorted.idex.%03ld",iPatchs);
           sprintf(nameout_size,"sorted.size.%03ld",iPatchs);
           printf("Read size file ...\n");
           
	   if ((fpSIZE = fopen(nameout_size, "rb")) == NULL) {
              printf("Error: can't open file %s.\n", nameout_size);
              exit(1);
	   }

           offset = 0;
           fseek(fpSIZE,offset,0);
           sizesm=3;
           fread(nsize,sizeof(B64_long),sizesm,fpSIZE);
	   tKmer = nsize[2];
           fclose(fpSIZE);

	   if ((fpHASH = fopen(nameout_hash, "rb")) == NULL) {
              printf("Error: can't open file %s\n.", nameout_hash);
              exit(1);
	   }

           offset = 0;
           fseek(fpHASH,offset,0);
           sizesm=tKmer;
           fread(patch_array,sizeof(B64_long),sizesm,fpHASH);
           fclose(fpHASH);
           
	   if ((fpIDEX = fopen(nameout_idex, "rb")) == NULL) {
              printf("Error: can't open file %s\n.", nameout_idex);
              exit(1);
	   }

           offset = 0;
           fseek(fpIDEX,offset,0);
           sizesm=tKmer;
           fread(patch_index,sizeof(int),sizesm,fpIDEX);
           fclose(fpIDEX);
           n_st = nsize[0];
	   n_ed = nsize[1];
           printf("file: %ld %ld %ld %d %d %s\n",nsize[0],nsize[1],tKmer,mat_low,mat_hig,nameout_size);
	   for(i=n_st;i<n_ed;i++)
	   {

	      if(patch_list[i]<1)
	        continue;
	      for(k=0;k<patch_list[i];k++)
	      {
                 B64_long kmer_long;
                 int kk,ki;
                 j = k+1;
                 num_sect = 1;
		 offset = patch_head[i]-patch_head[n_st];
                 while((j<patch_list[i])&&(patch_array[offset+j]==patch_array[offset+k]))
                 {
                   num_sect++;
                   j++;
                 }
                 kmer_long = patch_array[offset+k];
                 
                 if((iMat==i_patch)&&(num_sect<127))
                   hist[num_sect]++;
                 if(num_sect <= num_cut)
                 {
		   kk = patch_index[offset+k];
		   read_mask[kk]++;
                   for(ki=0;ki<num_sect;ki++)
                   {
		      kk = patch_index[offset+k+ki];
		      read_mask2[kk]++;
                   }
                 }
                 k = j-1;
              }
	   }
        }
     }
/*     for(i=0;i<tRead;i++)
     {
        printf("%5d %12d\n",i,read_mask[i]);
     }   */
     n_reads = 0;
     printf("mask done!\n");
     for(ac=args;ac<argc;ac++)
     {
        printf("read file: sss %s\n",argv[ac]);
	
        if ((fp = fopen(argv[ac], "rb")) == NULL) {
	    printf("Error: Cannot open file %s\n",argv[ac]);
	    exit(1);
	}
	
        fseek(fp, 0, SEEK_END);
        Size_q_pdata = ftell(fp) + 1;
        fclose(fp);

        if ((pdata = (char*)calloc(Size_q_pdata, sizeof(char))) == NULL) {
            printf("Error: calloc pdata\n");
	    exit(1);
	}
	
        num_seqque = extractFastq(argv[ac],pdata,Size_q_pdata);
	
        if ((segg = (fasta*)calloc((num_seqque+1), sizeof(fasta))) == NULL) {
            printf("Error: calloc segg\n");
	    exit(1);
	}
	
        if ((seq = decodeFastq(argv[ac], &num_seqque, &totalBases, pdata,Size_q_pdata, segg)) == NULL) {
            printf("Error: no query data found.\n");
	    exit(1);
	}
	
        nRead = num_seqque;
	memset(nameout_hash,'\0',100);
	sprintf(nameout_hash,"%s.screen",argv[ac]);
	
	if ((fpHASH = fopen(nameout_hash, "w")) == NULL) {
	    printf("Error: can't write to %s.\n", nameout_hash);
	    exit(1);
	}

        for(i=0;i<nRead;i++)
        {
	   if((read_mask[i+n_reads] <= num_ukmer)&&(read_mask2[i+n_reads] <= num_copy))
	   {
	     int rc,seq_len;
	     seqp = seq+i;
	     seq_len = seqp->length;
	     fprintf(fpHASH,"@%s\n",seqp->name);
	     for(rc=0;rc<seq_len;rc++)
	        fprintf(fpHASH,"%c",seqp->data[rc]);
	     fprintf(fpHASH,"\n");
	     fprintf(fpHASH,"+\n");
	     for(rc=0;rc<seq_len;rc++)
	        putc(seqp->qual[rc]+041,fpHASH);
	     fprintf(fpHASH,"\n");
	   }
        }
	fclose(fpHASH);
	n_reads = n_reads + nRead;
        free(pdata);
        free(segg);
     }
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Edge_Process(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i,k,n_Sbase=SEG_LEN;
     B64_long mhistc,num_maxmers,offset,sizesm;
     B64_long mhistcc = 0,tKmer,nPatchs,iPatchs,nsize[10];
     B64_long *patch_array,*patch_head;
     int **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     unsigned char **cmatrix(long nrl,long nrh,long ncl,long nch);
     char nameout_hash[100],nameout_idex[100],nameout_list[100],nameout_size[100];
     int *patch_index,*patch_list,*patch_list2,*n_list; 
     int id_read[64],num_sect; //id_read will never exceed mhist
     int nshift = (n_Sbase<<1)-10;//2^10=nsorts
     int j,ac,nsorts,nclip,*map_sort,tRead,iMat,nMat;
     int *r_matrix;
     unsigned char *c_matrix;
     FILE *fpHASH,*fpIDEX,*fpLIST,*fpSIZE;
     B64_long l_size,all_size,r_size;
     void ArraySort_Mix(int n, B64_long *arr, int *brr);

     mhistc = 0;
     l_size = lsize*1000000L;
     nsorts = num_reads;
     for(i=0;nsorts > 10;i++)
     {     
        nsorts = nsorts>>1; 
     }        
     nsorts = 1L<<i;
     nshift = (n_Sbase<<1)-i;//2^i=nsorts
     if((map_sort= (int *)calloc(1000,sizeof(int))) == NULL)
     {        
       printf("phusion: calloc - map_sort\n");
       exit(1); 
     }        
     if((patch_list= (int *)calloc(nsorts,sizeof(int))) == NULL)
     {        
       printf("phusion: calloc - patch_list\n");
       exit(1); 
     }        
     if((patch_list2= (int *)calloc(nsorts,sizeof(int))) == NULL)
     {        
       printf("phusion: calloc - patch_list\n");
       exit(1); 
     }        
     if((patch_head= (B64_long *)calloc(nsorts,sizeof(B64_long))) == NULL)
     {        
       printf("phusion: calloc - patch_head\n");
       exit(1);
     }
     if((hist = (B64_long *)calloc(400,sizeof(long))) == NULL)
     {
       printf("ERROR ssaha_init: calloc - hist\n");
       exit(1);
     }

     num_maxmers = 0;
     tKmer = 0;
     tRead = 0;
     for(ac=args;ac<argc;ac++)
     {
	memset(nameout_size,'\0',100);
	sprintf(nameout_size,"%s.size",argv[ac]);
	
	if ((fpSIZE = fopen(nameout_size, "rb")) == NULL) {
           printf("Error: can't open file %s.\n", nameout_size);
           exit(1);
	}

	offset = 0;
	fseek(fpSIZE,offset,0);
        sizesm=2;
        fread(nsize,sizeof(B64_long),sizesm,fpSIZE);
        printf("size: %ld %ld %ld\n",nsize[0],nsize[1],l_size);
        mhistc = nsize[1];
        tKmer = tKmer + mhistc;
        if(mhistc > num_maxmers)
          num_maxmers = mhistc;
        fclose(fpSIZE);
        tRead = tRead + nsize[0];
     }
     mhistcc = num_maxmers + 10;

     for(ac=args;ac<argc;ac++)
     {
        memset(nameout_list,'\0',100); 
	sprintf(nameout_list,"%s.list",argv[ac]); 
        
	if ((fpLIST = fopen(nameout_list, "rb")) == NULL) {
           printf("Error: can't open file %s.\n", nameout_list);
           exit(1);
	}

        offset = 0;
	fseek(fpLIST,offset,0);
        sizesm=nsorts;
        fread(patch_list2,sizeof(int),sizesm,fpLIST);
	fclose(fpLIST);

        printf("size2: %ld\n",mhistcc);

        for(i=0;i<nsorts;i++)
        {
           patch_list[i] = patch_list[i] + patch_list2[i]; 
        }
     }
     nPatchs = tKmer;
     nPatchs = nPatchs/l_size;
     nPatchs = nPatchs+1; 
     patch_head[0]=0;
     iPatchs = 0;
     map_sort[iPatchs] = 0;
     for(i=0;i<=nsorts;i++) 
     {
        if(i>0)
          patch_head[i] = patch_head[i-1] + patch_list[i-1];
	if((patch_head[i]-patch_head[map_sort[iPatchs]])>l_size)
	{
	  iPatchs++;
	  map_sort[iPatchs] = i;
	printf("got: %ld %d %ld %d %d %ld\n",i,patch_list[i],patch_head[i],map_sort[iPatchs],nsorts,patch_head[i]-patch_head[map_sort[iPatchs]]);
	}
     }
     map_sort[nPatchs] = nsorts; 

     all_size = 0;
     for(j=0;j<nPatchs;j++)
     {
        int n_st,n_ed;
	B64_long gsize;
	
	n_st = map_sort[j];
	n_ed = map_sort[j+1];
	gsize = patch_head[n_ed-1] - patch_head[n_st];
	if(gsize > all_size)
	  all_size = gsize;
	printf("size3: %ld %ld %d %d %d %ld %ld\n",gsize,all_size,n_st,n_ed,map_sort[j],tKmer,patch_head[nsorts]);
     }

     if((patch_array= (B64_long *)calloc(all_size+5000000,sizeof(B64_long))) == NULL)
     {
       printf("ssaha: calloc - patch_array\n");
       exit(1);
     }
     if((patch_index= (int *)calloc(all_size+5000000,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }

     /*   allocate memery for arrays to calculate the link number   */


     nRead = i_edge*1000000 + 50000;
     r_size = nRead;
     r_size = r_size*max_edge;
     printf("reads: %d %d %ld\n",tRead,nRead,nPatchs);
     if((r_matrix= (int *)calloc(r_size,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - r_matrix\n");
       exit(1);
     }
     if((c_matrix= (unsigned char *)calloc(r_size,sizeof(unsigned char))) == NULL)
     {
       printf("ssaha: calloc - c_matrix\n");
       exit(1);
     }
     if((n_list= (int *)calloc(tRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - n_list\n");
       exit(1);
     }
     if(tRead%(i_edge*1000000)==0)
       nMat = tRead/(i_edge*1000000);
     else
       nMat = tRead/(i_edge*1000000) + 1;
     for(iMat=i_patch;iMat<nMat;iMat++)
     {
        int mat_low;
        int mat_hig;

        mat_low = iMat*i_edge*1000000;
        mat_hig = (iMat+1)*i_edge*1000000;
        if(iMat == (nMat-1))
          mat_hig = tRead;
        nclip = 0;
        printf("mat: %d %d %d %d\n",iMat,nMat,mat_low,mat_hig);
        for(iPatchs=0;iPatchs<nPatchs;iPatchs++)
        {
           int n_st,n_ed;
	   n_st = map_sort[iPatchs];
	   n_ed = map_sort[iPatchs+1];
           memset(nameout_hash,'\0',100); 
           memset(nameout_idex,'\0',100);
           memset(nameout_size,'\0',100);
           sprintf(nameout_hash,"sorted.hash.%03ld",iPatchs); 
           sprintf(nameout_idex,"sorted.idex.%03ld",iPatchs);
           sprintf(nameout_size,"sorted.size.%03ld",iPatchs);
           printf("Read size file ...\n");
           
	   if ((fpSIZE = fopen(nameout_size, "rb")) == NULL) {
              printf("Error: can't open file %s.\n", nameout_size);
              exit(1);
	   }

           offset = 0;
           fseek(fpSIZE,offset,0);
           sizesm=3;
           fread(nsize,sizeof(B64_long),sizesm,fpSIZE);
	   tKmer = nsize[2];
           fclose(fpSIZE);

	   if ((fpHASH = fopen(nameout_hash, "rb")) == NULL) {
              printf("Error: can't open file %s\n.", nameout_hash);
              exit(1);
	   }

           offset = 0;
           fseek(fpHASH,offset,0);
           sizesm=tKmer;
           fread(patch_array,sizeof(B64_long),sizesm,fpHASH);
           fclose(fpHASH);

 	   if ((fpIDEX = fopen(nameout_idex, "rb")) == NULL) {
              printf("Error: can't open file %s\n.", nameout_idex);
              exit(1);
	   }

           offset = 0;
           fseek(fpIDEX,offset,0);
           sizesm=tKmer;
           fread(patch_index,sizeof(int),sizesm,fpIDEX);
           fclose(fpIDEX);
           n_st = nsize[0];
	   n_ed = nsize[1];
           printf("file: %ld %ld %ld %d %d %s\n",nsize[0],nsize[1],tKmer,mat_low,mat_hig,nameout_size);
	   for(i=n_st;i<n_ed;i++)
	   {

	      if(patch_list[i]<1)
	        continue;
	      for(k=0;k<patch_list[i];k++)
	      {
                 B64_long kmer_long;
                 int kk,ki;
                 j = k+1;
                 num_sect = 1;
		 offset = patch_head[i]-patch_head[n_st];
                 while((j<patch_list[i])&&(patch_array[offset+j]==patch_array[offset+k]))
                 {
                   num_sect++;
                   j++;
                 }
                 kmer_long = patch_array[offset+k];
                 if((iMat==i_patch)&&(num_sect<127))
                   hist[num_sect]++;
                 if((num_sect <= mhist)&&(num_sect>=mhist2))
                 {
                   for(kk=0;kk<num_sect;kk++)
		      id_read[kk] = patch_index[offset+kk+k];
                   for(kk=0;kk<num_sect;kk++)
                   {
                      int idk=id_read[kk];;
                      int cidk=n_list[idk];
		      B64_long k_offset;

		      if((idk<mat_low)||(idk>=mat_hig))
		        continue;
		      k_offset = (idk-mat_low);
		      k_offset = k_offset*max_edge;
                      if(cidk<(max_edge-2))
                      {
                        for(ki=0;ki<num_sect;ki++)
                        {
                           if(ki!=kk)
                           {
                             int idi=id_read[ki];
                             int match=0;
                             int km=cidk;
                             for(km=0;km<cidk;km++)
			     {
			        if(r_matrix[k_offset+km]==idi)
			        {
			          if(c_matrix[k_offset+km]<63)
			            c_matrix[k_offset+km]++;
			           match=1;
			           break;
			        }
			     }
                             if(match==0)
                             {
                               r_matrix[k_offset+cidk]=idi;
                               c_matrix[k_offset+cidk] = 1;
                               cidk++;
                               n_list[idk]++;
                             }
                           }
                        }
                      } 
                      else 
                        nclip++;
                   }
		   num_sect=0;
                 }
                 k = j-1;
              }
	   }
           printf("matrix here: %ld %ld %ld\n",nsize[0],nsize[1],tKmer);
        }

//        for(i=mat_low;i<mat_hig;i++)
//           printf("%d %5d %12ld\n",iMat,i,n_list[i]);
        tKmer = max_edge;
        tKmer = tKmer*(mat_hig-mat_low);
        printf("matrix out: %d %d %d %ld\n",max_edge,mat_low,mat_hig,tKmer);
        memset(nameout_hash,'\0',100);
        sprintf(nameout_hash,"matrix.link.%03d",iMat); 
        
	if ((fpHASH = fopen(nameout_hash, "wb")) == NULL) {
           printf("Error: can't write to file %s\n.", nameout_hash);
           exit(1);
	}

	offset = 0;
        fseek(fpHASH,offset,0);
        sizesm=tKmer;
        fwrite(r_matrix,sizeof(int),sizesm,fpHASH);
        fclose(fpHASH);

        memset(nameout_idex,'\0',100);
        sprintf(nameout_idex,"matrix.read.%03d",iMat); 
        
	if ((fpIDEX = fopen(nameout_idex, "wb")) == NULL) {
           printf("Error: can't write to file %s\n.", nameout_idex);
           exit(1);
	}

        offset = 0;
        fseek(fpIDEX,offset,0);
        sizesm=tKmer;
        fwrite(c_matrix,sizeof(unsigned char),sizesm,fpIDEX);
        fclose(fpIDEX);
     }
     for(i=0;i<128;i++)
        printf("hist: %5ld %12ld\n",i,hist[i]);
     tKmer = 0;
     for(i=0;i<tRead;i++)
     {
        tKmer = tKmer + n_list[i];
//        printf("%5d %12ld\n",i,n_list[i]);
     }

     memset(nameout_list,'\0',100);
     printf("Output list file ... %d %ld\n",tRead,tKmer);
     sprintf(nameout_list,"matrix.list");
     
     if ((fpLIST = fopen(nameout_list, "wb")) == NULL) {
     	printf("Error: can't write to %s.\n", nameout_list);
	exit(1);
     }
     
     offset = 0;
     fseek(fpLIST,offset,0);
     sizesm=tKmer;
     fwrite(n_list,sizeof(int),sizesm,fpLIST);
     fclose(fpLIST);

     nsize[0] = tKmer;
     nsize[1] = tRead;
     memset(nameout_size,'\0',100);
     printf("Output size file ... %d %ld\n",tRead,tKmer);
     sprintf(nameout_size,"matrix.size");
     
     if ((fpSIZE = fopen(nameout_size, "wb")) == NULL) {
        printf("Error: can't write to file %s\n.", nameout_size);
        exit(1);
     }

     offset = 0;
     fseek(fpSIZE,offset,0);
     sizesm=2;
     fwrite(nsize,sizeof(B64_long),sizesm,fpSIZE);
     fclose(fpSIZE);
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Sort_Process(int nRead,int SEG_LEN,int n_depth,int mat_st,int mat_ed,int s_len,int pmod, char **argv, int argc, int args)
/* =============================================  */
{
     B64_long i,k,n_Sbase=SEG_LEN;
     B64_long mhistc,num_maxmers,offset,sizesm;
     B64_long mhistcc = 0,tKmer,nPatchs,iPatchs,nsize[10];
     B64_long *patch_array,*patch_s_array,*patch_head,*patch_head2,*ray;
     char nameout_hash[100],nameout_idex[100],nameout_list[100],nameout_size[100];
     int *patch_s_index,*patch_list,*patch_list2,*patch_index,*dex; 
     int nshift = (n_Sbase<<1)-10;//2^10=nsorts
     int j,ac,nsorts,*map_sort;
     FILE *fpHASH,*fpIDEX,*fpLIST,*fpSIZE;
     B64_long l_size,all_size;
     void ArraySort_Mix(int n, B64_long *arr, int *brr);

     mhistc = 0;
     l_size = lsize*1000000L;
     nsorts = num_reads;
     for(i=0;nsorts > 10;i++)
     {     
        nsorts = nsorts>>1; 
     }        
     nsorts = 1L<<i;
     nshift = (n_Sbase<<1)-i;//2^i=nsorts
     if((map_sort= (int *)calloc(1000,sizeof(int))) == NULL)
     {        
       printf("phusion: calloc - map_sort\n");
       exit(1); 
     }        
     if((patch_list= (int *)calloc(nsorts + 1,sizeof(int))) == NULL)
     {        
       printf("phusion: calloc - patch_list\n");
       exit(1); 
     }        
     if((patch_list2= (int *)calloc(nsorts + 1,sizeof(int))) == NULL)
     {        
       printf("phusion: calloc - patch_list2\n");
       exit(1); 
     }        
     if((patch_head= (B64_long *)calloc(nsorts + 1,sizeof(B64_long))) == NULL)
     {        
       printf("phusion: calloc - patch_head\n");
       exit(1);
     }
     if((patch_head2= (B64_long *)calloc(nsorts + 1,sizeof(B64_long))) == NULL)
     {        
       printf("phusion: calloc - patch_head2\n");
       exit(1);
     }

     num_maxmers = 0;
     tKmer = 0;
     for(ac=args;ac<argc;ac++)
     {
	memset(nameout_size,'\0',100);
	sprintf(nameout_size,"%s.size",argv[ac]);
	
	if ((fpSIZE = fopen(nameout_size, "rb")) == NULL) {
           printf("Error: can't open file %s.\n", nameout_size);
           exit(1);
	}

	offset = 0;
	fseek(fpSIZE,offset,0);
        sizesm=2;
        fread(nsize,sizeof(B64_long),sizesm,fpSIZE);
        printf("size: %ld %ld %ld\n",nsize[0],nsize[1],l_size);
        mhistc = nsize[1];
        tKmer = tKmer + mhistc;
        if(mhistc > num_maxmers)
          num_maxmers = mhistc;
        fclose(fpSIZE);
     }
     mhistcc = num_maxmers + 10;
     if((patch_array= (B64_long *)calloc(mhistcc,sizeof(B64_long))) == NULL)
     {
       printf("ssaha: calloc - patch_array\n");
       exit(1);
     }
     if((patch_index= (int *)calloc(mhistcc,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }

     for(ac=args;ac<argc;ac++)
     {
        memset(nameout_list,'\0',100); 
	sprintf(nameout_list,"%s.list",argv[ac]); 
        
	if ((fpLIST = fopen(nameout_list, "rb")) == NULL) {
           printf("Error: can't open file %s.\n", nameout_list);
           exit(1);
	}

        offset = 0;
	fseek(fpLIST,offset,0);
        sizesm=nsorts;
        fread(patch_list2,sizeof(int),sizesm,fpLIST);
	fclose(fpLIST);

        printf("size2: %ld\n",mhistcc);

        for(i=0;i<nsorts;i++)
        {
           patch_list[i] = patch_list[i] + patch_list2[i]; 
        }
     }
     nPatchs = tKmer;
     nPatchs = nPatchs/l_size;
     nPatchs = nPatchs+1; 
     patch_head[0]=0;
     iPatchs = 0;
     map_sort[iPatchs] = 0;
     for(i=0;i<=nsorts;i++) 
     {
        if(i>0)
          patch_head[i] = patch_head[i-1] + patch_list[i-1];
	if((patch_head[i]-patch_head[map_sort[iPatchs]])>l_size)
	{
	  iPatchs++;
	  map_sort[iPatchs] = i-1;
//	printf("got: %d %d %ld %d %d %ld\n",i,patch_list[i],patch_head[i],map_sort[iPatchs],nsorts,patch_head[i]-patch_head[map_sort[iPatchs]]);
	}
     }
     map_sort[nPatchs] = nsorts; 
     memset(patch_list,0,4*nsorts); 

     all_size = 0;
     for(j=0;j<nPatchs;j++)
     {
        int n_st,n_ed;
	B64_long gsize;
	
	n_st = map_sort[j];
	n_ed = map_sort[j+1];
	gsize = patch_head[n_ed-1] - patch_head[n_st];
	if(gsize > all_size)
	  all_size = gsize;
	printf("size3: %ld %ld %d %d %d %ld %ld\n",gsize,all_size,n_st,n_ed,map_sort[j],tKmer,patch_head[nsorts]);
     }

     if((patch_s_array= (B64_long *)calloc(all_size+5000000,sizeof(B64_long))) == NULL)
     {
       printf("ssaha: calloc - patch_array\n");
       exit(1);
     }
     if((patch_s_index= (int *)calloc(all_size+5000000,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }
     for(j=i_patch;j<nPatchs;j++)
     {
        int n_st,n_ed,tRead;

	n_st = map_sort[j];
	n_ed = map_sort[j+1];
        tRead = 0;
        for(ac=args;ac<argc;ac++)
        {
           memset(nameout_hash,'\0',100); 
	   memset(nameout_idex,'\0',100);
	   memset(nameout_size,'\0',100);
	   memset(nameout_list,'\0',100);
	   sprintf(nameout_hash,"%s.hash",argv[ac]); 
	   sprintf(nameout_idex,"%s.idex",argv[ac]);
           sprintf(nameout_size,"%s.size",argv[ac]);
	   sprintf(nameout_list,"%s.list",argv[ac]);

	   
	   if ((fpSIZE = fopen(nameout_size, "rb")) == NULL) {
              printf("Error: can't open file %s.\n", nameout_size);
              exit(1);
	   }

	   offset = 0;
	   fseek(fpSIZE,offset,0);
	   sizesm=2;
	   fread(nsize,sizeof(B64_long),sizesm,fpSIZE);
	   fclose(fpSIZE);

           
	   if ((fpLIST = fopen(nameout_list, "rb")) == NULL) {
              printf("Error: can't open file %s.\n", nameout_list);
              exit(1);
	   }

           offset = 0;
  	   fseek(fpLIST,offset,0);
           sizesm=nsorts;
           fread(patch_list2,sizeof(int),sizesm,fpLIST);
	   fclose(fpLIST);

           patch_head2[0] = 0;
           for(i=0;i<=nsorts;i++) 
           {
              if(i>0)
                patch_head2[i] = patch_head2[i-1] + patch_list2[i-1];
	   }

	   sizesm = nsize[1];
           printf("get here: %d %ld %s %d %d %ld %ld\n",j,nPatchs,nameout_hash,n_st,n_ed,sizesm,patch_head[n_ed-1]);
           
	   if ((fpHASH = fopen(nameout_hash, "rb")) == NULL) {
              printf("Error: can't open file %s\n.", nameout_hash);
              exit(1);
	   }

           offset = 0;
           fseek(fpHASH,offset,0);
           fread(patch_array,sizeof(B64_long),sizesm,fpHASH);
           fclose(fpHASH);
	   
	   if ((fpIDEX = fopen(nameout_idex, "rb")) == NULL) {
              printf("Error: can't open file %s\n.", nameout_idex);
              exit(1);
	   }

           offset = 0; 
           fseek(fpIDEX,offset,0);
           fread(patch_index,sizeof(int),sizesm,fpIDEX);
           fclose(fpIDEX);

           for(i=n_st;i<n_ed;i++)
	   {
	      for(k=0;k<patch_list2[i];k++)
	      {
                B64_long patch_pos = patch_head[i]-patch_head[n_st]+patch_list[i];
                B64_long patch_pos2 = patch_head2[i];

                patch_s_array[patch_pos] = patch_array[patch_pos2+k];
                patch_s_index[patch_pos] = patch_index[patch_pos2+k]+tRead;
                patch_list[i]++; 
	      }
	   }
	   tRead = tRead + nsize[0];
	   printf("kmer: %d %d %d %d %ld\n",n_st,tRead,patch_list[0],patch_list2[0],patch_head[1]);
        }
        ray = patch_s_array;
        dex = patch_s_index;
	tKmer = 0;
        for(i=n_st;i<n_ed;i++) 
	{
	   tKmer = tKmer + patch_list[i];
	   if(patch_list[i] > 1)
	   {
	     B64_long patch_offset = patch_head[i]-patch_head[n_st];
	       
	     ArraySort_Mix(patch_list[i],ray+patch_offset,dex+patch_offset);
	   }
	}
        memset(nameout_hash,'\0',100); 
        memset(nameout_idex,'\0',100);
        memset(nameout_size,'\0',100);
        sprintf(nameout_hash,"sorted.hash.%03d",j); 
        sprintf(nameout_idex,"sorted.idex.%03d",j);
        sprintf(nameout_size,"sorted.size.%03d",j);
        
	if ((fpHASH = fopen(nameout_hash, "wb")) == NULL) {
           printf("Error: can't write to file %s\n.", nameout_hash);
           exit(1);
	}

        printf("Output hashtable file ... %d %d %ld %ld kmers\n",n_st,n_ed,tKmer,patch_head[n_ed]-patch_head[n_st]);
        offset = 0;
        fseek(fpHASH,offset,0);
        sizesm=tKmer;
        fwrite(patch_s_array,sizeof(B64_long),sizesm,fpHASH);
        fclose(fpHASH);

        printf("Output read index file ...\n");
        
	if ((fpIDEX = fopen(nameout_idex, "wb")) == NULL) {
           printf("Error: can't write to file %s\n.", nameout_idex);
           exit(1);
	}

        offset = 0;
        fseek(fpIDEX,offset,0);
        sizesm=tKmer;
        fwrite(patch_s_index,sizeof(int),sizesm,fpIDEX);
        fclose(fpIDEX);

        printf("Output size file ... %d\n",tRead);
        
	if ((fpSIZE = fopen(nameout_size, "wb")) == NULL) {
           printf("Error: can't write to file %s\n.", nameout_size);
           exit(1);
	}

        offset = 0;
        fseek(fpSIZE,offset,0);
        sizesm=3;
	nsize[0] = n_st;
	nsize[1] = n_ed;
	nsize[2] = tKmer;
        fwrite(nsize,sizeof(B64_long),sizesm,fpSIZE);
        fclose(fpSIZE);
     }
}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int  **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
unsigned int     **uimatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=((nrh-nrl+1)/100 + 1)*100,ncol=nch-ncl+1;
        unsigned int  **m;
	B64_long nri,nrn=nrow/100;

        /* allocate pointers to rows        */
        if((m=(unsigned int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }

        /* allocate rows and set pointers to them        */
	/* allocate in 100 batches to use freed memory */
	nrl = 0;
	for(nri=0;nri<100;nri++,nrl+=nrn) {
           if((m[nrl]=(unsigned int *)calloc(nrn*ncol,sizeof(int)))==NULL)
           {
              printf("error imatrix: calloc error No. 2 \n");
              return(NULL);
           }

           for(i=1;i<nrn;i++)
              m[i+nrl]=m[i+nrl-1]+ncol;
	}
       /* return pointer to array of pointers to rows   */
        return m;
}

/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
B64_long     **limatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        B64_long  **m;

        /* allocate pointers to rows        */
        if((m=(B64_long **)calloc(nrow,sizeof(B64_long*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(B64_long *)calloc(nrow*ncol,sizeof(B64_long)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C   
*/
/* =============================== */
void ArraySort_Mix(int n, B64_long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
unsigned char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        unsigned char **cm;

        /* allocate pointers to rows        */
        if((cm=(unsigned char **)calloc(nrow,sizeof(unsigned char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(unsigned char *)calloc(nrow*ncol,sizeof(unsigned char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}

 
