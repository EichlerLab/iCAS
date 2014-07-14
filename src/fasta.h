/* 
* Copyright (c) 2012, 2014 Genome Research Ltd.
*
* Author: Andrew Whitwham <aw7@sanger.ac.uk>
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
* fasta.h - header file for the fasta/fastq handling functions
*/

#ifndef _SANGER_FASTA_H
#define _SANGER_FASTA_H 1

typedef struct
{
    char *name;
    char *name2;
    char *path;
    char *SCFname;
    int  length;
    char *data;
    char *qual;
    int  finished;
} fasta;

#define B64_long long int

fasta *decodeFastq (char *fname, int *nContigs, B64_long *tB, char* pdata, B64_long Size_pdata,fasta *segg);
int extractFastq(char *fname, char *pdata, B64_long Size_pdata);
fasta *splitFastq ( fasta *iseg, int inSeg, int tB, int *nReads, int length, int step);
void fastaLC (fasta *seg, int nSeg);
void fastaUC (fasta *seg, int nSeg);
int extract_fastq_from_file(char *fname, char **pd, B64_long *size_pd);
fasta *decode_fastq(char *pdata, B64_long size_pdata, int *nContigs, B64_long *tB);

#endif
