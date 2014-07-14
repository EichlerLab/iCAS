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
* matrix.c - Code modified for Illumina Clone Assembly Pipeline                  *
*/

#include <stdlib.h>
#include <stdio.h>

#include "matrix.h"


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int **imatrix(long nrl, long nrh, long ncl, long nch) {
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    int  **m;

    /* allocate pointers to rows        */
    if ((m = (int **)calloc(nrow, sizeof(int*))) == NULL) {
       printf("error imatrix: calloc error No. 1 \n");
       return(NULL);
    }
    
    m += 0;
    m -= ncl;

    /* allocate rows and set pointers to them        */
    if((m[nrl] = (int *)calloc(nrow * ncol, sizeof(int))) == NULL) {
       printf("error imatrix: calloc error No. 2 \n");
       return(NULL);
    }
    
    m[nrl] += 0;
    m[nrl] -= nrl;

    for(i = nrl + 1; i <= nrh; i++) {
    	m[i] = m[i-1] + ncol;
    }
    
    /* return pointer to array of pointers to rows   */
    return m;
}


void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch) {
    free(m[nrl] + ncl);
    free(m + nrl);
}


/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
char **cmatrix(long nrl, long nrh, long ncl, long nch) {
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    char **cm;

    /* allocate pointers to rows        */
    if((cm = (char **)calloc(nrow, sizeof(char*))) == NULL) {
       printf("error cmatrix: calloc error No. 1 \n");
       return(NULL);
    }

    cm += 0;
    cm -= ncl;

    /* allocate rows and set pointers to them        */
    if((cm[nrl] = (char *)calloc(nrow * ncol, sizeof(char))) == NULL) {
       printf("error cmatrix: calloc error No. 2 \n");
       return(NULL);
    }
    
    cm[nrl] += 0;
    cm[nrl] -= nrl;

    for(i = nrl + 1; i <= nrh; i++) {
    	cm[i] = cm[i - 1] + ncol;
    }
    
    /* return pointer to array of pointers to rows   */
    return cm;
}


void free_cmatrix(char **c, long nrl, long nrh, long ncl, long nch) {
    free(c[nrl] + ncl);
    free(c + nrl);
}
