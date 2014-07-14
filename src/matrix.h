/*
* matrix.h - header file for the matrix functions
*
* Andrew Whitwham (aw7@sanger.ac.uk)
* Wellcome Trust Sanger Institute
* April 2012
*
* This software has been created by Genome Research Limited (GRL).     
* GRL hereby grants permission to use, copy, modify and distribute     
* this software and its documentation for non-commercial purposes      
* without fee at the user's own risk on the basis set out below.      
* GRL neither undertakes nor accepts any duty whether contractual or   
* otherwise in connection with the software, its use or the use of     
* any derivative, and makes no representations or warranties, express 
* or implied, concerning the software, its suitability, fitness for   
* a particular purpose or non-infringement.                           
* In no event shall the authors of the software or GRL be responsible  
* or liable for any loss or damage whatsoever arising in any way       
* directly or indirectly out of the use of this software or its        
* derivatives, even if advised of the possibility of such damage.     
* Our software can be freely distributed under the conditions set out  
* above, and must contain this copyright notice. 
*/

#ifndef _SANGER_MATRIX_H
#define _SANGER_MATRIX_H

int  **imatrix(long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
char **cmatrix(long nrl, long nrh, long ncl, long nch);
void free_cmatrix(char **c, long nrl, long nrh, long ncl, long nch);

#endif

