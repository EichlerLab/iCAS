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
* array_sort.h - header file for the array_sort funtions
*/

#ifndef _SANGER_ARRAY_SORT_H
#define _SANGER_ARRAY_SORT_H

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;
#define Max_N_NameBase 100

void ArraySort_Long(int n, long *arr);
void ArraySort_Int(int n, int *arr);
void ArraySort_Mix(int n, long *arr, int *brr);
void ArraySort_Int2(int n, int *arr, int *brr);
void ArraySort2_Int2(int n, int *arr, int *brr);
void ArraySort_Mix3(int n, long *arr, int *brr, int *crr);
void s_swap(char **Pair_Name, int i, int j);
void ArraySort_String(int n, char **Pair_Name, int *brr);

#endif
