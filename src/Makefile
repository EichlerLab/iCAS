# Copyright (c) 2011-2012, 2014 Genome Research Ltd.
#
# Author: Andrew Whitwham <aw7@sanger.ac.uk>
#
# This file is part of iCAS.
#
# iCAS is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>. 


# Illumina Clone Assembly System - C programs

CC = cc
CFLAGS = -O4 -Wall

LIB = -lm

OBJS =	convert_fastq.o \
	fastq2fasta.o \
	get_excreads.o \
	parse_cross_genome.o \
	fast.o \
	array_sort.o \
	matrix.o \
	phusion2_main.o \
	phusion2.o \
	map_screen_pcg.o \
	map_screen2_pcg.o \
	merge_pcg.o \
	screen_get_seqreads.o \
	merge_get_seqreads.o \
	rename_fastq.o \
	trim_fastq.o \
	ssaha_reads.o \
	ssaha_merge.o \
	get_seqreads.o
	

SRCS= $(OBJS:.o=.c)

INCS=	fasta.h \
	array_sort.h
	
		
PROGRAMS = convert_fastq fastq2fasta get_excreads parse_cross_genome \
 		phusion2 map_screen_pcg map_screen2_pcg merge_pcg screen_get_seqreads \
		merge_get_seqreads rename_fastq trim_fastq ssaha_reads ssaha_merge get_seqreads

all : $(PROGRAMS)

convert_fastq : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o convert_fastq convert_fastq.o fast.o $(LIB)

fastq2fasta : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o fastq2fasta fastq2fasta.o fast.o $(LIB)
	
get_excreads : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o get_excreads get_excreads.o fast.o array_sort.o matrix.o $(LIB)
	
parse_cross_genome : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o parse_cross_genome parse_cross_genome.o fast.o array_sort.o matrix.o $(LIB)
	
phusion2 : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o phusion2 phusion2.o fast.o phusion2_main.o $(LIB)

map_screen_pcg : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o map_screen_pcg map_screen_pcg.o fast.o array_sort.o matrix.o $(LIB)

map_screen2_pcg : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o map_screen2_pcg map_screen2_pcg.o fast.o array_sort.o matrix.o $(LIB)

merge_pcg : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o merge_pcg merge_pcg.o fast.o array_sort.o matrix.o $(LIB)

screen_get_seqreads : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o screen_get_seqreads screen_get_seqreads.o fast.o array_sort.o matrix.o $(LIB)

merge_get_seqreads : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o merge_get_seqreads merge_get_seqreads.o fast.o array_sort.o matrix.o $(LIB)

rename_fastq : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o rename_fastq rename_fastq.o fast.o array_sort.o matrix.o $(LIB)

trim_fastq : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o trim_fastq trim_fastq.o fast.o $(LIB)

ssaha_reads : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o ssaha_reads ssaha_reads.o fast.o $(LIB)
	
ssaha_merge : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o ssaha_merge ssaha_merge.o fast.o $(LIB)

get_seqreads : $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o get_seqreads get_seqreads.o fast.o array_sort.o matrix.o $(LIB)

clean:
	rm -f $(PROGRAMS) *.o
