# Copyright (c) 2012, 2014 Genome Research Ltd.
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


# Illumina Clone Assembly System Makefile - an installation
# script for icas.

use ExtUtils::MakeMaker;

my @c_binaries = map {"src/$_"} qw{
convert_fastq
fastq2fasta 
get_excreads 
parse_cross_genome
phusion2
map_screen_pcg
map_screen2_pcg
merge_pcg
screen_get_seqreads
merge_get_seqreads
rename_fastq
trim_fastq
ssaha_reads
ssaha_merge
get_seqreads
};

my @other_binaries = map {"bin/$_"} qw{
icas
cross_genome
};

my @binaries = (@c_binaries, @other_binaries);

WriteMakefile(  NAME	  => 'icas',
    	    	EXE_FILES => \@binaries,
		ABSTRACT  => 'Illumina Clone Assembly System');
