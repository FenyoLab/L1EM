import pysam
import sys

"""
Extract reads or pairs of reads that overlap a bed file.

Part of the L1-EM package.

Copyright (C) 2019 Wilson McKerrow

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

def main():
	bedfile = sys.argv[1]
	bamfile = sys.argv[2]
	outbamfile = sys.argv[3]
	if len(sys.argv) > 4:
		flanking = int(sys.argv[4])
	else:
		flanking = 400
	if len(sys.argv) > 5:
		maxNM = int(sys.argv[5])
	else:
		maxNM = 4
	
	inbam = pysam.AlignmentFile(bamfile,'rb')
	outbam = pysam.AlignmentFile(outbamfile,'wb',template=inbam)
	
	read_ids = set()
	for line in open(bedfile):
		chrom,start,stop = line.strip().split('\t')[:3]
		start = int(start)+flanking
		stop = int(stop)-flanking
		if chrom in inbam.references:
			for read in inbam.fetch(chrom,start,stop):
				if not read.is_unmapped:
					if not read.is_secondary and not read.is_supplementary and 'S' not in read.cigarstring and 'N' not in read.cigarstring and (not read.has_tag('NM') or read.get_tag('NM')<=maxNM):
						read_ids.add(read.query_name)
# 		if chrom[3:] in inbam.references:
# 			for read in inbam.fetch(chrom[3:],start,stop):
# 				if not read.is_secondary and not read.is_supplementary and 'S' not in read.cigarstring and 'N' not in read.cigarstring and read.get_tag('NM')<=3:
# 						read_ids.add(read.query_name)
# 		if '_' in chrom and chrom.split('_')[1].upper()+'.1' in inbam.references:
# 			for read in inbam.fetch(chrom.split('_')[1].upper()+'.1',start,stop):
# 				if not read.is_secondary and not read.is_supplementary and 'S' not in read.cigarstring and 'N' not in read.cigarstring and read.get_tag('NM')<=3:
# 					read_ids.add(read.query_name)
	
	inbam.close()
	inbam = pysam.AlignmentFile(bamfile,'rb')
	
	for read in inbam:
		if read.query_name in read_ids:
			if not read.is_secondary and not read.is_supplementary:
				outbam.write(read)
	
	inbam.close()
	outbam.close()

if __name__ == '__main__':
	main()
