import sys
import pysam
import random
import numpy

"""
Estimate median template length of a bam file.

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

bamfile = sys.argv[1]
fraction = float(sys.argv[2])

tlens = list()

for read in pysam.AlignmentFile(bamfile):
	if not read.is_unmapped and random.random() < fraction:
		tlens.append(read.template_length)

print numpy.median(numpy.abs(tlens))
