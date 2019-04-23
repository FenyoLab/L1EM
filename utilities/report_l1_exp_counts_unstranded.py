import cPickle as pickle
import sys

"""
Extract the LINE-1 transcript estimates from L1EM. This version is intended for use on CGC
to analyze TCGA data.

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

total = 0
for line in open('G_of_R_list.txt'):
	G_of_R = pickle.load(open(line.strip()))
	if G_of_R != None:
		total += pickle.load(open(line.strip())).shape[0]

X_est = dict(zip(pickle.load(open('names_final.pkl')),pickle.load(open('X_final.pkl'))))

written_seqs = set([])

print "family.category.locus.strand\tonly\t3prunon\tpassive"

names = X_est.keys()

for name in names:
	if 'exon' not in name:
		seq_name = '_'.join(name.split('_')[:-1])
		if seq_name in written_seqs:
			continue
		written_seqs.add(seq_name)
		print_string = seq_name.split('(')[0]
		only_name = seq_name+'_only'
		if only_name not in X_est:
			X_est[only_name]=0.0
		print_string += '\t'+str(total*X_est[only_name])
		runon_name = seq_name+'_3prunon'		
		if runon_name not in X_est:
			X_est[runon_name]=0.0
		print_string += '\t'+str(total*X_est[runon_name])
		runthrough_name = seq_name+'_runthrough'
		if runthrough_name not in X_est:
			X_est[runthrough_name]=0.0
		print_string += '\t'+str(total*X_est[runthrough_name])
		print print_string
