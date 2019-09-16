import cPickle as pickle

"""
Extract the estimate of proper transcription of L1HS elements.

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

print "family.category.locus.strand\sesne\tantisense"

names = X_est.keys()

for name in names:
	if 'exon' not in name:
		seq_name = '_'.join(name.split('_')[:-1])
		if seq_name in written_seqs:
			continue
		written_seqs.add(seq_name)
		print_string = seq_name.split('(')[0]
		sense_name = seq_name+'_sense'
		if sense_name not in X_est:
			X_est[sense_name]=0.0
		print_string += '\t'+str(total*X_est[sense_name])
		antisense_name = seq_name+'_antisense'		
		if antisense_name not in X_est:
			X_est[antisense_name]=0.0
		print_string += '\t'+str(total*X_est[antisense_name])
		print print_string
