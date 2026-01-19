# On Python2 import cPickle for performance improvement, else import pickle (available to both Py2 and Py3).
try:
    import cPickle as pickle
except ImportError:
    import pickle
import sys

"""
Extract the LINE-1 transcript estimates from L1EM. This version is intended for use on CGC
to analyze TCGA data.

Copyright (C) 2021 Wilson McKerrow

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

X_est = dict(zip(pickle.load(open(sys.argv[1],'rb')),pickle.load(open(sys.argv[2],'rb'))))

proper_pairs_in_original_bam = float(sys.argv[3])

total = float(sys.argv[4])

written_seqs = set([])

print("family.category.locus.strand\tonly\t3prunon")

names = list(X_est.keys())

for name in names:
	if 'L1MdTf_' in name or 'L1MdGf_' in name or 'L1MdA_I' in name or 'L1MdA_II' in name or 'L1MdA_III' in name:
		seq_name = '_'.join(name.split('_')[:-1])
		if seq_name in written_seqs:
			continue
		written_seqs.add(seq_name)
		print_string = seq_name.split('(')[0]
		only_name = seq_name+'_only'
		if only_name not in X_est:
			X_est[only_name]=0.0
		only_pairs = total*X_est[only_name]
		runon_name = seq_name+'_3prunon'
		if runon_name not in X_est:
			X_est[runon_name]=0.0
		runon_pairs = total*X_est[runon_name]
		runthrough_name = seq_name+'_runthrough'
		if runthrough_name not in X_est:
			X_est[runthrough_name]=0.0
		runthrough_pairs = total*X_est[runthrough_name]
		if (only_pairs+runon_pairs > 10*runthrough_pairs) & (only_pairs+runon_pairs>100):
			print(seq_name.split('(')[0]+'\t'+str(only_pairs/proper_pairs_in_original_bam*10**6)+'\t'+str(runon_pairs/proper_pairs_in_original_bam*10**6))
