import pickle
import numpy
import sys
import datetime
from scipy import sparse
from multiprocessing import Pool
import argparse

"""
This code takes as input the output of G_of_R.py and runs the EM algorithm to estimate
transcript abundances.

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

# Main calculation for the E step
def calculate_expcounts(G_of_R_pkl,X):
	G_of_R_file = open(G_of_R_pkl)
	G_of_R = pickle.load(G_of_R_file)
	G_of_R_file.close()
	if G_of_R == None:
		return 0.0,0.0
	L_of_R_mat = G_of_R.multiply(X)
	L_of_R = numpy.array(L_of_R_mat.sum(1))
	L_of_R_mat = L_of_R_mat[L_of_R[:,0]>=10**-200,:]
	L_of_R = L_of_R[L_of_R>=10**-200]
	L_of_R_inv = sparse.csr_matrix(1.0/L_of_R).transpose()
	exp_counts = L_of_R_mat.multiply(L_of_R_inv).sum(0)
	loglik = numpy.sum(numpy.log(L_of_R))
	if numpy.isfinite(loglik):
		return exp_counts,loglik
	else:
		return numpy.zeros(G_of_R.shape[1]),0.0

# Divide send each thread a chunk of the G_of_R pkl files.
def calculate_expcounts_chunk(input):
	G_of_R_pkl_list,X_len = input
	exp_counts = numpy.zeros(X_len.shape,dtype=numpy.float64)
	loglik = 0.0
	for G_of_R_pkl in G_of_R_pkl_list:
		this_exp_counts,this_loglik = calculate_expcounts(G_of_R_pkl,X_len)
		exp_counts += this_exp_counts
		loglik += this_loglik
	return exp_counts,loglik

# Parse commandline arguments
def GetArgs():

	def ParseArgs(parser):
		class Parser(argparse.ArgumentParser):
			def error(self, message):
				sys.stderr.write('error: %s\n' % message)
				self.print_help()
				sys.exit(2)

		parser.add_argument('-g', '--G_of_R_list',
							type=str,
							required=True,
							help='Text file listing paths to chunks of the G(R) matrix.')
		parser.add_argument('-l', '--TE_list',
							required=True,
							type=str,
							help='Text file listing the names of all transcripts. Output of G_of_R.py.')
		parser.add_argument('-s', '--stop_thresh',
							required=False,
							default=10**-7,
							type=float,
							help='Continue EM iterations until no transcription expression fraction (X_i) changes by more than this value. [1e-7]')
		parser.add_argument('-r', '--report_every',
							required=False,
							default=100,
							type=int,
							help='Write X every 100 steps. [100]')
		parser.add_argument('-m', '--max_nEMsteps',
							required=False,
							default=10000,
							type=int,
							help='Terminate if threshold has not been reached after this many EM steps [10000]')
		parser.add_argument('-t', '--nThreads',
							required=False,
							default=16,
							type=int,
							help='Divide E step into this many threads. [16]')
		parser.add_argument('-p', '--prefix',
							required=False,
							type=str,
							default='',
							help='If specified, this prefix will be used for output files.')
		return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.G_of_R_list, args.TE_list, args.stop_thresh, args.report_every, args.max_nEMsteps, args.nThreads, args.prefix


def main():
	G_of_R_list, TE_list, stop_thresh, report_every, max_nEMsteps, nThreads, prefix = GetArgs()

	# All the transcripts names in the same order as the G_of_R matrix columns
	TE_names = list()
	for name in open(TE_list):
		TE_names.append(name.strip().split('\t')[0])

	# Intial guess
	X = sparse.csr_matrix(numpy.ones((1,len(TE_names)),dtype=numpy.float64)/len(TE_names))

	# Split up the pickle files into a set for each thread.
	G_of_R_pkl_fulllist = list()
	for G_of_R_pkl in open(G_of_R_list):
		G_of_R_pkl_fulllist.append(G_of_R_pkl.strip())
	G_of_R_pkl_lists = list()
	listsize = len(G_of_R_pkl_fulllist)/nThreads
	nlistsp1 = len(G_of_R_pkl_fulllist)%nThreads
	k = 0
	for i in range(nlistsp1):
		G_of_R_pkl_lists.append(G_of_R_pkl_fulllist[k:k+listsize+1])
		k+=listsize+1
	for i in range(nlistsp1,nThreads):
		G_of_R_pkl_lists.append(G_of_R_pkl_fulllist[k:k+listsize])
		k+=listsize

	masterPool = Pool(processes = nThreads)

	# Run the EM steps
	for step in range(max_nEMsteps):
		starttime = datetime.datetime.now()
		exp_counts = numpy.zeros((1,len(TE_names)),dtype=numpy.float64)
		loglik = 0.0

		outputs = masterPool.map(calculate_expcounts_chunk,zip(G_of_R_pkl_lists,[X]*nThreads))
		for output in outputs:
			this_exp_counts,this_loglik = output
			exp_counts += this_exp_counts
			loglik += this_loglik

		last_X = X.copy()
		X = sparse.csr_matrix(exp_counts/numpy.sum(exp_counts))
		print(str(step)+" "+str(numpy.max(numpy.abs(X.toarray()-last_X.toarray())))+" "+str(loglik)+" "+str(datetime.datetime.now()-starttime))

		if (step+1) % report_every == 0:
			pickle.dump(X.toarray()[X.toarray() > 10**-10],open(prefix+'X_step_'+str(step+1)+'.pkl','w'),protocol=pickle.HIGHEST_PROTOCOL)
			pickle.dump(numpy.array(TE_names)[X.toarray()[0,:] > 10**-10],open(prefix+'names_step_'+str(step+1)+'.pkl','w'),protocol=pickle.HIGHEST_PROTOCOL)

		if numpy.max(numpy.abs(X.toarray()-last_X.toarray())) < stop_thresh:
			break

	# Output the final results
	pickle.dump(X.toarray()[X.toarray() > 10**-10],open(prefix+'X_final.pkl','w'),protocol=pickle.HIGHEST_PROTOCOL)
	pickle.dump(numpy.array(TE_names)[X.toarray()[0,:] > 10**-10],open(prefix+'names_final.pkl','w'),protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
	main()
