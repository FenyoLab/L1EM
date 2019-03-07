import pysam
import sys
import numpy
import cPickle
from scipy import sparse
import datetime
import argparse

"""
This script reads through a bam file resulting from a bwa aln alignment to the L1EM reference.
The output is a sparse matrix in which the rows are reads, the columns are transcripts
and the entries are the likelihood of that read arising from that transcript.
The matrix is pickled and saved. The column names are writted to a text file.

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

"""
This class stores relevant information about a read's potential alignment as a dictionary
with references names as keys and as list of potential alignments to that reference name
as values.
"""
class read_alignments(object):
	def __init__(self, alignment,rnames,P):
		self.alignments = dict()
		self.alignments[rnames[alignment.rname]] = [alignment_at_name(alignment.reference_start,alignment.is_reverse,P)]
	# Add a new alignment, passing a pysam aligned_segnment object.
	def add(self, alignment,rnames,P):
		if rnames[alignment.rname] not in self.alignments:
			self.alignments[rnames[alignment.rname]] = [alignment_at_name(alignment.reference_start,alignment.is_reverse,P)]
		else:
			self.alignments[rnames[alignment.rname]].append(alignment_at_name(alignment.reference_start,alignment.is_reverse,P))
	# Add a new alignment, passing the output of parseXA.
	def addXA(self,refname,start,is_reverse,P):
		if refname not in self.alignments:
			self.alignments[refname] = [alignment_at_name(start,is_reverse,P)]
		else:
			self.alignments[refname].append(alignment_at_name(start,is_reverse,P))

# Stores position, strand and likelihood for an alignment.
class alignment_at_name(object):
	def __init__(self,start,is_reverse,P):
		self.start = start
		self.is_reverse = is_reverse
		self.P = P

# Read command line arguments
def GetArgs():

	def ParseArgs(parser):
		class Parser(argparse.ArgumentParser):
			def error(self, message):
				sys.stderr.write('error: %s\n' % message)
				self.print_help()
				sys.exit(2)

		parser.add_argument('-b', '--bamfile',
							type=str,
							required=True,
							help='Bam to generate alignments from. Required.')
		parser.add_argument('-e', '--error_prob',
							required=False,
							default=0.01,
							type=float,
							help='Probability of an alignment mismatch. [0.01]')
		parser.add_argument('-m', '--max_start2start_len',
							required=False,
							default=500,
							type=int,
							help='Maximium distance between read starts to be considered concordant. [500]')
		parser.add_argument('-r', '--reads_per_pickle',
							required=False,
							default=12500,
							type=int,
							help='Split output into chunks of this many reads. [12500]')
		parser.add_argument('-p', '--prefix',
							required=False,
							default='G_of_R',
							type=str,
							help='Prefix for output file(s) [G_of_R]')
		parser.add_argument('-n', '--NMdiff',
							required=False,
							default=2,
							type=int,
							help='Ignore alignments with edit distance that exceed the best alignment by more than this number. [2]')
		parser.add_argument('-i', '--insert_mean',
							required=True,
							type=float,
							help='Median template length. Required.')
		parser.add_argument('--flanking',
							required=False,
							default=400,
							type=int,
							help='Number of flanking bases included on each end of repeats in reference fasta. [400]')
		parser.add_argument('--as_start',
							required=False,
							default=500,
							type=int,
							help='Position of the antisense TSS in L1. [500]')
		parser.add_argument('-w', '--wiggle',
							required=False,
							default=20,
							type=int,
							help='Extend L1 annotation this many bases in both directions. [20]')
		parser.add_argument('--min_len',
							required=False,
							default=500,
							type=int,
							help='When alignments probabilities are normalized for element length take max of elements length and this value. [500]')
		parser.add_argument('--min_exon_len',
							required=False,
							default=100,
							type=int,
							help='When alignments probabilities are normalized for exon length take max of elements length and this value. [100]')
		return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.bamfile, args.error_prob, args.max_start2start_len, args.reads_per_pickle, args.prefix, args.NMdiff, args.insert_mean, args.flanking, args.as_start,args.wiggle, args.min_len, args.min_exon_len

"""
Takes as input alignments (read_alignments class) of two paired reads and returns a sparse
row matrix with the likelihoods of all properly paired alignments.
"""
def get_concardant_alignments(alignments1,alignments2,max_start2start_len,rnames_index,rlens,insert_mean,nreps,read_length,flanking,as_start,wiggle,min_len,min_exon_len):
	this_G_of_R = numpy.zeros(5*nreps)
	for refname in alignments1.alignments:
		if refname not in alignments2.alignments:
			continue
		for aln1 in alignments1.alignments[refname]:
			for aln2 in alignments2.alignments[refname]:
				if aln1.is_reverse == aln2.is_reverse:
					continue
				if max(aln1.start,aln2.start)-min(aln1.start,aln2.start) <= max_start2start_len:
					has_5pUTR = refname.split('.')[1]=='1'
					if refname.split('.')[1]=='2':
						this_G_of_R[2*nreps+rnames_index[refname]] += aln1.P*aln2.P/(max(rlens[rnames_index[refname]]-insert_mean,min_exon_len))
						continue
					is_sense = not aln2.is_reverse
					within_5p = min(aln1.start,aln2.start) > flanking -wiggle
					within_3p = max(aln1.start,aln2.start)+read_length < rlens[rnames_index[refname]]-flanking +wiggle
					overlap_element = max(aln1.start,aln2.start)+read_length > flanking and min(aln1.start,aln2.start) < rlens[rnames_index[refname]]-flanking
					if not overlap_element:
						continue
					if is_sense:
						this_G_of_R[rnames_index[refname]] += aln1.P*aln2.P/(rlens[rnames_index[refname]]-2*flanking+insert_mean+2*wiggle)
					if not is_sense:
						this_G_of_R[nreps+rnames_index[refname]] += aln1.P*aln2.P/(rlens[rnames_index[refname]]-2*flanking+insert_mean+2*wiggle)
					if  within_5p and within_3p and is_sense and has_5pUTR:
						this_G_of_R[2*nreps+rnames_index[refname]] += aln1.P*aln2.P/max(rlens[rnames_index[refname]]-2*flanking-insert_mean+2*wiggle,min_len)
					if  within_5p and is_sense and has_5pUTR:
						this_G_of_R[3*nreps+rnames_index[refname]] += aln1.P*aln2.P/(rlens[rnames_index[refname]]-2*flanking+2*wiggle)
					if has_5pUTR and rlens[rnames_index[refname]] > flanking+as_start and max(aln1.start,aln2.start)+read_length < flanking+as_start and (not is_sense) and rlens[rnames_index[refname]] > flanking+as_start:
						this_G_of_R[4*nreps+rnames_index[refname]] += aln1.P*aln2.P/(as_start+insert_mean+wiggle)
	return sparse.csr_matrix(this_G_of_R)

# Parse secondary alignments in the XA tag from bwa aln.
def parseXA(alignments,XAtagdict,error_prob,maxNM,reversed):
	for aln in [x.split(',') for x in XAtagdict.split(';')[:-1]]:
		refname = aln[0]
		#if not reversed:
		#	is_reverse = aln[1][0] == '-'
		#else:
		#	is_reverse = aln[1][0] == '+'
		is_reverse = aln[1][0] == '-'
		start = int(aln[1][1:])
		cigarstring = aln[2]
		NM = int(aln[3])
		if NM <= maxNM and 'S' not in cigarstring and 'H' not in cigarstring:
			P = error_prob**NM
			alignments.addXA(refname,start,is_reverse,P)
	return alignments

def main():
	bamfile, error_prob, max_start2start_len, reads_per_pickle, prefix, NMdiff, insert_mean, flanking, as_start, wiggle, min_len, min_exon_len = GetArgs()
	
	pickle_num = 0
	
	bam = pysam.Samfile(bamfile, "rb")
	rnames = bam.references
	rlens = bam.lengths
	nreps = len(rnames)
	rnames_index = dict()
	for i in range(nreps):
		rnames_index[rnames[i]] = i
	
	# Write transcript (column) names
	TEnamefile = open(prefix+'_TE_list.txt','w')
	for i in range(nreps):
		TEnamefile.write(rnames[i]+'_senserunthrough'+'\t'+str(rlens[i]+2*flanking)+'\n')
	for i in range(nreps):
		TEnamefile.write(rnames[i]+'_antisenserunthrough'+'\t'+str(rlens[i]+2*flanking)+'\n')
	for i in range(nreps):
		TEnamefile.write(rnames[i]+'_only'+'\t'+str(rlens[i])+'\n')
	for i in range(nreps):
		TEnamefile.write(rnames[i]+'_3prunon'+'\t'+str(rlens[i]+flanking)+'\n')
	for i in range(nreps):
		TEnamefile.write(rnames[i]+'_antisense'+'\t'+str(flanking+as_start)+'\n')
	TEnamefile.close()
	
	read_id = None
	
	G_of_R = None
	G_of_R_list_file = open(prefix+'_list.txt','w')
	G_of_R_row = 0
	
	starttime = datetime.datetime.now()
	
	# Read through the name sorted bam file
	for alignment in bam:
		read_length = alignment.query_length
		# Throw out alignments that are unmapped, clipped or low quality
		if alignment.is_unmapped:
			continue
		if 'N' in alignment.cigarstring or 'S' in alignment.cigarstring or 'H' in alignment.cigarstring or 'P' in alignment.cigarstring or '=' in alignment.cigarstring or 'X' in alignment.cigarstring:
			continue
		if numpy.mean(alignment.query_qualities) < 30:
			continue
			
		if not read_id:
			read_id = alignment.qname
			new_read_id1 = True
			new_read_id2 = True
		
		# Once we have read all entries for a given query name, create a row for that fragment
		if read_id != alignment.qname:
			if not (new_read_id1 or new_read_id2):
				this_G_of_R = get_concardant_alignments(alignments1,alignments2,max_start2start_len,rnames_index,rlens,insert_mean,nreps,read_length,flanking,as_start,wiggle,min_len,min_exon_len)
				if this_G_of_R.nnz > 0:
					if G_of_R_row > 0:
						G_of_R = sparse.vstack([G_of_R,this_G_of_R])
					else:
						G_of_R = this_G_of_R
					G_of_R_row += 1
				# If necessary, break up matrix into multiple pickle files.
				if G_of_R_row >= reads_per_pickle:
					cPickle.dump(G_of_R,open(prefix+'.'+str(pickle_num)+'.pk2','wb'),protocol=cPickle.HIGHEST_PROTOCOL)
					G_of_R_list_file.write(prefix+'.'+str(pickle_num)+'.pk2\n')
					pickle_num += 1
					G_of_R_row = 0
					G_of_R = None
					print 'wrote',reads_per_pickle,'reads in',datetime.datetime.now()-starttime
					starttime = datetime.datetime.now()
			
			read_id = alignment.qname
			new_read_id1 = True
			new_read_id2 = True
		
		# Parse primary alignment
		# There's a bug in bwa samse (0.7.17) when writing NM tag for overlapping read pairs
		NMtag = dict(alignment.tags)['XM']
		for pair in alignment.cigartuples:
			NMtag += (pair[0]>0)*pair[1]
		P = error_prob**NMtag
		
		if alignment.is_read1:
			if new_read_id1:
				alignments1 = read_alignments(alignment,rnames,P)
				new_read_id1 = False
			else:
				alignments1.add(alignment,rnames,P)
		else:
			if new_read_id2:
				alignments2 = read_alignments(alignment,rnames,P)
				new_read_id2 = False
			else:
				alignments2.add(alignment,rnames,P)
		
		# Parse secondary alignments
		if 'XA' in dict(alignment.tags):
			if alignment.is_read1:
				alignments1 = parseXA(alignments1,dict(alignment.tags)['XA'],error_prob,NMtag+NMdiff,alignment.is_reverse)
			else:
				alignments2 = parseXA(alignments2,dict(alignment.tags)['XA'],error_prob,NMtag+NMdiff,alignment.is_reverse)
	
	# Make row for last read
	if not (new_read_id1 or new_read_id2):
		this_G_of_R = get_concardant_alignments(alignments1,alignments2,max_start2start_len,rnames_index,rlens,insert_mean,nreps,read_length,flanking,as_start,wiggle,min_len,min_exon_len)
		if this_G_of_R.nnz > 0:
			if G_of_R_row > 0:
				G_of_R = sparse.vstack([G_of_R,this_G_of_R])
			else:
				G_of_R = this_G_of_R
	
	# Write matrix to disk.
	cPickle.dump(G_of_R,open(prefix+'.'+str(pickle_num)+'.pk2','wb'),protocol=cPickle.HIGHEST_PROTOCOL)
	G_of_R_list_file.write(prefix+'.'+str(pickle_num)+'.pk2\n')
	print G_of_R_row
	
if __name__ == '__main__':
	main()
