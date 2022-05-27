import sys
import cPickle as pickle

exp_prob_pkls_list = sys.argv[1]
bam_info_list = sys.argv[2]
orf1_list = sys.argv[3]
orf2_list = sys.argv[4]
min_FPM = float(sys.argv[5])
allowed_runthrough_fraction = float(sys.argv[6])

l1pa_pairs = dict()
mapped_pairs = dict()

orf1_intact = set()
for line in open(orf1_list):
	orf1_intact.add(line.strip())
orf2_intact = set()
for line in open(orf2_list):
	orf2_intact.add(line.strip())

for line in open(bam_info_list):
	name = line.strip().split('/')[-1][:-4]
	baminfo = open(line.strip()).readlines()
	mapped_pairs[name] = int(baminfo[1])
	l1pa_pairs[name] = int(baminfo[2])

print('name\torf1_FPM\tORF2_FPM\tboth_FPM\tL1HS_expression_FPM\tL1HS_all_FPM')

for line in open(exp_prob_pkls_list):
	names_file, X_file = line.strip().split('\t')
	sample_name = names_file.split('/')[-1][:-16]
	exp_prob = dict(zip(pickle.load(open(names_file)),pickle.load(open(X_file))))
	orf1 = 0.0
	orf2 = 0.0
	both = 0.0
	L1HS_exp = 0.0
	L1HS_all = 0.0
	for transcript in exp_prob:
		if 'L1HS' in transcript:
			L1HS_all += exp_prob[transcript]*l1pa_pairs[sample_name]/mapped_pairs[sample_name]*10**6
		if 'only' not in transcript:
			continue
		seq_name = '_'.join(transcript.split('_')[:-1])
		only_name = seq_name+'_only'
		runon_name = seq_name+'_3prunon'
		runthrough_name = seq_name+'_senserunthrough'
		FPM = 0.0
		FPM += exp_prob[only_name]*l1pa_pairs[sample_name]/mapped_pairs[sample_name]*10**6
		if runon_name in exp_prob:
			FPM += exp_prob[runon_name]*l1pa_pairs[sample_name]/mapped_pairs[sample_name]*10**6
		if runthrough_name in exp_prob:
			runthrough_FPM = exp_prob[runthrough_name]*l1pa_pairs[sample_name]/mapped_pairs[sample_name]*10**6
		else:
			runthrough_FPM = 0.0
		FPM *= FPM >= min_FPM and runthrough_FPM/(runthrough_FPM+FPM) <= allowed_runthrough_fraction
		if seq_name.split('(')[0][:-2] in orf1_intact:
			orf1 += FPM
		if seq_name.split('(')[0][:-2] in orf2_intact:
			orf2 += FPM
		if seq_name.split('(')[0][:-2] in orf1_intact and seq_name.split('(')[0][:-2] in orf2_intact:
			both += FPM
		if 'L1HS' in seq_name:
			L1HS_exp += FPM
	print(sample_name +'\t'+ str(orf1) +'\t'+ str(orf2) +'\t'+ str(both) +'\t'+ str(L1HS_exp) +'\t'+ str(L1HS_all))
