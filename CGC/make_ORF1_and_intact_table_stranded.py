import sys
import pickle

exp_prob_pkls_list = sys.argv[1]
bam_info_list = sys.argv[2]
orf1_list = sys.argv[3]
orf2_list = sys.argv[4]
allowed_runthrough_fraction = float(sys.argv[5])

output_orf1_name = sys.argv[6]
output_intact_name = sys.argv[7]

orf1_intact = set()
for line in open(orf1_list):
	orf1_intact.add(line.strip())
orf2_intact = set()
for line in open(orf2_list):
	orf2_intact.add(line.strip())

exp_probs = dict()
seqs = set([])

for line in open(exp_prob_pkls_list):
	names_file, X_file = line.strip().split('\t')
	name = names_file.split('/')[-1][:-16]
	exp_probs[name] = dict(zip(pickle.load(open(names_file)),pickle.load(open(X_file))))
	seqs = seqs | set(exp_probs[name].keys())

l1pa_pairs = dict()
mapped_pairs = dict()

for line in open(bam_info_list):
	name = line.strip().split('/')[-1][:-4]
	baminfo = open(line.strip()).readlines()
	mapped_pairs[name] = int(baminfo[1])
	l1pa_pairs[name] = int(baminfo[2])

output_orf1 = open(output_orf1_name,'w')
output_intact = open(output_intact_name,'w')

print_string = "locus"
for name in exp_probs:
	print_string += "\t"+name

output_orf1.write (print_string+'\n')
output_intact.write (print_string+'\n')

completed = set()

for name in seqs:
	seq_name = '_'.join(name.split('_')[:-1])
	if seq_name in completed:
		continue
	else:
		completed.add(seq_name)
	print_string = seq_name.split('(')[0]
	only_name = seq_name+'_only'
	runon_name = seq_name+'_3prunon'
	senserunthrough_name = seq_name+'_senserunthrough'
	antisenserunthrough_name = seq_name+'_antisenserunthrough'
	for name in exp_probs:
		FPM = 0.0
		runthrough_FPM = 0.0
		if only_name in exp_probs[name]:
			FPM += exp_probs[name][only_name]*l1pa_pairs[name]/mapped_pairs[name]*10**6
		if runon_name in exp_probs[name]:
			FPM += exp_probs[name][runon_name]*l1pa_pairs[name]/mapped_pairs[name]*10**6
		if senserunthrough_name in exp_probs[name]:
			runthrough_FPM += exp_probs[name][runthrough_name]*l1pa_pairs[name]/mapped_pairs[name]*10**6
		if FPM>0 and FPM/(FPM+runthrough_FPM) > allowed_runthrough_fraction:
			print_string += '\t'+str(FPM)
		else:
			print_string += '\t0.0'
	if seq_name.split('(')[0][:-2] in orf1_intact:
		output_orf1.write(print_string+'\n')
	if seq_name.split('(')[0][:-2] in orf1_intact and seq_name.split('(')[0][:-2] in orf2_intact:
		output_intact.write(print_string+'\n')

output_orf1.close()
output_intact.close()
