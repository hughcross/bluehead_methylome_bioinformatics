#!/usr/bin/env python


## add /1 and /2 to reads in sam file (that were removed by bwa) 
## determine if it affects output of SSPACE perl script sam_bam2tab.pl


#infile = open('sortN_subL1A_blue6k_mp.sam', 'r')
#infile = open('sortN_blue_350bp_insert.sam', 'r')
#infile = open('sortN_blue_550bp_insert.sam', 'r')
#infile = open('sortN_blue_6k_mate_pair.sam', 'r')
infile = open('sortN_blue_8k_mate_pair.sam', 'r')
#infile = open('sortN_sub100k_blue350bp.sam', 'r')
#outfile = open('modnew_sortN_sub100k_blue350bp.sam', 'w')# this one using lists
#outfile = open('modset_sortN_subL1A_blue6k_mp.sam', 'w')# this one using sets
#outfile = open('mod_sortN_blue_350bp_insert.sam', 'w')# this one using sets
#outfile = open('mod_sortN_blue_550bp_insert.sam', 'w')# this one using sets
#outfile = open('mod_sortN_blue_6k_mate_pair.sam', 'w')# this one using sets
outfile = open('mod_sortN_blue_8k_mate_pair.sam', 'w')# this one using sets

#samples = [] # make this a set, does not have to be indexed
samples = set() ## try this too
for line in infile:
	if line.startswith('@'):
		outfile.write(line)
	else:
		samplename = line.split('\t')[0]
		if samplename in samples:
			newsamplename = samplename+'/2'
		else:
			newsamplename = samplename+'/1'
			#samples.append(samplename)
			samples.add(samplename)
		newline = line.replace(samplename,newsamplename)
		outfile.write(newline)


