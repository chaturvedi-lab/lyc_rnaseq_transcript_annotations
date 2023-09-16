#!/usr/bin/python3
import re
import sys
import copy
import argparse

import pandas as pd
from tqdm import tqdm
import goatools
from goatools import obo_parser
goont = obo_parser.GODag('./go-basic.obo')
import datetime


## -------------------------------------------- ##
## Author: Samridhi Chaturvedi
## Copyright: Copyright 2023
## Credits: [""]
## License: GPL3
## Version: 0.0.1
## Maintainer: Samridhi Chaturvedi
## Email: samridhi.chaturvedi@*.com
## Status: dev
## -------------------------------------------- ##

# -------------- #
# Pre-requisites #
# -------------- #
'''
Repository:
 - https://github.com/karwaan/lyc_rnaseq_transcript_annotations

Packages:
 - pip install git+git://github.com/tanghaibao/goatools.git

Data: 
 - wget http://purl.obolibrary.org/obo/go/go-basic.obo

Usage: 
 - python create_transcript_annotations.py --pos result1_all.txt --ann annot1631_sub.txt --out sub_transcript_annotations_table_1631.out
 - python create_transcript_annotations.py --help
'''

# --------------------- #
# Global variables conf #
# --------------------- #
''' 
Here we are getting SNPs in 1000bp of each start 
and stop position in the annotation table
'''
transcript_window = 10000			# snp position near or on this window range
sep_merge = '~'

# --------------------- #
# Function Definitions  #
# --------------------- #
''' Utility to test integer '''
def is_integer(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def printf(log):
	current_time = datetime.datetime.now()
	formatted_time = current_time.strftime('%Y-%m-%d %H:%M:%S.%f')
	print("{}: {}".format(formatted_time, log))

''' Load scaffold position to a dataframe '''
def load_scafpos(input_position_file):
	pos_df =  pd.read_csv(input_position_file, sep='\t', skipinitialspace=True)
	pos_df["scaffold"] = pos_df["scaffold"].astype(int)
	pos_df["start"]    = pos_df["start"].astype(int)
	pos_df["stop"]     = pos_df["stop"].astype(int)
	pos_df.sort_values(['cuffids', 'scaffold', 'start', 'stop'], axis=0, ascending=True, inplace=True, na_position="first")
	return pos_df

''' Load annotation table and sort the scaffolds and positions '''
def load_annot(input_annotations_file):
	col_names = ["seq_id","source","type","start","end","score","strand","phase","attributes"]
	gene_df = pd.read_csv(input_annotations_file, sep='\t', comment="#", names=col_names, header=None, skipinitialspace=True)
	gene_df["scaffold"] = gene_df["seq_id"].str.split('_')
	gene_df["scaffold"] = gene_df["scaffold"].str[1]
	gene_df["scaffold"] = gene_df["scaffold"].astype(int)
	gene_df.sort_values(['scaffold', 'start', 'end'], axis=0, ascending=True, inplace=True, na_position="first")
	return gene_df



''' Fetch the IPR and GO terms from the GOA file '''
def fetch_ipr_go_info(ann_info):
	#print ann_info
	ipr_match = re.findall(r'IPR[0-9]+', ann_info) #19,20 (concat IPR using ;)
	ipr_match_str = ';'.join(map(str,ipr_match))
	go_match =  re.findall(r'GO:[0-9]+', ann_info) #21,22 (concat GO using ;)
	go_match_str = ';'.join(map(str,go_match))
	#print "IPR: ", ipr_match_str
	#print "GO : ", go_match_str
	go_term = []
	for g in go_match:
		go_term.append(goont.get(g))
	go_term_str = ';'.join(map(str,go_term))  # Go term for corresponding id
	return ipr_match_str, go_match_str, go_term_str


''' Flatten and merge the annotations with the same snp position '''
def compute_format_result(scaffold, snppos, annotations):
	out_array = [int(0)] * 23
	out_array[0] = scaffold
	out_array[1] = snppos
	for ann in annotations:
		if ann is not None:
			#type,start stop,distance_type,annotation
			ann_array = ann.split(sep_snppos_tmp)
			# Check for ongene only once
			if ann_array[3] == 'on': 
				if ann_array[0] == 'gene':
					# Use the annotation start/stop of 
					out_array[2] = ann_array[1] # annstart
					out_array[3] = ann_array[2] # annstop
					out_array[4] = 1			# ongene
					ann_info = str(ann_array[4])
					ipr_match_str, go_match_str, go_term_str = fetch_ipr_go_info(ann_info)
					out_array[19] = ipr_match_str  	# ipr
					out_array[21] = go_match_str	# go number
					out_array[22] = go_term_str		# go text
				elif ann_array[0] == 'CDS':
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[5] = 1 			# oncds
				elif ann_array[0] == 'mRNA':
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[6] = 1
				elif ann_array[0] == 'exon':
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[7] = 1
				elif ann_array[0] == 'five_prime_UTR' or ann_array[0] == 'three_prime_UTR':
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[8] = 1
				elif ann_array[0] == 'protein_match' or ann_array[0] == 'expressed_sequence_match':
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop					
					out_array[9] = 1
				elif ann_array[0] == 'match' or ann_array[0] == 'match_part' :
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[10] = 1
			# Go here only if immediately previous annotation is near and not on
			elif ann_array[3] == 'left' or ann_array[3] == 'right':
				# Go here only if type matches and is not on for that type
				if ann_array[0] == 'gene' and out_array[4] == 0:
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					# Keep the IPR/GO information for near gene
					ann_info = str(ann_array[4])
					ipr_match_str, go_match_str, go_term_str = fetch_ipr_go_info(ann_info)
					# Merge logic
					if out_array[19] != 0:
						#print "\nAttempt Merge, snp: ",snppos," - ", out_array, '\t-\t', ann_array
						out_array[2] = out_array[2] + sep_merge + ann_array[1] # annstart 
						out_array[3] = out_array[3] + sep_merge + ann_array[2] # annstop 
						out_array[19] = out_array[19] + sep_merge + ipr_match_str 	# ipr text
						out_array[21] = out_array[21] + sep_merge + go_match_str  # go number
						out_array[22] = out_array[22] + sep_merge + go_term_str 	# go text
					else:
						# Persist the current 
						out_array[19] = ipr_match_str  	# ipr
						out_array[21] = go_match_str	# go number
						out_array[22] = go_term_str		# go text
					out_array[11] = 1 				# neargene
					out_array[18] = out_array[18] + 1
				elif ann_array[0] == 'CDS' and out_array[5] == 0:	  # not onCDS
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[12] = 1 			# oncds
				elif ann_array[0] == 'mRNA' and out_array[6] == 0:
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[13] = 1
				elif ann_array[0] == 'exon' and out_array[7] == 0:
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[14] = 1
				elif (ann_array[0] == 'five_prime_UTR' or ann_array[0] == 'three_prime_UTR') and (out_array[8] == 0):
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[15] = 1
				elif (ann_array[0] == 'protein_match' or ann_array[0] == 'expressed_sequence_match') and (out_array[9] == 0):
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop					
					out_array[16] = 1
				elif (ann_array[0] == 'match' or ann_array[0] == 'match_part') and (out_array[10] == 0) :
					if out_array[2] == 0: out_array[2] = ann_array[1] # annstart
					if out_array[3] == 0: out_array[3] = ann_array[2] # annstop
					out_array[17] = 1				

	#print "current out ---------->  SNP=", snppos, "\t", out_array
	# check the previous annot array item's position for concatenation
	# checking if type is gene, on & near distance type; only then proceed to concatenation
	resultarray.append(out_array)
	




''' Check if the overlap falls in ON or NEAR category based on position for the same scaffold '''
def check_overlap_window(cuff_df, ann_item_df):
	pos_start, pos_end = cuff_df.start, cuff_df.stop
	ann_start, ann_end = ann_item_df.start, ann_item_df.end
	# ON	: position window overlap on (a.) start, (b.) stop, (c.) both ann, (d.) within ann window
	if ((pos_start >= ann_start and pos_start <= ann_end) or 	# (a.)
		(pos_end   >= ann_start and pos_end   <= ann_end) or	# (b.)
		(pos_start >= ann_start and pos_end   <= ann_end) or	# (c.)
		(pos_start <= ann_start and pos_end   >= ann_end) 		# (d.)
		):
		return "ON"
	# NEAR	: position window (e.) upstream, (f.) downstream within the acceptable annotation window of transcript_window=10K 
	elif ((pos_start < ann_start and pos_end < ann_start and abs(pos_end - ann_start) <= transcript_window) or 	# (e.)
			(pos_start > ann_end and pos_end > ann_end and abs(pos_start - ann_end) <= transcript_window)	     	# (f.)
		):
		return "NEAR"
	# None
	return None


''' Compute cuff position window overlap in the annotation GFF data'''
def compute_position_overlap_annotation(cuff_df, gff_df, debug_mode=False):
	out_row = {
				"cuffids"		: None,
				"scaffold"		: None,
				"cuffstart"		: None,
				"cuffstop"		: None,
				"annstart"		: None,
				"annstop"		: None,
				"ongene"		: 0,
				"oncds"			: 0,
				"onmRNA"		: 0,
				"onexon"		: 0,
				"onte"			: 0,
				"onprotein"		: 0,
				"onmatch"		: 0,
				"neargene"		: 0,
				"nearcds"		: 0,
				"nearmRNA"		: 0,
				"nearexon"		: 0,
				"nearte"		: 0,
				"nearprotein"	: 0,
				"nearmatch"		: 0,
				"numneargenes"	: 0,
				"IPRnumber"		: None,
				"IPRtext"		: None,
				"GOnumber"		: None,
				"GOtext"		: None
			}
	out_row["cuffids"] 	 = cuff_df.cuffids
	out_row["scaffold"]  = cuff_df.scaffold
	out_row["cuffstart"] = cuff_df.start
	out_row["cuffstop"]  = cuff_df.stop
	
	#printf("Processing Transcript - CuffId:{}, ScaffId:{}, Pos({}, {})".format(cuff_df["cuffids"], cuff_df["scaffold"], cuff_df["start"], cuff_df["stop"]))

	for ann in gff_df.itertuples():	# Pandas Dataframe iterrows is very slow and needs to be changed
		idx = ann.Index
		_type = ann.type
		if _type not in ['gene', 'CDS', 'mRNA', 'exon', 
							'five_prime_UTR', 'three_prime_UTR', 'protein_match', 
							'expressed_sequence_match', 'match', 'match_part']:
			continue # Skip these types

		# Compare if cuff_window is (1.) ON or (2.) NEAR annotations for the same scaffold
		pos_scaff, ann_scaff = cuff_df.scaffold, ann.scaffold
		overlap_state = None
		if pos_scaff == ann_scaff:
			overlap_state = check_overlap_window(cuff_df, ann)

		if overlap_state is not None and overlap_state == "ON":
			if _type == "gene":
				out_row["ongene"] = 1											   					# ongene
				# Use the start/stop of ON gene annotation	
				out_row["annstart"] = str(ann.start) 												# annstart
				out_row["annstop"]  = str(ann.end)  												# annstop
				# Keep the IPR/GO information for near gene	
				ipr_match_str, go_match_str, go_term_str = fetch_ipr_go_info(ann.attributes)	
				out_row["IPRnumber"] = ipr_match_str  												# ipr
				#TODO out_row["IPRtext"] = Read this from entry.list file 
				out_row["GOnumber"]  = go_match_str													# go numb
				out_row["GOtext"]    = go_term_str													# go text

		elif overlap_state is not None and overlap_state == "NEAR":
			# Go here only if type matches and is not ON for that type
			if _type == "gene" and out_row["ongene"] == 0:
				out_row["neargene"] = 1											   					# neargene
				out_row["numneargenes"] = out_row["numneargenes"] + 1			   					# numneargenes
				# Keep the IPR/GO information for near gene
				ipr_match_str, go_match_str, go_term_str = fetch_ipr_go_info(ann.attributes)
				if out_row["IPRnumber"] is not None: 												# Merge logic
						out_row["annstart"]  = out_row["annstart"] + sep_merge + str(ann.start)  	# annstart
						out_row["annstop"]   = out_row["annstop"] + sep_merge + str(ann.end) 		# annstop
						out_row["IPRnumber"] = out_row["IPRnumber"] + sep_merge + ipr_match_str 	# ipr num
						out_row["GOnumber"]  = out_row["GOnumber"] + sep_merge + go_match_str		# go numb
						out_row["GOtext"]    = out_row["GOtext"] + sep_merge + go_term_str			# go text						
				else:																				# Persist the current 
					if out_row["annstart"] is None: out_row["annstart"] = str(ann.start)			# annstar
					if out_row["annstop"] is None:  out_row["annstop"]  = str(ann.end) 				# annstop
					out_row["IPRnumber"] = ipr_match_str  											# ipr num
					out_row["GOnumber"]  = go_match_str												# go numb
					out_row["GOtext"]    = go_term_str												# go text

		if debug_mode and idx%10000 == 0:
			printf("\t{} Annotation Item: {}, {}".format(idx, cuff_df.start, cuff_df.stop))
			printf("\tMatching position window overlap:{}".format(out_row))
			
	return out_row

''' Search annotations for cuff window position '''
def search_cuffposition_in_annotation(pos_df, gff_df):
	result_list = []
	# Operate on each cuff/transcript id 
	with tqdm(total = pos_df.shape[0]) as pbar: 
		for idx, cuffid in pos_df.iterrows():
			pbar.set_description("CuffId({}), ScaffId({})".format(cuffid["cuffids"], cuffid["scaffold"]))
			pbar.update(1)
			out_row_dict = compute_position_overlap_annotation(cuffid, gff_df, debug_mode=False)
			result_list.append(pd.DataFrame([out_row_dict]))

	# Accumulate and Report result on each cuffid
	result_df = pd.concat(result_list)
	return result_df




if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--pos", help="Add transcript scaffold positions file full path", default="result1_all.txt")
	parser.add_argument("--ann", help="Add genome annotations file full path", default="annot1631_sub.txt")
	parser.add_argument("--out", help="Add transcript annotation output file full path", default="sub_transcript_annotations_table_1631.out")
	args = parser.parse_args()
	if args.pos:  input_position_file = args.pos
	if args.ann:  input_annotations_file = args.ann
	if args.out:  output_file_name = args.out

	printf("Configs: \n\t POS: {} \n\t ANN: {} \n\t OUT:  {} \n".format(input_position_file, input_annotations_file, output_file_name))
	# Load the input data
	pos_df = load_scafpos(input_position_file)
	#sample_tr_df = pos_df[pos_df["scaffold"] == 1631]
	printf("Loaded scaffold position data: \n {}".format(pos_df.head(4)))
	gff_df = load_annot(input_annotations_file)
	printf("Loaded GFF data: \n {}".format(gff_df.head(4)))

	# compute the positions
	result_df = search_cuffposition_in_annotation(pos_df, gff_df)
	printf("Result data: \n {}".format(result_df.head()))

	# Wriite results to file
	result_df.to_csv(output_file_name, sep='\t', header=True, index=False)
	printf("Result writtent to file:  {}".format(output_file_name))