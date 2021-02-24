#!/usr/bin/env python3

import sys, argparse

#This script compares contact tracing data with our groundtruths and calculates the relative statistics.
#True detected represents close contacts that were found positive.


def createlist(infile, a):								#function that creates a list

	with open(infile, 'r') as fh:
		for line in fh:
			if line.startswith('>'):
				continue
			else:
				c = line.rstrip().split('\t')
				a.append(c)
		return a




def statistics(infile,outfile,GT):						#function that calculates True Detected, False Positive, False Negative, Sensitivity and Precision
    fout = open(outfile,'w')
    positive = set()
    close_contact = set()

    for line in infile:									#create a set for each of the following cases: Positive case, Close contact and Traveler from contact tracing data
        if line[3] == "POSITIVE CASE":
            positive.add(line[0])

        elif line[3] == "CLOSE CONTACT":
            close_contact.add(line[0])

        elif line[3] == "TRAVELER":
            continue


    column_GT_Baseline=set()							#create a set for each of the following cases: positive Swab, Groundtruth Baseline, Groundtruth Direct, Groundtruth Indirect
    column_GT_Direct=set()
    column_GT_Indirect=set()
    column_pos_Swab=set()

    for l in GT:

        if l[0] != ".":
            column_pos_Swab.add(l[0])

        if l[1] != ".":
            column_GT_Baseline.add(l[1])

        if l[2] != ".":
        	column_GT_Direct.add(l[2])

        if l[3] != ".":
            column_GT_Indirect.add(l[3])


    True_Detected = close_contact & positive										 #calculates all the statistics: true detected, false positive and false negative for each of the following cases: positive Swab, Groundtruth Baseline, Groundtruth Direct, Groundtruth Indirect
    False_Positive = close_contact - positive
    False_Negative_A = column_GT_Baseline-(column_GT_Baseline & close_contact)
    False_Negative_B = column_GT_Direct-(column_GT_Direct & close_contact)
    False_Negative_C = column_GT_Indirect-(column_GT_Indirect & close_contact)
    False_Negative_swab = column_pos_Swab-(column_pos_Swab & close_contact)

    TD = len(True_Detected)
    FP = len(False_Positive)
    FN_A = len(False_Negative_A)
    FN_B = len(False_Negative_B)
    FN_C = len(False_Negative_C)
    FN_Swab = len(False_Negative_swab)

    Sens_A = TD/(TD+FN_A)															#calculates Sensitivity for each case
    Sens_B = TD/(TD+FN_B)
    Sens_C = TD/(TD+FN_C)
    Sens_Swab = TD/(TD+FN_Swab)

    Precision = TD/(TD+FP)															#calculates Precision

    fout.write("\t"+"GT_A"+"\t"+"GT_B"+"\t"+"GT_C"+"\t"+"Swab")
    fout.write("\n")
    fout.write("True_Detected"+"\t"+str(TD)+"\t"+str(TD)+"\t"+str(TD)+"\t"+str(TD))
    fout.write("\n")
    fout.write("False_Positive"+"\t"+str(FP)+"\t"+str(FP)+"\t"+str(FP)+"\t"+str(FP))
    fout.write("\n")
    fout.write("False_Negative"+"\t"+str(FN_A)+"\t"+str(FN_B)+"\t"+str(FN_C)+"\t"+str(FN_Swab))
    fout.write("\n")
    fout.write("Sensitivity"+"\t"+str(Sens_A)+"\t"+str(Sens_B)+"\t"+str(Sens_C)+"\t"+str(Sens_Swab))
    fout.write("\n")
    fout.write("Precision"+"\t"+str(Precision)+"\t"+str(Precision)+"\t"+str(Precision)+"\t"+str(Precision))
    fout.write("\n")



def main(args):
	x = []
	y = []
	cluster = createlist(args['in'], x)
	GT = createlist(args['GT'],y)
	stat = statistics(cluster, args['out'], GT)



if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Parser')
	parser.add_argument('-in', metavar='Contact tracing data',  help='contact tracing data')
	parser.add_argument('-out', metavar='Table with Statistics', help='Table with Statistics')
	parser.add_argument('-GT', metavar='GroundTruth', help='File with four columns: Swab_pos, GT_A, GT_B, GT_C')
	args = vars(parser.parse_args())

	main(args)



