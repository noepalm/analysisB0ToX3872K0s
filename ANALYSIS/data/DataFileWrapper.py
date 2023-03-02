import ROOT
import sys
import os
import argparse
import glob


parser = argparse.ArgumentParser (description = 'make corrected template')
parser.add_argument('-y', '--year', help='year of RUN 2 to be processed', type = int)
parser.add_argument('--sub', type = str , default = 'A', help='set to be processed')
parser.add_argument('--path', default = './data')

args = parser.parse_args ()

path2data = args.path + "/Data"+str(args.year)+"_"+args.sub+".txt"
print(path2data)
data_files = []
with open(path2data) as file2data:
    for line in file2data:
        print(line)
        if (line.startswith('#') or line.startswith('\n')): continue
        data_files.append(line.rstrip('\n'))

print(data_files)


files = glob.glob(data_files[0] + "/xNANO_data_*.root")
print(" ... in {0} there are {1} files".format(data_files[0], len(files)))
        

