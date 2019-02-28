#!/usr/bin/python3
import sys
import re
import os

'''
This file extracts the Gene-name and product (gene description) from 
the pBIC_1a.gff3 file in the data-directory and writes the results 
to a tab-separated file in the Data/Reference_data/Plasmid folder.
'''

# Exit if script executed from wrong directory
if os.path.split(os.getcwd())[1] != "Python":
    print("The script isn't run from the Python directory")
    sys.exit(1)

# If file exists exit program 
filePathWrite = "../../Data/Reference_data/Plasmid/Annotation_plasmid.tsv"
if os.path.isfile(filePathWrite):
    print("Annotation file already exists")
    sys.exit(0)

# Open the annotation file, note the file most be run from Code/Python dir
fileName = "../../Data/Reference_data/Plasmid/pBIC_1a.gff3"

try:
    fp = open(fileName, "r")
    fileText = fp.read()
    fp.close()
except IOError:
    print("Couldn't open {}".format(fileName))
    sys.exit(1)


# Regex to match Gene-name and description from GFF3-file3
regex = re.compile(r"ID=gene-CHX41_[\d]*;Name=(.*?);.*?\n.*?product=(.*?);")
allMatches = re.findall(regex, fileText)

# Create string to export
stringToExport = "Gene_name\tDescription\n"
for entries in allMatches:
    stringToExport += entries[0] + "\t" + entries[1] + "\n"


try:
    fpWrite = open(filePathWrite, "w")
    fpWrite.write(stringToExport)
    fpWrite.close()

except IOError:
    print("Error writing the file")
    sys.exit(1)


# Everything went as expected
sys.exit(0)

