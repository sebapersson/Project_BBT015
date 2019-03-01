#!/bin/bash
# This script will run the entire analysis and produce all tables and figures,
# the result is stored in Results/Figures_copy and Results/Tables_copy.

# The different steps performed by the scripts are:
# Downloading the sample data (and storing in the Data-folder)
# Creating an extra annotation file for the plasmid. 
# Downloading bowtie-0.12.7 to the bin-directory.
# Aligning the reads to and storing the result in Intermediate/Count_mat
# Checking that the required R-packages are installed. 
# Creating count_matrix for plasmid and E. coli and storing in intermediate.
# Performing DESeq analysis. 

# If for example the sample-data already is present the script will skip
# that step. Due to relative paths being used the script has to be run
# from the Code/Shell directory. If some step fails the code will print
# an error message to stderr. 

# Check that the script is run from Shell directory
currentDir=${PWD##*/}

if [ ! $currentDir == "Shell" ]; then
    >&2 echo "The script most be run from Code/Shell directory"
    exit 1
fi 

# ------------------------ # Download sample data # ------------------------ #
# Current directory Code/Shell
# Run download data
echo "Downloading data"
./Download_data.sh

if [ $? != 0 ];then
   >&2 echo "Error when downloading data"
   exit 1
fi

# ------------------------ # Create annotation # -------------------------- #
# Move to Python directory
cd ../Python/

echo ""
echo "Creating extra annotation file"
./Annotate_plasmid.py

if [ $? != 0 ];then
   >&2 echo "Error when running Annotate_plasmid.py"
   exit 1
fi

# ------------------------ # Downloading bowtie # -------------------------- #
# Move to shell directory
cd ../Shell/

echo ""
echo "Downloading bowtie and setting up indices"
./Download_bowtie-0_12_7.sh 

if [ $? != 0 ];then
   >&2 echo "Error when running Download_bowtie-0_12_7.sh"
   exit 1
fi

# ---------------------- # Performing alignment # ------------------------- #
# Current directory Code/Shell
echo ""
echo "Aligning the reads"
./Alignment.sh 

if [ $? != 0 ];then
   >&2 echo "Error when running Alignment.sh"
   exit 1
fi

# -------------------------- # Check libraries # -------------------------- #
# Move to R directory
cd ../R
echo ""
echo "Checking libraries"
Rscript ./Check_libraries_installed.R 2> /dev/null 

if [ $? != 0 ];then
    >&2 echo "All packages aren't installed"
    >&2 echo "Exit code $?"
   exit 1
fi

echo "All libraries are installed"

# ----------------------- # Create count matrix # ------------------------- #
# Current directory Code/R
echo ""
echo "Creating count matrices"

# Don't display loading of libraries 
Rscript ./Create_count_mat_plasmid.R 2> /dev/null

if [ $? != 0 ];then
   >&2 echo "Error when running Create_count_mat_plasmid.R"
   exit 1
fi

# Don't display loading of libraries 
Rscript ./Create_count_mat_E_coli.R 2> /dev/null

if [ $? != 0 ];then
   >&2 echo "Error when running Create_count_mat_E_coli.R"
   exit 1
fi

# -------------------- # Differential analysis # ------------------------- #
# Current directory Code/R
echo ""
echo "Differential analysis E. coli"

Rscript ./Differential_analysis_E_coli.R 2> /dev/null

if [ $? != 0 ];then
   >&2 echo "Error when running ./Differential_analysis_E_coli.R"
   exit 1
fi

echo "Done, result can be found Results/Tables_copy and Results/Figures_copy"

echo ""
echo "Differential analysis plasmid"

Rscript ./Differential_analysis_plasmid.R 2> /dev/null

if [ $? != 0 ];then
   >&2 echo "Error when running Differential_analysis_plasmid.R"
   exit 1
fi

echo "Done, result can be found Results/Tables_copy and Results/Figures_copy"

# ------------------------ # Analysis done # ------------------------- #
echo ""
echo "The analysis ran with zero errors"
exit 0
