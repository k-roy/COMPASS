#!/bin/bash

######################################### INSTALL PROGRAMS BEGIN ############################################
### EXTERNAL PROGRAMS ARE INSTALLED AS PART OF THE COMPASS CONDA ENVIRONMENT, EXCEPT SAMFIXCIGAR WHICH MUST BE DOWNLOADED FROM JVARKIT.GIT ###

# genome fasta downloaded from SGD (chromosome names reformatted to standard chrI, chrII etc, chrMito manually)
# genome gff downloaded from SGD (chromosome names reformatted to standard chrI, chrII etc, chrMito manually)

# download and install samfixcigar
SCRIPTS="/path/to/scripts/"
cd $SCRIPTS
git clone "https://github.com/lindenb/jvarkit.git"
cd jvarkit
./gradlew samfixcigar

conda env create -f COMPASS_environment.yml


######################################### INSTALL PROGRAMS END ##############################################
