#!/bin/sh

# RunAllCutPastePATables.sh
#     written by H. Jeffries, July 31, 2006, Version 1.0

#  SCRIPT::CutPastePATables.sh
#      --Merge selected columns from all *.{ctb|ptb} files in subdirectories
#        in given top level directories

# expects an ENVIRONMENT Variable 'NET_BAL_HOME' to point to location
#   of a 'bin' directory where 'CutPastePATables.sh' script lives.

# Run this from a top level directory where PA output files
#  live, ie, where b1b.psito2n2  and  emis2.ptsrtceqag1  are present.
# It creates for each model PA focus output tree, text table of columns for PA focus regions
#    two tables of chemical PA parameters (37 lines each); 
#         c2 table is daily c3 table is peak ozone hr
#    two tables of physical PA parameters (26 lines each);
#         c2 table is daily c3 table is peak ozone hr
# The 37 lines is default, the c2 is default..

output=$(pwd)"/summaries/"

if [ ! -d ${output} ]
  then
    echo "ERROR::output summary directory ${output} does not exist."
    exit 3
fi

# first emissions scenario...
# chemical tables, daily and peak ozone...
# create file b1b_allctb2.txt -- b1b EI, last 37 lines from col 2 (daily) of PA physcial table 
$NET_BAL_HOME/bin/CutPastePATables.sh -o $output -pb1b_ b1b.psito2n2

# create file b1b_allctb3.txt -- b1b EI, last 37 lines from col 3 (peak ozone hr) of PA physcial table 
$NET_BAL_HOME/bin/CutPastePATables.sh -o $output -pb1b_ -c3 b1b.psito2n2

# physical tables, daily and peak ozone...
# create file b1b_allptb2.txt -- b1b EI, last 26 lines from col 2 (daily) of PA physcial table 
$NET_BAL_HOME/bin/CutPastePATables.sh -o $output -pb1b_ -tptb -l26 b1b.psito2n2

# create file b1b_allptb3.txt -- b1b EI, last 26 lines from col 3 (peak ozone hr) of PA physcial table 
$NET_BAL_HOME/bin/CutPastePATables.sh -o $output -pb1b_ -tptb -l26 -c3 b1b.psito2n2


# second emissions scenario...
# chemical tables, daily and peak ozone...
# create file em2_allctb2.txt -- AG EI, last 37 lines from col 2 (daily) of PA physcial table 
$NET_BAL_HOME/bin/CutPastePATables.sh -o $output -pem2_ emis2.ptsrtceqag1

# create file em2_allctb3.txt -- AG EI, last 26 lines from col 3 (peak ozone hr) of PA physcial table 
$NET_BAL_HOME/bin/CutPastePATables.sh -o $output -pem2_ -c3 emis2.ptsrtceqag1

# physical tables, daily and peak ozone...
# create file em2_allptb2.txt -- AG EI, last 26 lines from col 2 (daily) of PA physcial table 
$NET_BAL_HOME/bin/CutPastePATables.sh -o $output -pem2_ -tptb -l26 emis2.ptsrtceqag1

# create file em2_allptb3.txt -- AG EI, last 26 lines from col 3 (peak ozone hr) of PA physcial table 
$NET_BAL_HOME/bin/CutPastePATables.sh -o $output -pem2_ -tptb -l26 -c3 emis2.ptsrtceqag1


# third emissions scenario...
# chemical tables, daily and peak ozone...
# create file FORMeqVOC_allctb2.txt -- FORMeqVOC EI, last 37 lines from col 2 (daily) of PA physcial table 
$NET_BAL_HOME/bin/CutPastePATables.sh -o $output -pFORMeqVOC_ b1b.psito2n2.FORMeqVOC

# create file FORMeqVOC_allctb3.txt -- FORMeqVOC EI, last 37 lines from col 3 (peak ozone hr) of PA physcial table 
$NET_BAL_HOME/bin/CutPastePATables.sh -o $output -pFORMeqVOC_ -c3 b1b.psito2n2.FORMeqVOC

# physical tables, daily and peak ozone...
# create file FORMeqVOC_allptb2.txt -- FORMeqVOC EI, last 26 lines from col 2 (daily) of PA physcial table 
$NET_BAL_HOME/bin/CutPastePATables.sh -o $output -pFORMeqVOC_ -tptb -l26 b1b.psito2n2.FORMeqVOC

# create file FORMeqVOC_allptb3.txt -- FORMeqVOC EI, last 26 lines from col 3 (peak ozone hr) of PA physcial table 
$NET_BAL_HOME/bin/CutPastePATables.sh -o $output -pFORMeqVOC_ -tptb -l26 -c3 b1b.psito2n2.FORMeqVOC


