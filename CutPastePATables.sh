#! /bin/sh

# CutPastePATables.sh
#     written by H. Jeffries, July 31, 2006, Version 1.0

# Use 'tail', 'cut' and 'paste' to create a merged file of 
#   col '2' {daily} or col '3' {peak hour} PA parameters.

# Arguments Needed:
#   -p outfilename prefix, eg, 'b1b' or 'emis2', default = ''
#   -t 'ctb' or 'ptb';   default is 'ctb'
#   -c n : select which column of *.ctb to cut; default = 2 == daily 3 == hourly
#   -l m : number of tail lines to keep in each table file; default = 37 = daily
#   -o outputpath : path of output file; default = $pwd

#   at least one "top level" dir having subdirs with *.ctb files
#   Multiple such "top level" dirs can be given

# needs no ENVIRONMENT Variables 

echo "SCRIPT::CutPastePATables.sh
      --Merge selected columns from all *.{ctb|ptb} files in subdirectories
        in given top level directories
      ++ Version 1.0, July 31, 2006, Jeffries"
echo

# process cmdline options...
args=`getopt p:t:c:l:o: $*`
if [ $? != 0 ]
then
	   echo 'Usage: PROG [-p PREFIX] [-c 3] [-l 26] [-o OUT]  dir1 [dir2 [dir 3]] '
	   echo 'Where:'
	   echo "  -p outfilename prefix, eg, 'b1b' or 'emis2', default = ''"
	   echo "  -t 'ctb' or 'ptb';   default is 'ctb'"
	   echo "  -c n : select which column of *.ctb to cut; default = 2"
	   echo "  -l n : number of tail lines to keep in each table file; default = 37"
	   echo "            for physical parameters use 26 lines"
	   echo "  -o outputpath : path of output file; default = $pwd"
	   exit 2
fi

set -- $args

# PROCESS OPTIONS:

# the default outfileprefix is ''
#   check args to see if different table ('ptb') is needed
# the default type of table is chemical, eg, 'ctb'
#   check args to see if different table ('ptb') is needed
# the default column to cut is 'daily' or col 2
#   check args to see if different column is requested
# the default outputpath is pwd
#   check args to see if different path is requested
prefix=''
xtb='ctb'
column=2
tlines=37
outputpath=$(pwd)"/"
for i
do
 case "$i"
   in
   -p)
     prefix="$2"; shift;
     shift;;
   -t)
     xtb="$2"; shift;
     shift;;
   -c)
     column="$2"; shift;
     shift;;
   -l)
     tlines="$2"; shift;
     shift;;
   -o)
     outputpath="$2"; shift;
     shift;;
   --)
     shift; break;;
 esac
done

echo "Prefix is $prefix"
echo "Table  is $xtb"
echo "Column is $column"
echo "Lines  is $tlines"
echo "Output is $outputpath"
echo

# set up for actions only done the first time executes.
firsttime=0

# now get all cmdline arguments as a list
topdirs=$*
# change to a toplevel dir and run rest of script
for topdir in ${topdirs}
do
 if [ ! -d ${topdir} ]
  then
    echo "ERROR::starting top level directory ${topdir} does not exist."
  else
	pushd ${topdir} 1>/dev/null
	echo  "Top Directory is:"
	pwd
	echo
	
	dirlist=$(ls)
	
	echo $dirlist
	echo
	
	for subdir in $dirlist
	do
	  if [ -d ${subdir} ]
		then
			echo "CHANGING TO ${subdir}"
			pushd ${subdir}/pypa_v2  1>/dev/null
			if [ $? != 0 ]
			 then
			   echo ERROR:: no 'pypa_v2' sub-dir in ${subdir}
			   exit 3
			fi
			# get list of *.xtb files in this subdir
			tblist=$(ls *.${xtb})
			for tb in $tblist
				do
				 # create a column of parameter names and paste into output file
				 if [ $firsttime -eq 0 ]
				  then
<<<<<<< .mine
				   tail -n ${tlines} $tb | cut -f 1 - > ${outputpath}col1.txt
				   paste ${outputpath}col1.txt > ${outputpath}${prefix}all${xtb}${column}.txt
				   rm -f ${outputpath}col1.txt
=======
				   tail -n ${tlines} $tb | cut -f 1 - > ${outputpath}${prefix}all${xtb}${column}.txt
>>>>>>> .r122
				   firsttime=1
				   echo "firsttime"
				 fi
				 ls $tb
				 tail -n ${tlines} $tb | cut -f $column - > ${outputpath}colx.txt
				 paste ${outputpath}${prefix}all${xtb}${column}.txt ${outputpath}colx.txt > ${outputpath}temp.txt
				 mv -f ${outputpath}temp.txt  ${outputpath}${prefix}all${xtb}${column}.txt
				 rm -f ${outputpath}colx.txt
				done
			popd  1>/dev/null
<<<<<<< .mine
			# add a first line to file that gives the files full name to id file on import
			echo "${prefix}all${xtb}${column}.txt"  > ${outputpath}head.txt
			cat  ${outputpath}head.txt  ${outputpath}${prefix}all${xtb}${column}.txt > ${outputpath}temp.txt
			mv -f ${outputpath}temp.txt   ${outputpath}${prefix}all${xtb}${column}.txt
			rm -f ${outputpath}head.txt
			echo
=======
>>>>>>> .r122
		fi
	done
<<<<<<< .mine
	rm -f ${outputpath}col1.txt
	
=======
>>>>>>> .r122
	echo "Finished All Subdirs in ${topdir}"
	
	popd  1>/dev/null
	fi
done
# add a first line to file that gives the files full name to id file on import
echo "${prefix}all${xtb}${column}.txt"  > ${outputpath}head.txt
cat  ${outputpath}head.txt  ${outputpath}${prefix}all${xtb}${column}.txt > ${outputpath}temp.txt
mv -f ${outputpath}temp.txt   ${outputpath}${prefix}all${xtb}${column}.txt
rm -f ${outputpath}head.txt
echo

echo "Finished All Top Directories"

