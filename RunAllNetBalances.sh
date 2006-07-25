#! /bin/sh

# RunAllNetBalances.sh
#     written by H. Jeffries, July 15, 2006, Version 1.0

# Run the "NetBalance_CB4.py script on each  *.ext file in
#   each ${subdir}/pypa_v2 directory of the list of top directories
#   given on the cmd line..

# Arguments Needed:
#   at least one "top level" dir having subdirs with *.ext files
#   Multiple such "top level" dirs can be given

# expects an ENVIRONMENT Variable 'NET_BAL_HOME' to point to location
#   of a 'bin' directory where 'net_balance_CB4.py' script lives.

echo "SCRIPT::RunAllNetBalance.sh
      --Run NetBalance_CB4.py on all *.ext in subdirectories
        in given top level directories
      ++ Version 1.0, July 15, 2006, Jeffries"
echo

if [ $# -lt 1 ]
 then
 	echo "ERROR:: need 1 or more top-level directory names on cmd line."
 	exit 1
fi
# get all cmdline arguments as a list
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
			# get list of *.ext files in this subdir
			extlist=$(ls *.ext)
			for ext in $extlist
				do
					# find the file name part of *.ext file
					srcname=${ext%*.ext}
					# make output file the same names
					outname=${srcname}.txt
					# run the python scrip on these files...
					python $NET_BAL_HOME/bin/net_balance_CB4.py $ext $outname
				done
			popd  1>/dev/null
			echo
		fi
	done
	echo "Finished All Subdirs in ${topdir}"
	
	popd  1>/dev/null
	fi
done

echo "Finished All Top Directories"

