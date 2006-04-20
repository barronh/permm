#!/usr/bin/env python2.3
"""Make Files of Net Reactions and Chemical Parameters from IRR/IPR Data
Version 0.9 HEJ April 9, 2006, 
This file is
	$HeadURL$
	$LastChangedDate: 2006-04-19 14:19:58 -0400 (Wed, 19 Apr 2006) $
	$LastChangedRevision: 106 $
	$LastChangedBy: jeffries $


net_balance.py [options] EXTFILE OUTFILE

This version is hard-coded for CB4 Mechanism in CAMx

use -h or --help for details on options

"""

import os.path
import sys

import re
import copy   # for deepcopy
	
from optparse import OptionParser
	


SCRIPT_ID_STRING = "Net_Reaction_CB4.py, 2006-04-09 (c) 2006-04-09 (m)"

__version__ = "R69"


# create a cmd line parser...
usage  = "usage: %prog [options] EXTFILE OUTFILE"
version = "%%prog %s" % __version__
parser = OptionParser(usage=usage, version=version)

parser.add_option("-s", "--show",
	dest="show",
	default=False,
	action="store_true",
	help="Display graphs on screen")

parser.add_option("-o", "--outdir",
	dest="outdir",
	default=False,
	help="prepend output dir to outfilename")

# extent of output
parser.add_option("-q", "--quiet",
	action="store_false", 
	dest="verbose", 
	default=True,
	help="don't print status messages to stdout")

# output format
parser.add_option("-x", "--xml",
	action="store_const", 
	dest="output_format",
	const="xml",
	default='fixed_width',
	help="Format output as xml")
	
# get options and arguments
(options, args) = parser.parse_args()

# requires filenames as argument
if len(args) != 2:
	parser.error("Invalid number of arguments")
	
input_filename  = os.path.abspath(args[0])
output_filename = os.path.abspath(args[1])

# requires a valid ext file
if not os.path.exists(input_filename):
	parser.error("Input file does not exist")


# document all these options...
if options.verbose:
	print
	print "======================================================"
	print "Script ",  SCRIPT_ID_STRING
	print "This is        version ",  __version__
	print ""
	print "Reading from:"
	print "%s" % input_filename
	print ""
	print "Writing to  :"
	print "%s" % output_filename
	print ""

# ----- options now processed


# create regular expressions to read the input file..
time_re  = re.compile('Time =[0-9]{6}', re.IGNORECASE)
irr_re   = re.compile('\{\s*\d+\}\s+\d+', re.IGNORECASE)
ipr_re   = re.compile('"\w+\s*"\s*\d+\s*', re.IGNORECASE)
split_re = re.compile('[ ]+')


# function to get an hour's worth of data..
def get_irr_data(f):
	"""
	f -- opened irr&ipr file
	return (time (as integer), ir_rates (as list), ip_rates (as list of lists)
	if time < 0, end of file
	"""
	# initialize ir and ip storage...
	ir_time  = -1
	ir_rates = [0,]
	ip_rates = []
	#for i in range(0,23):  # make space for process rates for CB4's 24 species
		#ip_rates.append([])
	
	line = f.readline()
	if line[0] == '|':
		return (ir_time, ir_rates, ip_rates)
		
	# next line is first Time =
	if time_re.match(line) != None:   # if this is a 'Time =" line...
		time = line.split('=')[1][0:2] # pick 'HH' out of 'HH0000'
		ir_time = int(time)           # store hh
	else :
		print "ERROR:: did not find a time."
		sys.exit(1)
	
	# read in the !"Rxn No"    "Int Rate" header
	line = f.readline()
	
	# read the irr values
	while line[0] != ';' :
		line = f.readline()
		if irr_re.match(line) != None:    # if this is a { n} n.nnnn line ...
			ir_value = float(split_re.split(line.replace('{','').replace('}','').strip())[1])
			ir_rates.append(ir_value)
	

	line = f.readline()
	while line[0] != ';' :
		line = f.readline()
		if ipr_re.match(line) != None:
			tmp=[]
			for ip_value in split_re.split(line.replace('"','').replace('"','').strip())[1:]:
				tmp.append(float(ip_value))
			ip_rates.append(tmp)
	
	return (ir_time, ir_rates, ip_rates)





#  GLOBAL DATA FOR CAMX CB4 Mechanism 3 ...

SPC_Names = [
  'NO  ', 'NO2 ', 'O3  ', 'OLE ', 'PAN ', 'N2O5', 'PAR ', 'TOL ', 'XYL ',\
  'FORM', 'ALD2', 'ETH ', 'CRES', 'MGLY', 'OPEN', 'PNA ', 'CO  ', 'HONO',\
  'H2O2', 'HNO3', 'ISOP', 'MEOH', 'ETOH', 'CH4 ', 'O   ', 'OH  ', 'HO2 ',\
  'NO3 ', 'C2O3', 'XO2 ', 'XO2N', 'NTR ', 'CRO ', 'ISPD', 'TO2 ', 'ROR ',\
  'SO2 ', 'xHO2', 'H2O ', '-OOX', 'VOC', 'O3_old' ]

# Use pHO2 to track 'prompt HO2' and separate it from direct HO2.
# Use H2O as a product in some reactions to track a term pathway
# Use tracking species -OOX to collect the XO2+XO2 type reaction
#    products;  this is not a real species in CB4.

# only the first 24 species are included in the physical processes...
Max_Process_SPC = 24

# integer indexing for CB4 species
iNO   =  0
iNO2  =  1
iO3   =  2
iOLE  =  3
iPAN  =  4
iN2O5 =  5
iPAR  =  6
iTOL  =  7
iXYL  =  8
iFORM =  9
iALD2 = 10 
iETH  = 11 
iCRES = 12 
iMGLY = 13 
iOPEN = 14 
iPNA  = 15 
iCO   = 16
iHONO = 17 
iH2O2 = 18 
iHNO3 = 19 
iISOP = 20 
iMEOH = 21 
iETOH = 22 
iCH4  = 23
iO    = 24
iOH   = 25
iHO2  = 26
iNO3  = 27
iC2O3 = 28
iXO2  = 29
iXO2N = 30
iNTR  = 31
iCRO  = 32
iISPD = 33
iTO2  = 34
iROR  = 35
iSO2  = 36
ixHO2 = 37
iH2O  = 38
iOOX  = 39
iVOC  = 40
iO3_old = 41

max_i_spc = iO3_old + 1

# integer indexing for processes
iInitial    =  0
iChemistry  =  1
iEmiss_Area =  2
iEmiss_Pnt  =  3
iEmiss_PiG  =  4
iAdv_W      =  5
iAdv_E      =  6
iAdv_S      =  7
iAdv_N      =  8
iAdv_B      =  9
iAdv_T      = 10 
iDil_V      = 11 
iDif_W      = 12 
iDif_E      = 13 
iDif_S      = 14 
iDif_N      = 15 
iDif_B      = 16
iDif_T      = 17 
iDep_D      = 18 
iDep_W      = 19 
iChem_Aero  = 20 
iDilut      = 21 
iTrain      = 22 
iFinal      = 23 
iUn_Con     = 24
iAv_Cel_Vol = 25

max_i_process = iAv_Cel_Vol + 1



# Set up data representation for the net reactions

# Start with info about the net_rxns sets...
#   allocate vector of names of the net_rxns
net_rxn_names = []
#   this is indexed by next_net_rxn_set

# allocate a look up table between global species ID num (called iSPC)
#    and net_rxn_masses and net_rxn_spcname vector element 
#   [where SPC is any species name in the above global list].
net_rxn_species = []
#   this is indexed by next_net_rxn_set
#   each row of net_rxn_species will be a vector of len(max_i_spc)
#    that has in its iSPC position, the jindex number for the species
#    that are involved in this particular net reaction; iSPC that are not
#    involved in this net_rxn will have a jNONE value inserted

# the corresponding jindex number for a spc that is NOT in a net_rxn set
#     is negative to indicate that this iSPC is not used in this net_rxn
jNONE = -1

# allocate a vector of starting jindex in net_rxn_masses for each net_rxn set
net_rxn_jindex = []
#   these are also indexed by next_net_rxn_set


# the number for the next net_rxn set to be added...
next_net_rxn_set = 0
num_net_rxn_sets = 0


# allocate net rxn masses vector
net_rxn_masses = []
#   adjacent elements in this numeric vector hold and accumulate
#    the net masses for species iSPC in the set of rxns being lumped

daily_net_rxn_masses = []
#   accumulate the net_rxn_masses over all time.

# allocate net rxn species names index vector
net_rxn_spcname = []
#   adjacent elements in this integer vector hold the iSPC value
#    for the net masses in net_rxn_masses

# next location in net_rxn_masses vector to append another net rxn set
#   and in the net_rxn_spcname vector to insert the iSPC number
jstart_next_nr_set = 0

# collect each hour's net_rxn_masses in a list...
hour_number = []
hourly_net_rxn_masses = []


# create a function to return the j-index associated with species i in net 
#    reaction set k
def i2j(k,i):
	j = net_rxn_species[k][i]
	# print "i2j  k i  j", (k, i, j)
	return j
	

# Next add info about the net_processes sets...
#   allocate vector of names of the net_processes
net_processes_sets_names = []
#   this is indexed by next_net_proc_set

# allocate net processes masses array
net_process_masses = []

#   accumulate the net_process_masses over all time.
total_net_processes_masses = []

# the number for the next net_process set to be added...
next_net_proc_set = 0
num_net_proc_sets = 0

##### working here.....





# @@@@@@@@ N E T   R E A C T I O N S  @@@@@@@@@
#
#  1) ** O3+hv   Radical Source  **   IdNum = n_O3hvrad
#  2) ** HONO+hv Radical Source  **   IdNum = n_HONOhvrad
#  3) ** Ald+hv  Radical Source  **   IdNum = n_Aldhvrad
#  4) ** Ox+Org  Radical Source  **   IdNum = n_OxOrgrad
#  5) ** NO3+Org Radical Source  **   IdNum = n_NO3Orgrad
#  6) ** OH  + (organic+NO2)     **   IdNum = n_OHOrgOxid
#  7) ** C2O3     + NO Oxidation **   IdNum = n_C2O3NOOxid
#  8) ** XO2/XO2N + NO Oxidation **   IdNum = n_XO2NOOxid
#  9) ** PAN Production          **   IdNum = n_PANProd
# 10) ** HO2 to OH via Radical   **   IdNum = n_HO2toOHrad
# 11) ** HO2      + NO Oxidation **   IdNum = n_HO2NOOxid


## ==========================================================
# define the species storage for each set of net reactions

# init a counter for number of net reaction sets so far...
#    used to assign a name that refers to a given net_reaction.
nr_num = 0


## begin net_reaction storage allocation

#-------------------------------------------------------
### The Net Rxns for ** O3+hv Radical Source **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('O3+hv radical source')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iOH  ] = jj; jj += 1; indx_net_rxn.append(iOH  )


# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------
### The Net Rxns for ** HONO+hv Radical Source **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('HONO+hv radical source')
n_HONOhvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iOH  ] = jj; jj += 1; indx_net_rxn.append(iOH  )
jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species 
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------
### The Net Rxns for ** Ald+hv Radical Source **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('Ald+hv radical source')
n_Aldhvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc
indx_net_rxn = []


# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iISPD] = jj; jj += 1; indx_net_rxn.append(iISPD)
jndx_net_rxn[iHO2 ] = jj; jj += 1; indx_net_rxn.append(iHO2 )
jndx_net_rxn[iC2O3] = jj; jj += 1; indx_net_rxn.append(iC2O3)
jndx_net_rxn[iXO2 ] = jj; jj += 1; indx_net_rxn.append(iXO2 )
jndx_net_rxn[ixHO2] = jj; jj += 1; indx_net_rxn.append(ixHO2)

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species 
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------
### The Net Rxns for ** Ox+Org Radical Source **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('Ox+organic radical source')
n_OxOrgrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iO   ] = jj; jj += 1; indx_net_rxn.append(iO   )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iISPD] = jj; jj += 1; indx_net_rxn.append(iISPD)
jndx_net_rxn[iOH  ] = jj; jj += 1; indx_net_rxn.append(iOH  )
jndx_net_rxn[iHO2 ] = jj; jj += 1; indx_net_rxn.append(iHO2 )
jndx_net_rxn[iC2O3] = jj; jj += 1; indx_net_rxn.append(iC2O3)
jndx_net_rxn[iXO2 ] = jj; jj += 1; indx_net_rxn.append(iXO2 )
jndx_net_rxn[ixHO2] = jj; jj += 1; indx_net_rxn.append(ixHO2)
jndx_net_rxn[iXO2N] = jj; jj += 1; indx_net_rxn.append(iXO2N)

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species in the 
#   OH+organic net rxn
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------
### The Net Rxns for ** NO3+Org Radical Source **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('NO3+organic radical source')
n_NO3Orgrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)


# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc
indx_net_rxn = []


# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO3 ] = jj; jj += 1; indx_net_rxn.append(iNO3 )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iISPD] = jj; jj += 1; indx_net_rxn.append(iISPD)
jndx_net_rxn[iHO2 ] = jj; jj += 1; indx_net_rxn.append(iHO2 )
jndx_net_rxn[iC2O3] = jj; jj += 1; indx_net_rxn.append(iC2O3)
jndx_net_rxn[iXO2 ] = jj; jj += 1; indx_net_rxn.append(iXO2 )
jndx_net_rxn[ixHO2] = jj; jj += 1; indx_net_rxn.append(ixHO2)
jndx_net_rxn[iXO2N] = jj; jj += 1; indx_net_rxn.append(iXO2N)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iNTR ] = jj; jj += 1; indx_net_rxn.append(iNTR )

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species 
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------
### The Net Rxns for ** OH + (organic+NO2) **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('OH+(organic+NO2)')
n_OHOrgOxid = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)


# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc
indx_net_rxn = []


# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
# use H2O to catch mag of OH + HO2 = H2O + O2 reaction..
jndx_net_rxn[iH2O ] = jj; jj += 1; indx_net_rxn.append(iH2O )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH)
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iCH4 ] = jj; jj += 1; indx_net_rxn.append(iCH4 )
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iISPD] = jj; jj += 1; indx_net_rxn.append(iISPD)
jndx_net_rxn[iOH  ] = jj; jj += 1; indx_net_rxn.append(iOH  )
jndx_net_rxn[iHO2 ] = jj; jj += 1; indx_net_rxn.append(iHO2 )
jndx_net_rxn[iC2O3] = jj; jj += 1; indx_net_rxn.append(iC2O3)
jndx_net_rxn[iXO2 ] = jj; jj += 1; indx_net_rxn.append(iXO2 )
jndx_net_rxn[ixHO2] = jj; jj += 1; indx_net_rxn.append(ixHO2)
jndx_net_rxn[iXO2N] = jj; jj += 1; indx_net_rxn.append(iXO2N)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iVOC ] = jj; jj += 1; indx_net_rxn.append(iVOC )

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species in the 
#   OH+organic net rxn
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	

#-------------------------------------------------------
# The Net Rxns for ** C2O3 + NO Oxidation **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('C2O3 + NO Oxidation')
n_C2O3NOOxid = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)


# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart
jndx_net_rxn[iC2O3] = jj; jj += 1; indx_net_rxn.append(iC2O3)
jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iXO2 ] = jj; jj += 1; indx_net_rxn.append(iXO2 )
jndx_net_rxn[ixHO2] = jj; jj += 1; indx_net_rxn.append(ixHO2)
jndx_net_rxn[iOOX ] = jj; jj += 1; indx_net_rxn.append(iOOX )

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species 
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------
# The Net Rxns for ** XO2/XO2N + NO Oxidation **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('XO2/XO2N Radical NO Oxidation')
n_XO2NOOxid = nr_num  # provide an interger with a name for this nr
nr_num += 1


this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)


# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart
jndx_net_rxn[iXO2 ] = jj; jj += 1; indx_net_rxn.append(iXO2 )
jndx_net_rxn[iXO2N] = jj; jj += 1; indx_net_rxn.append(iXO2N)
jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iNTR ] = jj; jj += 1; indx_net_rxn.append(iNTR )
jndx_net_rxn[iOOX ] = jj; jj += 1; indx_net_rxn.append(iOOX )
jndx_net_rxn[ixHO2] = jj; jj += 1; indx_net_rxn.append(ixHO2)

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species 
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------
### The Net Rxns for ** PAN Production **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('PAN Production')
n_PANProd = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)


# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart
jndx_net_rxn[iC2O3] = jj; jj += 1; indx_net_rxn.append(iC2O3)
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN )

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species in the 
#   OH+organic net rxn
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)

#-------------------------------------------------------
### The Net Rxns for ** HO2 to OH via radical **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('HO2 to OH via Radical')
n_HO2toOHrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)


# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc
indx_net_rxn = []


# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iC2O3] = jj; jj += 1; indx_net_rxn.append(iC2O3)
jndx_net_rxn[iXO2 ] = jj; jj += 1; indx_net_rxn.append(iXO2 )
jndx_net_rxn[ixHO2] = jj; jj += 1; indx_net_rxn.append(ixHO2)
jndx_net_rxn[iXO2N] = jj; jj += 1; indx_net_rxn.append(iXO2N)
jndx_net_rxn[iHO2 ] = jj; jj += 1; indx_net_rxn.append(iHO2 )
jndx_net_rxn[iOH  ] = jj; jj += 1; indx_net_rxn.append(iOH  )
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iOOX ] = jj; jj += 1; indx_net_rxn.append(iOOX )

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species in the 
#   OH+organic net rxn
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------
### The Net Rxns for ** HO2+NO Oxidation **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('HO2 + NO Oxidation')
n_HO2NOOxid = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)


# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc
indx_net_rxn = []


# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart
jndx_net_rxn[iHO2 ] = jj; jj += 1; indx_net_rxn.append(iHO2 )
jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iOH  ] = jj; jj += 1; indx_net_rxn.append(iOH  )
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species in the 
#   OH+organic net rxn
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)

#-------------------------------------------------------
### The Net Rxns for ** NO cycle **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('NO CYCLE')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3_old] = jj; jj += 1; indx_net_rxn.append(iO3_old)
jndx_net_rxn[iNO3 ] = jj; jj += 1; indx_net_rxn.append(iNO3 )
jndx_net_rxn[iOH  ] = jj; jj += 1; indx_net_rxn.append(iOH  )
jndx_net_rxn[iHO2 ] = jj; jj += 1; indx_net_rxn.append(iHO2 )
jndx_net_rxn[iC2O3] = jj; jj += 1; indx_net_rxn.append(iC2O3)
jndx_net_rxn[iXO2 ] = jj; jj += 1; indx_net_rxn.append(iXO2 )
jndx_net_rxn[iTO2 ] = jj; jj += 1; indx_net_rxn.append(iTO2 )
jndx_net_rxn[iNTR ] = jj; jj += 1; indx_net_rxn.append(iNTR )
jndx_net_rxn[iXO2N] = jj; jj += 1; indx_net_rxn.append(iXO2N)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------

#-------------------------------------------------------
### The Net Rxns for ** NO2 cycle **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('NO2 CYCLE')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3) #### didn't include water...
jndx_net_rxn[iNO3 ] = jj; jj += 1; indx_net_rxn.append(iNO3 )
jndx_net_rxn[iOH  ] = jj; jj += 1; indx_net_rxn.append(iOH  )
jndx_net_rxn[iC2O3] = jj; jj += 1; indx_net_rxn.append(iC2O3)
jndx_net_rxn[iTO2 ] = jj; jj += 1; indx_net_rxn.append(iTO2 )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN )
jndx_net_rxn[iCRO ] = jj; jj += 1; indx_net_rxn.append(iCRO )
jndx_net_rxn[iNTR ] = jj; jj += 1; indx_net_rxn.append(iNTR )
jndx_net_rxn[iROR ] = jj; jj += 1; indx_net_rxn.append(iROR )
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iISPD] = jj; jj += 1; indx_net_rxn.append(iISPD)


# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------
################ PROCESSES ##############################

#IN. CONCENTRATIONS:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('INITIAL CONCENTRATIONS')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#CHEMISTRY:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('CHEMISTRY')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	

# AREA EMISSIONS:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('AREA EMISSIONS')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	
	# POINT EMISSIONS:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('POINT EMISSIONS')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# TOTAL EMISSIONS:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('TOTAL EMISSIONS')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	
	# PIG CHANGE:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('PiG CHANGE')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	
	
	
	# WEST B. ADVECTION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('WEST B. ADVECTION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	
	# EAST B. ADVECTION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('EAST B. ADVECTION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	
	# SOUTH B. ADVECTION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('SOUTH B. ADVECTION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	
	# NORTH B. ADVECTION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('NORTH B. ADVECTION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	
	# BOTTOM B. ADVETION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('BOTTOM B. ADVECTION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	
	# TOP B. ADVECTION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('TOP B. ADVECTION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	
	# DIL. IN THE VERT.:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('DILUTION IN THE VERTICAL')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	
	# WEST B. DIFFUSION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('WEST B. DIFFUSION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# EAST B. DIFFUSION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('EAST B. DIFFUSION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# SOUTH B. DIFFUSION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('SOUTH B. DIFFUSION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# NORTH B. DIFFUSION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('NORTH B. DIFFUSION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# BOTTOM B. DIFFUSION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('BOTTOM B. DIFFUSION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# TOP B. DIFFUSION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('TOP B. DIFFUSION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# DRY DEPOSITION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('DRY DEPOSITION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# WET DEPOSITION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('WET DEPOSITION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# AEROSOL CHEMISTRY:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('AEROSOL CHEMISTRY')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# DILUTION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('DILUTION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# EN(DE)TRAINMENT:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('EN(DE)TRAINMENT')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# HORIZONTAL TRANS:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('HORIZONTAL TRANS')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# VERTICAL TRANS:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('VERTICAL TRANS')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# FINAL CONCENTRATION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('FINAL CONCENTRATION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	
	# UNITS CONVERSION:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('UNITS CONVERSION')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	# AVERAGE CELL VOLUME:

this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('AVERAGE CELL VOLUME')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc + [jNONE]*max_i_process
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN ) 
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iPNA ] = jj; jj += 1; indx_net_rxn.append(iPNA )
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHONO] = jj; jj += 1; indx_net_rxn.append(iHONO)
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH) 
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)



# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
	




## >>>>>> finished setting up all net reaction storage <<<<<<
num_net_rxn_sets  = this_net_rxn_set
max_j  = jstart_next_nr_set


## >>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
##  P A R T   T W O  : C A L C U L A T E  T H E  N E T  M A S S E S

## ==========================================================
if options.verbose:
	print "Opening the files.."

# open the files...
fout = open(output_filename, 'w')

fin = open(input_filename, 'r')

# first two lines are 'doc' lines giving data source
doc1 = fin.readline()
line = fin.readline()

if options.verbose:
	print "IRR file doc line was"
	print doc1
	print


## read first hours' data... 
(time, ir, ip) = get_irr_data(fin);  


while ( time > 0 ) :
	if options.verbose:
		print "Time %d  hours" % time
	if len(ir) != 97:
		print "ERROR: no ir data!"
		sys.exit(1)
		
	# ... process the ir data for this time...
	#      ... zero out the net reaction masses vector
	for i in range(0,len(net_rxn_masses)):
		net_rxn_masses[i] = 0.0
	
	# initialize the counter for net reaction sets so
	# that it can be incremented in each reaction...
	kk = -1
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** O3+hv Radical Source **
	kk += 1  
	
	# new OH from O3+hv
	# { 11} O1D+H2O =2*OH
	# { 34} H2O2    =2*OH
	

	net_rxn_masses[i2j(kk,iO3  )] = -ir[11]
	net_rxn_masses[i2j(kk,iH2O2)] = -ir[34]
	net_rxn_masses[i2j(kk,iOH  )] = 2.0*ir[11] + 2.0*ir[34]
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** HONO+hv Radical Source **
	kk += 1  
	
	# new OH from HONO+hv
	# { 21} NO+NO2+H2O=2.000*HONO
	# { 22} NO+OH=HONO
	# { 23} HONO=NO+OH
	# { 24} OH+HONO=NO2
	# { 25} HONO+HONO=NO+NO2
	
	# reactant losses in HONO +hv...
	#   ... the NO, NO2 losses
	net_rxn_masses[i2j(kk,iNO  )] = -ir[21]-ir[22]-ir[59]-ir[78]-ir[94]\
										+ir[23]+ir[25]
	net_rxn_masses[i2j(kk,iNO2 )] = -ir[21]+ir[24]+ir[25]
		
	#   ... the HONO products
	net_rxn_masses[i2j(kk,iHONO)] =  2.0*ir[21]+ir[22]-ir[23]-ir[24]-2*ir[25]
	net_rxn_masses[i2j(kk,iOH  )] =  ir[23]-ir[22]-ir[24]
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** Ald+hv Radical Source **
	kk += 1  
	
	# { 38} FORM=2*HO2+CO
	# { 45} ALD2=FORM+HO2+CO + xHO2 + XO2
	# { 69} OPEN=C2O3+HO2+CO
	# { 74} MGLY=C2O3+HO2+CO
	# { 95} ISPD=0.333*CO+0.067*ALD2+0.9*FORM+0.832*PAR+1.033*HO2+0.7*XO2+0.967*C2O3
	
	# NB: in {95} the 1.033 HO2 and  0.7 XO2 are treated as 
	#                 0.333 HO2 and (0.7 xHO2 + 0.7 XO2)
	
	# reactant losses in ALD +hv...
	#   ... the alds
	net_rxn_masses[i2j(kk,iFORM)] = -ir[38]
	net_rxn_masses[i2j(kk,iALD2)] = -ir[45]
	net_rxn_masses[i2j(kk,iOPEN)] = -ir[69]
	net_rxn_masses[i2j(kk,iMGLY)] = -ir[74]
	net_rxn_masses[i2j(kk,iISPD)] = -ir[95]
	
	#   ... the radical products
	net_rxn_masses[i2j(kk,iHO2 )] =  2*ir[38]+ir[45]+ir[69]+ir[74]+0.333*ir[95]
	net_rxn_masses[i2j(kk,iC2O3)] =  ir[69]+ir[74]+0.967*ir[95]
	net_rxn_masses[i2j(kk,iXO2 )] =  ir[45]+0.7*ir[95]
	net_rxn_masses[i2j(kk,ixHO2)] =  ir[45]+0.7*ir[95]
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** Ox+Org Radical Source **
	kk += 1  
	
	# new OH, XO2, XO2N, and HO2 sources from O+org
	# { 40} FORM+O=OH+HO2+CO
	# { 42} ALD2+O=C2O3+OH
	# { 56} O+OLE=0.63*ALD2+0.38*HO2+0.28*XO2+0.3*CO+0.2*FORM+0.02*XO2N+0.22*PAR+0.2*OH
	#                      [0.10*HO2+0.28*xHO2+0.28*XO2]
	# { 60} O+ETH=FORM+1.7*HO2+CO+0.7*XO2+0.3*OH
	#                 [1.0*HO2   +0.7*xHO2+0.7*XO2]
	# { 75} O+ISOP=0.75*ISPD+0.5*FORM+0.25*XO2+0.25*HO2+0.25*C2O3+0.25*PAR
	#                 [0.25*xHO2+0.25*XO2]
	# 
	# new OH, XO2, XO2N, and HO2 sources from O3+org
	# { 58} O3+OLE=0.5*ALD2+0.74*FORM+0.44*HO2+0.22*XO2+0.1*OH+0.33*CO+-1*PAR
	#                   [0.22*HO2 + 0.22*xHO2 + 0.22*XO2]
	# { 62} O3+ETH=FORM+0.42*CO+0.12*HO2
	# { 71} OPEN+O3=0.03*ALD2+0.62*C2O3+0.7*FORM+0.76*HO2+0.03*XO2+0.69*CO+0.08*OH+0.2*MGLY
	#                   [0.73*HO2 + 0.03*xHO2 + 0.03*XO2]
	#                   
	# { 77} O3+ISOP=0.65*ISPD+0.6*FORM+0.066*HO2+0.2*XO2+
	#               0.266*OH+0.2*C2O3+0.15*ALD2+0.35*PAR+0.066*CO
	#                   [0.00*HO2 + 0.066*xHO2 + 0.2*XO2]
	# { 93} O3+ISPD=0.114*C2O3+0.15*FORM+0.154*HO2+0.064*XO2+
	#               0.85*MGLY+0.268*OH+0.02*ALD2+0.36*PAR+0.225*CO
	#                   [0.090*HO2 + 0.064*xHO2 + 0.064*XO2]
	
	# reactant losses in Ox + org...
	#   ... the O and O3 losses
	net_rxn_masses[i2j(kk,iO   )] = -ir[40]-ir[42]-ir[56]-ir[60]-ir[75]
	net_rxn_masses[i2j(kk,iO3  )] = -ir[58]-ir[62]-ir[71]-ir[77]-ir[93]
	
	#   ... the organics losses
	net_rxn_masses[i2j(kk,iOLE )] = -ir[56]-ir[58]
	net_rxn_masses[i2j(kk,iFORM)] = -ir[40]
	net_rxn_masses[i2j(kk,iALD2)] = -ir[42]
	net_rxn_masses[i2j(kk,iETH )] = -ir[60]-ir[62]
	net_rxn_masses[i2j(kk,iOPEN)] = -ir[71]
	net_rxn_masses[i2j(kk,iISOP)] = -ir[75]-ir[77]
	net_rxn_masses[i2j(kk,iISPD)] = -ir[93]
	
	#   ... the new radical products
	net_rxn_masses[i2j(kk,iOH  )] =  ir[40]+ir[42]+0.20*ir[56]+0.30*ir[60]\
										+0.1*ir[58]+0.08*ir[71]+0.266*ir[77]\
										+0.268*ir[93]
	net_rxn_masses[i2j(kk,iHO2 )] =  ir[40]+0.10*ir[56]+1.0*ir[60]\
										+0.22*ir[58]+0.12*ir[62]+0.73*ir[71]\
										+0.090*ir[93]
	net_rxn_masses[i2j(kk,iC2O3)] =  ir[42]+0.25*ir[75]+0.62*ir[71]+0.2*ir[77]\
										+0.114*ir[93]
	net_rxn_masses[i2j(kk,iXO2 )] =  0.28*ir[56]+0.7*ir[60]+0.25*ir[75]\
										+0.22*ir[58]+0.03*ir[71]+0.066*ir[77]\
										+0.064*ir[93]
	net_rxn_masses[i2j(kk,ixHO2)] =  0.28*ir[56]+0.7*ir[60]+0.25*ir[75]\
										+0.22*ir[58]+0.03*ir[71]+0.066*ir[77]\
										+0.066*ir[93]
	net_rxn_masses[i2j(kk,iXO2N)] =  0.02*ir[56]
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** NO3+Org Radical Source **
	kk += 1  
	
	# new OH, XO2, HO2 NTR, and HNO3 sources from NO3+org
	# { 41} FORM+NO3=HNO3+HO2+CO
	# { 44} ALD2+NO3=C2O3+HNO3
	# { 59} NO3+OLE=0.91*XO2+FORM+0.09*XO2N+ALD2+NO2+-1*PAR
	#               [ 0.0*xHO2 + 0.91*XO2]
	# { 78} NO3+ISOP=0.2*ISPD+0.8*NTR+XO2+0.8*HO2+0.2*NO2+0.8*ALD2+2.4*PAR
	#                [ 0.0 *HO2 + 0.8*xHO2 + 1.0*XO2]
	# { 94} NO3+ISPD=0.357*ALD2+0.282*FORM+0.925*HO2+0.075*XO2+
	#                1.282*PAR+0.643*CO+0.85*NTR+0.075*C2O3+0.15*HNO3
	#                [ 0.725*HO2 + 0.075*xHO2 + 0.075*XO2]
	
	# reactant losses in NO3 + org...
	#   ... the NO3 losses
	net_rxn_masses[i2j(kk,iNO3 )] = -ir[41]-ir[44]-ir[59]-ir[78]-ir[94]
	
	#   ... the organics losses
	net_rxn_masses[i2j(kk,iOLE )] = -ir[59]
	net_rxn_masses[i2j(kk,iFORM)] = -ir[41]
	net_rxn_masses[i2j(kk,iALD2)] = -ir[44]
	net_rxn_masses[i2j(kk,iISOP)] = -ir[78]
	net_rxn_masses[i2j(kk,iISPD)] = -ir[94]
	
	#   ... the new radical products
	net_rxn_masses[i2j(kk,iHO2 )] =  ir[41]+0.725*ir[94]
	net_rxn_masses[i2j(kk,iC2O3)] =  ir[44]+0.75*ir[75]
	net_rxn_masses[i2j(kk,iXO2 )] =  0.91*ir[59]+ir[78]+0.075*ir[94]
	net_rxn_masses[i2j(kk,ixHO2)] =  0.08*ir[78]+0.075*ir[94]
	net_rxn_masses[i2j(kk,iXO2N)] =  0.09*ir[59]
	net_rxn_masses[i2j(kk,iHNO3)] =  ir[41]+ir[44]+ir[94]
	net_rxn_masses[i2j(kk,iNO2 )] =  ir[59]+0.2*ir[78]
	net_rxn_masses[i2j(kk,iNTR )] =  0.80*ir[78]+0.85*ir[94]
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** OH + (organic+NO2) **
	kk += 1  
	
	# { 26} OH + NO2  = HNO3
	
	# { 90} OH + HO2  = H2O + O2
	
	# { 12} OH + O3   = HO2
	# { 35} OH + H2O2 = HO2
	# { 36} OH + CO   = HO2
	# { 37} OH + FORM = HO2  + CO
	# { 84} OH + MEOH = HO2  + FORM
	# { 85} OH + ETOH = HO2  + ALD2
	
	# { 43} OH + ALD2 = C2O3
	# { 73} OH + MGLY = C2O3 + XO2 
	# { 70} OH + OPEN = C2O3 + XO2 +xHO2 +2*CO +HO2 +FORM  
	#
	# { 51} OH + CH4  = FORM + XO2 +xHO2
	# { 52} OH + PAR  = 0.87*XO2 + 0.13*XO2N + (0.11*xHO2 + 0.76*ROR) +0.11*ALD2+-0.11*PAR
	#      { 53} ROR  = 0.96*XO2 + 0.04*XO2N +  0.94*xHO2 + 1.1*ALD2 + -2.1*PAR 
	#                                    
	#      { 54} ROR  = HO2  <direct HO2, not prompt>
	#      { 55} ROR + NO2 = NTR    ! OMITTED from this net reaction set
	#
	# { 61} OH + ETH  = XO2+1.56*FORM+0.22*ALD2 +xHO2
	# { 57} OH + OLE  = FORM+ALD2+-1*PAR+XO2 +xHO2
	#
	# { 63} OH + TOL  = 0.08*XO2+0.08*xHO2+ 0.36*HO2+0.36*CRES+0.56*TO2
	#      { 64} TO2 + NO = 0.9*NO2+0.9*xHO2+0.9*OPEN+0.1*NTR
	#                    [ code 0.9*NO2 as 0.9*XO2 + 0.9*xHO2 ]
	#                    [ code 0.1*NTR as 0.1*XNO2           ]
	#                    [ add NO loss and NO2 gain in XO2/XO2N net reactions]
	#      { 65} TO2  = CRES + HO2
	# { 66} OH + CRES = 0.4*CRO + 0.6*XO2 + 0.6*xHO2 + 0.3*OPEN
	#      { 68} CRO + NO2 = NTR
	# { 72} OH + XYL  = 0.5*XO2 + 0.5*xHO2 +0.2*HO2 + 0.2*CRES +0.8*MGLY+1.1*PAR+0.3*TO2
	#
	# { 76} OH + ISOP = 0.912*ISPD+0.629*FORM+0.991*XO2+0.912*xHO2+0.088*XO2N
	#                                                   
	# { 92} OH + ISPD = 1.565*PAR+0.167*FORM + 0.713*XO2 + 0.503*xHO2 + 
	#                     0.334*CO+0.168*MGLY+0.273*ALD2+0.498*C2O3
	
	# notes:: 
	#   3) for ir[55], ROR, NTR is from NO2+rad not NO+RO2, omitted here
	#                    but tracked in '?????' net reaction set.
	#   4) for ir[66], CRO, NTR is from NO2+rad not NO+RO2, omitted the CRO production.
	#                   but ir[68] production is tracked in '?????' net reaction set.
	#
	#   5) for ir[64], TO2, the 0.9*NO2 is coded as 0.9*XO2 with 0.9*xHO2
	#                       and the 0.1*NTR is coded as 0.1*XO2N
	#                       The NO is omitted as a a loss here but coded as
	#                       loss in the XO2 + NO net reaction below.
	
	
	#   ... the organics losses
	net_rxn_masses[i2j(kk,iNO2 )] = -ir[26]
	net_rxn_masses[i2j(kk,iH2O )] = -ir[90]
	net_rxn_masses[i2j(kk,iO3  )] = -ir[12]
	net_rxn_masses[i2j(kk,iH2O2)] = -ir[35]
	net_rxn_masses[i2j(kk,iCO  )] = -ir[36]
	net_rxn_masses[i2j(kk,iFORM)] = -ir[37]
	net_rxn_masses[i2j(kk,iMEOH)] = -ir[84]
	net_rxn_masses[i2j(kk,iETOH)] = -ir[85]
	net_rxn_masses[i2j(kk,iALD2)] = -ir[43]
	net_rxn_masses[i2j(kk,iMGLY)] = -ir[73]
	net_rxn_masses[i2j(kk,iCH4 )] = -ir[51]
	net_rxn_masses[i2j(kk,iPAR )] = -ir[52]
	net_rxn_masses[i2j(kk,iETH )] = -ir[61]
	net_rxn_masses[i2j(kk,iOLE )] = -ir[57]
	net_rxn_masses[i2j(kk,iTOL )] = -ir[63]
	net_rxn_masses[i2j(kk,iXYL )] = -ir[72]
	net_rxn_masses[i2j(kk,iCRES)] = -ir[66]
	net_rxn_masses[i2j(kk,iOPEN)] = -ir[70]
	net_rxn_masses[i2j(kk,iISOP)] = -ir[76]
	net_rxn_masses[i2j(kk,iISPD)] = -ir[92]
	#    ... the OH losses
	net_rxn_masses[i2j(kk,iOH  )] = -ir[26]-ir[90]-ir[12]-ir[35]-ir[36]\
	                                -ir[37]-ir[84]-ir[85]-ir[43]-ir[73]\
	                                -ir[70]-ir[51]-ir[52]-ir[61]-ir[57]\
	                                -ir[63]-ir[66]-ir[72]-ir[76]-ir[92]
	                                
	#    ... the radical products...
	net_rxn_masses[i2j(kk,iHO2 )] = +ir[12]+ir[35]+ir[36]+ir[37]\
	                                +ir[84]+ir[85]+ir[70]\
	                                +ir[54]+0.36*ir[63]\
	                                +0.2*ir[72]+ir[65]
	
	net_rxn_masses[i2j(kk,iC2O3)] =  ir[43]+ir[70]+ir[73]+0.498*ir[92]
	
	net_rxn_masses[i2j(kk,iXO2 )] =  ir[51]+0.87*ir[52]+ir[57]+ir[61]\
									+0.08*ir[63]+0.6*ir[66]\
									+ir[70]+0.5*ir[72]+ir[73]+0.991*ir[76]\
									+0.713*ir[92]+0.96*ir[53]+0.9*ir[64]
									
	net_rxn_masses[i2j(kk,ixHO2)] =  ir[70]+ir[51]+0.11*ir[52]+0.94*ir[53]\
	                                +ir[61]+ir[57]\
									+0.08*ir[63]+0.9*ir[64]+0.6*ir[66]\
									+0.5*ir[72]+0.912*ir[76]\
									+0.503*ir[92]
									
	#     ... radical losses via NO2 or NO->NO2 reactions ...
	net_rxn_masses[i2j(kk,iHNO3)] =  ir[26]
	net_rxn_masses[i2j(kk,iXO2N)] =  0.13*ir[52]+0.088*ir[76]+0.04*ir[53]\
										+0.1*ir[64]
	
	# sum of VOCs 
	net_rxn_masses[i2j(kk,iVOC )] =  -ir[37]-ir[84]-ir[85]-ir[43]-ir[73]\
                                     -ir[51]-ir[52]-ir[61]-ir[57]-ir[63]\
                                     -ir[72]-ir[66]-ir[70]-ir[76]-ir[92]


	
	#-------------------------------------------------------
	# The Net Rxns for ** C2O3 + NO Oxidation **
	kk += 1  
		
	# {46} C2O3 + NO = FORM + NO2 + xHO2 + XO2
	#
	# {49} C2O3 + C2O3 = 2FORM + 2XO2 + 2xHO2
	# {50} C2O3 + HO2 = 0.79FORM + 0.79XO2 + 0.79xHO2 +0.79OH
	#
	# notes:
	#    1) show the ir[50} OH production as a RO2 + HO2 --> OH 
	#       process in the 'HO2 to OH via radical' reaction set below.
	
	net_rxn_masses[i2j(kk,iC2O3)] = -ir[46] -2*ir[49] -ir[50]
	net_rxn_masses[i2j(kk,iNO  )] = -ir[46]
	net_rxn_masses[i2j(kk,iNO2 )] = +ir[46]
	net_rxn_masses[i2j(kk,iFORM)] = +ir[46] +2*ir[49] +0.79*ir[50]
	net_rxn_masses[i2j(kk,iXO2 )] = +ir[46] +2*ir[49] +0.79*ir[50]
	net_rxn_masses[i2j(kk,ixHO2)] = +ir[46] +2*ir[49] +0.79*ir[50]
	net_rxn_masses[i2j(kk,iOOX )] = +0.21*ir[50]
		
	
	#-------------------------------------------------------
	# The Net Rxns for ** XO2/XO2N + NO Oxidation **
	kk += 1  
		
	# {79} XO2  + NO   = NO2
	
	# {81} XO2N + NO   =  NTR
	#
	# {80} XO2  + XO2  = 2*-OOX
	# {88} XO2N + XO2N = 2*-OOX
	# {89} XO2  + XO2N = 2*-OOX
	
	# notes:: 
	#   1) -OOX is an added 'peroxide-like' species to track XO2 termination
	#   2)  ir[64] see note 5 under ** OH + (organic+NO2) **
	
	net_rxn_masses[i2j(kk,iXO2 )] = -ir[79] -2*ir[80] -ir[89]
	net_rxn_masses[i2j(kk,iXO2N)] = -ir[81] -2*ir[88] -ir[89] +0.1*ir[64]
	net_rxn_masses[i2j(kk,iNO  )] = -ir[79]   -ir[81] -ir[64]
	net_rxn_masses[i2j(kk,iNO2 )] = +ir[79] +0.9*ir[64]
	
	net_rxn_masses[i2j(kk,ixHO2)] = ir[45]+0.7*ir[95]\
	                                +0.28*ir[56]+0.7*ir[60]+0.25*ir[75]\
	                                +0.22*ir[58]+0.03*ir[71]+0.066*ir[77]\
	                                +0.066*ir[93]\
	                                +0.08*ir[78]+0.075*ir[94]\
	                                +ir[70]+ir[51]+0.11*ir[52]+0.94*ir[53]\
	                                +ir[61]+ir[57]\
	                                +0.08*ir[63]+0.9*ir[64]+0.6*ir[66]\
	                                +0.5*ir[72]+0.912*ir[76]\
	                                +0.503*ir[92]\
	                                +ir[46] +2*ir[49] +0.79*ir[50]\
	                                +0.79*ir[50]
	
	net_rxn_masses[i2j(kk,iNTR )] = +ir[81]
	net_rxn_masses[i2j(kk,iOOX )] = +2*ir[80] +2*ir[88] +2*ir[89]
	

	
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** PAN Production **
	kk += 1
	
	
	# {47} C2O3 + NO2 = PAN
	# {48} PAN = C2O3 + NO2
	
	net_rxn_masses[i2j(kk,iC2O3)] = -ir[47] +ir[48]
	net_rxn_masses[i2j(kk,iNO2 )] = -ir[47] +ir[48]
	net_rxn_masses[i2j(kk,iPAN )] = +ir[47] -ir[48]
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** HO2 to OH via radical **
	kk += 1  
	
	# {13} O3   + HO2 = OH
	# {50} C2O3 + HO2 = 0.79 OH + 0.79 FORM + 0.79 XO2 + 0.79 xHO2

	# {86} XO2  + HO2  =   -OOX  
	# {87} XO2N + HO2  =   -OOX  
	
	# {32} HO2 + HO2 = H2O2
	# {33} HO2 + HO2 + H2O = H2O2
	
	
	
	net_rxn_masses[i2j(kk,iO3  )] = -ir[13]
	net_rxn_masses[i2j(kk,iC2O3)] = -ir[50]
	net_rxn_masses[i2j(kk,iXO2 )] = -ir[86]+0.79*ir[50]
	net_rxn_masses[i2j(kk,ixHO2)] = +0.79*ir[50]
	net_rxn_masses[i2j(kk,iXO2N)] = -ir[87]
	net_rxn_masses[i2j(kk,iHO2 )] = -ir[13] -ir[50] -ir[86] -ir[87]\
	                                   -2.0*ir[32] -2.0*ir[33]
	                                   
	net_rxn_masses[i2j(kk,iOH  )] = +ir[13] +0.79*ir[50]
	net_rxn_masses[i2j(kk,iH2O2)] = +ir[32] +ir[33] 
	net_rxn_masses[i2j(kk,iOOX )] = +ir[86] +ir[87]
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** HO2+NO Oxidation **
	kk += 1  
	
	# {28} HO2 + NO = OH + NO2
	
	# {29} HO2 + NO2 = PNA
	# {30} PNA = HO2 + NO2
	# {31} OH + PNA = NO2
	
	
	net_rxn_masses[i2j(kk,iHO2 )] = -ir[28] -ir[29] +ir[30]
	net_rxn_masses[i2j(kk,iNO  )] = -ir[28]
	net_rxn_masses[i2j(kk,iNO2 )] = +ir[28] -ir[29] +ir[30] +ir[31]
	net_rxn_masses[i2j(kk,iOH  )] = +ir[28]
	net_rxn_masses[i2j(kk,iPNA )] = +ir[29] -ir[30] -ir[31]
	
	## NO Cycle
	
	kk += 1  
	
	

	chemNO = +ir[1] +ir[4] +0.11*ir[14] +ir[16] -ir[3] -ir[6] -ir[15]
	chem_NO = 0
	if chemNO > 0:
		chem_NO = chemNO
	NO_O3 = -ir[1] -ir[4] -0.11*ir[14] -ir[16] +ir[3] +ir[6] +ir[15]
	old_O3 = 0
	if NO_O3 > 0:
		old_O3 = NO_O3
	from_O3_old = +ir[3] +ir[6] +0.89*ir[14] +2*ir[15] +ir[16] +ir[17] -ir[1] -ir[4] -ir[5] -ir[7] -ir[16] -ir[19]
	rem_O3_old = 0
	if from_O3_old > 0:
		rem_O3_old = from_O3_old
	NOz = old_O3 - rem_O3_old
	
	
	
	net_rxn_masses[i2j(kk,iNO  )] = +ip[iNO][iInitial] +ip[iNO][iEmiss_Area] \
									+ip[iNO][iEmiss_Pnt] +ip[iNO][iEmiss_PiG] +ip[iNO][iAdv_W]\
									+ip[iNO][iAdv_E] +ip[iNO][iAdv_S] +ip[iNO][iAdv_N]\
									+ip[iNO][iAdv_B] +ip[iNO][iAdv_T] +ip[iNO][iDil_V]\
									+ip[iNO][iDif_W] +ip[iNO][iDif_E] +ip[iNO][iDif_S]\
									+ip[iNO][iDif_N] +ip[iNO][iDif_B] +ip[iNO][iDif_T]\
									+ip[iNO][iDep_D] +ip[iNO][iDep_W] +ip[iNO][iChem_Aero]\
									+ip[iNO][iDilut] +ip[iNO][iTrain] -ip[iNO][iFinal]\
									+chem_NO
	net_rxn_masses[i2j(kk,iNO2 )] = +ir[28] +ir[46] +ir[64] +ir[79] +rem_O3_old
	net_rxn_masses[i2j(kk,iO3_old)] = old_O3 ## this is not actually ozone, but cosumption of NO in the O3/NO2 chemistry
	net_rxn_masses[i2j(kk,iNO3 )] = NOz    #+ir[5] +ir[7] -ir[14] -ir[16] +ir[17] -ir[19] 
	net_rxn_masses[i2j(kk,iHONO)] = +0.2*ir[21] +ir[22] -ir[23] -ir[24] -2*ir[25]
	net_rxn_masses[i2j(kk,iOH  )] = +ir[28]
	net_rxn_masses[i2j(kk,iHO2 )] = -ir[28]
	net_rxn_masses[i2j(kk,iC2O3)] = -ir[46]
	net_rxn_masses[i2j(kk,iXO2 )] = -ir[79]
	net_rxn_masses[i2j(kk,iTO2 )] = -ir[64]
	net_rxn_masses[i2j(kk,iNTR )] = +ir[81]
	net_rxn_masses[i2j(kk,iXO2N)] = -ir[81]
	
	
	############## NO2 Cycle ################################
	
	kk += 1  
	
	O3 = +ir[2] -ir[3] -ir[7] -ir[8] -ir[9]
	prodO3 = 0
	if O3 > 0:
		prodO3 = O3
	
	net_rxn_masses[i2j(kk,iNO  )] = chem_NO
	net_rxn_masses[i2j(kk,iNO2 )] = -ir[26] -ir[47] +ir[48] -ir[55] +0.9*ir[64]\
									-ir[68] +0.2*ir[78] -ir[96] -chem_NO
	net_rxn_masses[i2j(kk,iO3  )] = prodO3
	net_rxn_masses[i2j(kk,iHNO3)] = +ir[26]
	net_rxn_masses[i2j(kk,iNO3 )] = -ir[94]   ### JUST FOR THE PRODUCTION OF NTR
	net_rxn_masses[i2j(kk,iOH  )] = -ir[26] 
	net_rxn_masses[i2j(kk,iC2O3)] = -ir[47] +ir[48]
	net_rxn_masses[i2j(kk,iTO2 )] = -ir[64]
	net_rxn_masses[i2j(kk,iPAN )] = +ir[47] -ir[48]
	net_rxn_masses[i2j(kk,iCRO )] = -ir[68]
	net_rxn_masses[i2j(kk,iNTR )] = +ir[55] +0.1*ir[64] +ir[68] +0.8*ir[78]\
									+0.85*ir[94] +0.8*ir[96]
	net_rxn_masses[i2j(kk,iROR )] = -ir[55]
	net_rxn_masses[i2j(kk,iISOP)] = -ir[96]
	net_rxn_masses[i2j(kk,iISPD)] = -ir[94]
	
	
	
	############## PROCESSES #############################################
	
	# Initial Concentrations:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = +ip[iNO  ][iInitial]
	net_rxn_masses[i2j(kk,iNO2 )] = +ip[iNO2 ][iInitial]
	net_rxn_masses[i2j(kk,iO3  )] = +ip[iO3  ][iInitial]
	net_rxn_masses[i2j(kk,iOLE )] = +ip[iOLE ][iInitial]
	net_rxn_masses[i2j(kk,iPAN )] = +ip[iPAN ][iInitial]
	net_rxn_masses[i2j(kk,iN2O5)] = +ip[iN2O5][iInitial]
	net_rxn_masses[i2j(kk,iPAR )] = +ip[iPAR ][iInitial]
	net_rxn_masses[i2j(kk,iTOL )] = +ip[iTOL ][iInitial]
	net_rxn_masses[i2j(kk,iXYL )] = +ip[iXYL ][iInitial]
	net_rxn_masses[i2j(kk,iFORM)] = +ip[iFORM][iInitial]
	net_rxn_masses[i2j(kk,iALD2)] = +ip[iALD2][iInitial]
	net_rxn_masses[i2j(kk,iETH )] = +ip[iETH ][iInitial]
	net_rxn_masses[i2j(kk,iCRES)] = +ip[iCRES][iInitial]
	net_rxn_masses[i2j(kk,iMGLY)] = +ip[iMGLY][iInitial]
	net_rxn_masses[i2j(kk,iOPEN)] = +ip[iOPEN][iInitial]
	net_rxn_masses[i2j(kk,iPNA )] = +ip[iPNA ][iInitial]
	net_rxn_masses[i2j(kk,iCO  )] = +ip[iCO  ][iInitial]
	net_rxn_masses[i2j(kk,iHONO)] = +ip[iHONO][iInitial]
	net_rxn_masses[i2j(kk,iH2O2)] = +ip[iH2O2][iInitial]
	net_rxn_masses[i2j(kk,iHNO3)] = +ip[iHNO3][iInitial]
	net_rxn_masses[i2j(kk,iISOP)] = +ip[iISOP][iInitial]
	net_rxn_masses[i2j(kk,iMEOH)] = +ip[iMEOH][iInitial]
	net_rxn_masses[i2j(kk,iETOH)] = +ip[iETOH][iInitial]
	
	
	# Chemistry:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iChemistry]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iChemistry]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iChemistry]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iChemistry]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iChemistry]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iChemistry]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iChemistry]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iChemistry]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iChemistry]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iChemistry]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iChemistry]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iChemistry]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iChemistry]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iChemistry]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iChemistry]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iChemistry]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iChemistry]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iChemistry]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iChemistry]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iChemistry]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iChemistry]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iChemistry]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iChemistry]
	
	# Area Emissions:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iEmiss_Area]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iEmiss_Area]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iEmiss_Area]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iEmiss_Area]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iEmiss_Area]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iEmiss_Area]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iEmiss_Area]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iEmiss_Area]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iEmiss_Area]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iEmiss_Area]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iEmiss_Area]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iEmiss_Area]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iEmiss_Area]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iEmiss_Area]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iEmiss_Area]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iEmiss_Area]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iEmiss_Area]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iEmiss_Area]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iEmiss_Area]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iEmiss_Area]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iEmiss_Area]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iEmiss_Area]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iEmiss_Area]
	
	# Point Emissions:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iEmiss_Pnt]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iEmiss_Pnt]
	
	
	# TOTAL EMISSIONS:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = +ip[iNO  ][iEmiss_Area] +ip[iNO  ][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iNO2 )] = +ip[iNO2 ][iEmiss_Area] +ip[iNO2 ][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iO3  )] = +ip[iO3  ][iEmiss_Area] +ip[iO3  ][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iOLE )] = +ip[iOLE ][iEmiss_Area] +ip[iOLE ][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iPAN )] = +ip[iPAN ][iEmiss_Area] +ip[iPAN ][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iN2O5)] = +ip[iN2O5][iEmiss_Area] +ip[iN2O5][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iPAR )] = +ip[iPAR ][iEmiss_Area] +ip[iPAR ][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iTOL )] = +ip[iTOL ][iEmiss_Area] +ip[iTOL ][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iXYL )] = +ip[iXYL ][iEmiss_Area] +ip[iXYL ][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iFORM)] = +ip[iFORM][iEmiss_Area] +ip[iFORM][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iALD2)] = +ip[iALD2][iEmiss_Area] +ip[iALD2][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iETH )] = +ip[iETH ][iEmiss_Area] +ip[iETH ][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iCRES)] = +ip[iCRES][iEmiss_Area] +ip[iCRES][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iMGLY)] = +ip[iMGLY][iEmiss_Area] +ip[iMGLY][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iOPEN)] = +ip[iOPEN][iEmiss_Area] +ip[iOPEN][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iPNA )] = +ip[iPNA ][iEmiss_Area] +ip[iPNA ][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iCO  )] = +ip[iCO  ][iEmiss_Area] +ip[iCO  ][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iHONO)] = +ip[iHONO][iEmiss_Area] +ip[iHONO][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iH2O2)] = +ip[iH2O2][iEmiss_Area] +ip[iH2O2][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iHNO3)] = +ip[iHNO3][iEmiss_Area] +ip[iHNO3][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iISOP)] = +ip[iISOP][iEmiss_Area] +ip[iISOP][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iMEOH)] = +ip[iMEOH][iEmiss_Area] +ip[iMEOH][iEmiss_Pnt] 
	net_rxn_masses[i2j(kk,iETOH)] = +ip[iETOH][iEmiss_Area] +ip[iETOH][iEmiss_Pnt] 
	
	# PiG:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iEmiss_PiG]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iEmiss_PiG]
	
	
	# west b. advection:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iAdv_W]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iAdv_W]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iAdv_W]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iAdv_W]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iAdv_W]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iAdv_W]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iAdv_W]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iAdv_W]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iAdv_W]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iAdv_W]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iAdv_W]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iAdv_W]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iAdv_W]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iAdv_W]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iAdv_W]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iAdv_W]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iAdv_W]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iAdv_W]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iAdv_W]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iAdv_W]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iAdv_W]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iAdv_W]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iAdv_W]
	
	# east b. advection:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iAdv_E]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iAdv_E]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iAdv_E]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iAdv_E]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iAdv_E]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iAdv_E]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iAdv_E]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iAdv_E]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iAdv_E]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iAdv_E]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iAdv_E]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iAdv_E]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iAdv_E]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iAdv_E]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iAdv_E]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iAdv_E]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iAdv_E]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iAdv_E]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iAdv_E]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iAdv_E]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iAdv_E]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iAdv_E]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iAdv_E]
	
	# south b. advection:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iAdv_S]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iAdv_S]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iAdv_S]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iAdv_S]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iAdv_S]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iAdv_S]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iAdv_S]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iAdv_S]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iAdv_S]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iAdv_S]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iAdv_S]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iAdv_S]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iAdv_S]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iAdv_S]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iAdv_S]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iAdv_S]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iAdv_S]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iAdv_S]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iAdv_S]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iAdv_S]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iAdv_S]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iAdv_S]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iAdv_S]
	
	# north b. advection:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iAdv_N]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iAdv_N]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iAdv_N]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iAdv_N]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iAdv_N]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iAdv_N]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iAdv_N]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iAdv_N]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iAdv_N]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iAdv_N]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iAdv_N]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iAdv_N]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iAdv_N]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iAdv_N]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iAdv_N]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iAdv_N]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iAdv_N]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iAdv_N]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iAdv_N]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iAdv_N]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iAdv_N]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iAdv_N]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iAdv_N]
	
	# bottom b. advection:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iAdv_B]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iAdv_B]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iAdv_B]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iAdv_B]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iAdv_B]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iAdv_B]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iAdv_B]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iAdv_B]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iAdv_B]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iAdv_B]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iAdv_B]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iAdv_B]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iAdv_B]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iAdv_B]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iAdv_B]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iAdv_B]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iAdv_B]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iAdv_B]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iAdv_B]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iAdv_B]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iAdv_B]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iAdv_B]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iAdv_B]
	
	# top b. advection:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iAdv_T]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iAdv_T]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iAdv_T]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iAdv_T]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iAdv_T]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iAdv_T]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iAdv_T]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iAdv_T]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iAdv_T]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iAdv_T]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iAdv_T]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iAdv_T]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iAdv_T]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iAdv_T]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iAdv_T]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iAdv_T]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iAdv_T]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iAdv_T]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iAdv_T]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iAdv_T]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iAdv_T]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iAdv_T]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iAdv_T]
	
	# dil. in the vert.:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iDil_V]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iDil_V]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iDil_V]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iDil_V]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iDil_V]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iDil_V]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iDil_V]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iDil_V]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iDil_V]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iDil_V]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iDil_V]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iDil_V]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iDil_V]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iDil_V]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iDil_V]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iDil_V]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iDil_V]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iDil_V]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iDil_V]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iDil_V]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iDil_V]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iDil_V]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iDil_V]
	
	# west b. diffusion:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iDif_W]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iDif_W]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iDif_W]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iDif_W]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iDif_W]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iDif_W]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iDif_W]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iDif_W]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iDif_W]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iDif_W]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iDif_W]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iDif_W]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iDif_W]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iDif_W]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iDif_W]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iDif_W]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iDif_W]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iDif_W]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iDif_W]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iDif_W]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iDif_W]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iDif_W]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iDif_W]
	
	# east b. diffusion:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iDif_E]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iDif_E]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iDif_E]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iDif_E]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iDif_E]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iDif_E]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iDif_E]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iDif_E]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iDif_E]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iDif_E]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iDif_E]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iDif_E]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iDif_E]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iDif_E]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iDif_E]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iDif_E]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iDif_E]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iDif_E]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iDif_E]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iDif_E]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iDif_E]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iDif_E]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iDif_E]
	
	# south b. diffusion:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iDif_S]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iDif_S]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iDif_S]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iDif_S]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iDif_S]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iDif_S]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iDif_S]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iDif_S]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iDif_S]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iDif_S]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iDif_S]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iDif_S]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iDif_S]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iDif_S]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iDif_S]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iDif_S]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iDif_S]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iDif_S]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iDif_S]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iDif_S]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iDif_S]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iDif_S]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iDif_S]
	
	# north b. diffusion:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iDif_N]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iDif_N]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iDif_N]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iDif_N]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iDif_N]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iDif_N]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iDif_N]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iDif_N]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iDif_N]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iDif_N]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iDif_N]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iDif_N]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iDif_N]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iDif_N]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iDif_N]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iDif_N]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iDif_N]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iDif_N]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iDif_N]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iDif_N]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iDif_N]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iDif_N]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iDif_N]
	
	# bottom b. diffusion:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iDif_B]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iDif_B]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iDif_B]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iDif_B]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iDif_B]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iDif_B]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iDif_B]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iDif_B]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iDif_B]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iDif_B]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iDif_B]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iDif_B]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iDif_B]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iDif_B]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iDif_B]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iDif_B]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iDif_B]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iDif_B]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iDif_B]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iDif_B]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iDif_B]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iDif_B]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iDif_B]
	
	# top b. diffusion:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iDif_T]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iDif_T]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iDif_T]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iDif_T]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iDif_T]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iDif_T]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iDif_T]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iDif_T]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iDif_T]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iDif_T]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iDif_T]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iDif_T]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iDif_T]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iDif_T]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iDif_T]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iDif_T]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iDif_T]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iDif_T]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iDif_T]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iDif_T]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iDif_T]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iDif_T]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iDif_T]
	
	# dry deposition:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iDep_D]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iDep_D]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iDep_D]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iDep_D]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iDep_D]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iDep_D]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iDep_D]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iDep_D]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iDep_D]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iDep_D]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iDep_D]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iDep_D]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iDep_D]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iDep_D]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iDep_D]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iDep_D]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iDep_D]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iDep_D]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iDep_D]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iDep_D]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iDep_D]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iDep_D]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iDep_D]
	
	# wet deposition:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iDep_W]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iDep_W]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iDep_W]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iDep_W]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iDep_W]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iDep_W]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iDep_W]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iDep_W]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iDep_W]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iDep_W]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iDep_W]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iDep_W]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iDep_W]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iDep_W]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iDep_W]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iDep_W]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iDep_W]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iDep_W]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iDep_W]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iDep_W]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iDep_W]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iDep_W]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iDep_W]
	
	# aerosol chemistry:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iChem_Aero]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iChem_Aero]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iChem_Aero]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iChem_Aero]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iChem_Aero]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iChem_Aero]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iChem_Aero]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iChem_Aero]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iChem_Aero]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iChem_Aero]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iChem_Aero]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iChem_Aero]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iChem_Aero]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iChem_Aero]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iChem_Aero]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iChem_Aero]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iChem_Aero]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iChem_Aero]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iChem_Aero]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iChem_Aero]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iChem_Aero]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iChem_Aero]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iChem_Aero]
	
	# dilution:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iDilut]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iDilut]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iDilut]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iDilut]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iDilut]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iDilut]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iDilut]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iDilut]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iDilut]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iDilut]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iDilut]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iDilut]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iDilut]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iDilut]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iDilut]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iDilut]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iDilut]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iDilut]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iDilut]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iDilut]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iDilut]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iDilut]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iDilut]
	
	# en(de)trainment:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iTrain]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iTrain]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iTrain]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iTrain]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iTrain]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iTrain]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iTrain]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iTrain]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iTrain]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iTrain]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iTrain]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iTrain]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iTrain]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iTrain]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iTrain]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iTrain]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iTrain]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iTrain]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iTrain]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iTrain]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iTrain]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iTrain]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iTrain]
	
	# HORIZONTAL TRANSPORTATION:
	
	kk += 1
	i
	net_rxn_masses[i2j(kk,iNO  )] = +ip[iNO  ][iAdv_W] +ip[iNO  ][iAdv_E] +ip[iNO  ][iAdv_S] +ip[iNO  ][iAdv_N] +ip[iNO  ][iDif_W] +ip[iNO  ][iDif_E] +ip[iNO  ][iDif_S] +ip[iNO  ][iDif_N]
	net_rxn_masses[i2j(kk,iNO2 )] = +ip[iNO2 ][iAdv_W] +ip[iNO2 ][iAdv_E] +ip[iNO2 ][iAdv_S] +ip[iNO2 ][iAdv_N] +ip[iNO2 ][iDif_W] +ip[iNO2 ][iDif_E] +ip[iNO2 ][iDif_S] +ip[iNO2 ][iDif_N]
	net_rxn_masses[i2j(kk,iO3  )] = +ip[iO3  ][iAdv_W] +ip[iO3  ][iAdv_E] +ip[iO3  ][iAdv_S] +ip[iO3  ][iAdv_N] +ip[iO3  ][iDif_W] +ip[iO3  ][iDif_E] +ip[iO3  ][iDif_S] +ip[iO3  ][iDif_N]
	net_rxn_masses[i2j(kk,iOLE )] = +ip[iOLE ][iAdv_W] +ip[iOLE ][iAdv_E] +ip[iOLE ][iAdv_S] +ip[iOLE ][iAdv_N] +ip[iOLE ][iDif_W] +ip[iOLE ][iDif_E] +ip[iOLE ][iDif_S] +ip[iOLE ][iDif_N]
	net_rxn_masses[i2j(kk,iPAN )] = +ip[iPAN ][iAdv_W] +ip[iPAN ][iAdv_E] +ip[iPAN ][iAdv_S] +ip[iPAN ][iAdv_N] +ip[iPAN ][iDif_W] +ip[iPAN ][iDif_E] +ip[iPAN ][iDif_S] +ip[iPAN ][iDif_N]
	net_rxn_masses[i2j(kk,iN2O5)] = +ip[iN2O5][iAdv_W] +ip[iN2O5][iAdv_E] +ip[iN2O5][iAdv_S] +ip[iN2O5][iAdv_N] +ip[iN2O5][iDif_W] +ip[iN2O5][iDif_E] +ip[iN2O5][iDif_S] +ip[iN2O5][iDif_N]
	net_rxn_masses[i2j(kk,iPAR )] = +ip[iPAR ][iAdv_W] +ip[iPAR ][iAdv_E] +ip[iPAR ][iAdv_S] +ip[iPAR ][iAdv_N] +ip[iPAR ][iDif_W] +ip[iPAR ][iDif_E] +ip[iPAR ][iDif_S] +ip[iPAR ][iDif_N]
	net_rxn_masses[i2j(kk,iTOL )] = +ip[iTOL ][iAdv_W] +ip[iTOL ][iAdv_E] +ip[iTOL ][iAdv_S] +ip[iTOL ][iAdv_N] +ip[iTOL ][iDif_W] +ip[iTOL ][iDif_E] +ip[iTOL ][iDif_S] +ip[iTOL ][iDif_N]
	net_rxn_masses[i2j(kk,iXYL )] = +ip[iXYL ][iAdv_W] +ip[iXYL ][iAdv_E] +ip[iXYL ][iAdv_S] +ip[iXYL ][iAdv_N] +ip[iXYL ][iDif_W] +ip[iXYL ][iDif_E] +ip[iXYL ][iDif_S] +ip[iXYL ][iDif_N]
	net_rxn_masses[i2j(kk,iFORM)] = +ip[iFORM][iAdv_W] +ip[iFORM][iAdv_E] +ip[iFORM][iAdv_S] +ip[iFORM][iAdv_N] +ip[iFORM][iDif_W] +ip[iFORM][iDif_E] +ip[iFORM][iDif_S] +ip[iFORM][iDif_N]
	net_rxn_masses[i2j(kk,iALD2)] = +ip[iALD2][iAdv_W] +ip[iALD2][iAdv_E] +ip[iALD2][iAdv_S] +ip[iALD2][iAdv_N] +ip[iALD2][iDif_W] +ip[iALD2][iDif_E] +ip[iALD2][iDif_S] +ip[iALD2][iDif_N]
	net_rxn_masses[i2j(kk,iETH )] = +ip[iETH ][iAdv_W] +ip[iETH ][iAdv_E] +ip[iETH ][iAdv_S] +ip[iETH ][iAdv_N] +ip[iETH ][iDif_W] +ip[iETH ][iDif_E] +ip[iETH ][iDif_S] +ip[iETH ][iDif_N]
	net_rxn_masses[i2j(kk,iCRES)] = +ip[iCRES][iAdv_W] +ip[iCRES][iAdv_E] +ip[iCRES][iAdv_S] +ip[iCRES][iAdv_N] +ip[iCRES][iDif_W] +ip[iCRES][iDif_E] +ip[iCRES][iDif_S] +ip[iCRES][iDif_N]
	net_rxn_masses[i2j(kk,iMGLY)] = +ip[iMGLY][iAdv_W] +ip[iMGLY][iAdv_E] +ip[iMGLY][iAdv_S] +ip[iMGLY][iAdv_N] +ip[iMGLY][iDif_W] +ip[iMGLY][iDif_E] +ip[iMGLY][iDif_S] +ip[iMGLY][iDif_N]
	net_rxn_masses[i2j(kk,iOPEN)] = +ip[iOPEN][iAdv_W] +ip[iOPEN][iAdv_E] +ip[iOPEN][iAdv_S] +ip[iOPEN][iAdv_N] +ip[iOPEN][iDif_W] +ip[iOPEN][iDif_E] +ip[iOPEN][iDif_S] +ip[iOPEN][iDif_N]
	net_rxn_masses[i2j(kk,iPNA )] = +ip[iPNA ][iAdv_W] +ip[iPNA ][iAdv_E] +ip[iPNA ][iAdv_S] +ip[iPNA ][iAdv_N] +ip[iPNA ][iDif_W] +ip[iPNA ][iDif_E] +ip[iPNA ][iDif_S] +ip[iPNA ][iDif_N]
	net_rxn_masses[i2j(kk,iCO  )] = +ip[iCO  ][iAdv_W] +ip[iCO  ][iAdv_E] +ip[iCO  ][iAdv_S] +ip[iCO  ][iAdv_N] +ip[iCO  ][iDif_W] +ip[iCO  ][iDif_E] +ip[iCO  ][iDif_S] +ip[iCO  ][iDif_N]
	net_rxn_masses[i2j(kk,iHONO)] = +ip[iHONO][iAdv_W] +ip[iHONO][iAdv_E] +ip[iHONO][iAdv_S] +ip[iHONO][iAdv_N] +ip[iHONO][iDif_W] +ip[iHONO][iDif_E] +ip[iHONO][iDif_S] +ip[iHONO][iDif_N]
	net_rxn_masses[i2j(kk,iH2O2)] = +ip[iH2O2][iAdv_W] +ip[iH2O2][iAdv_E] +ip[iH2O2][iAdv_S] +ip[iH2O2][iAdv_N] +ip[iH2O2][iDif_W] +ip[iH2O2][iDif_E] +ip[iH2O2][iDif_S] +ip[iH2O2][iDif_N]
	net_rxn_masses[i2j(kk,iHNO3)] = +ip[iHNO3][iAdv_W] +ip[iHNO3][iAdv_E] +ip[iHNO3][iAdv_S] +ip[iHNO3][iAdv_N] +ip[iHNO3][iDif_W] +ip[iHNO3][iDif_E] +ip[iHNO3][iDif_S] +ip[iHNO3][iDif_N]
	net_rxn_masses[i2j(kk,iISOP)] = +ip[iISOP][iAdv_W] +ip[iISOP][iAdv_E] +ip[iISOP][iAdv_S] +ip[iISOP][iAdv_N] +ip[iISOP][iDif_W] +ip[iISOP][iDif_E] +ip[iISOP][iDif_S] +ip[iISOP][iDif_N]
	net_rxn_masses[i2j(kk,iMEOH)] = +ip[iMEOH][iAdv_W] +ip[iMEOH][iAdv_E] +ip[iMEOH][iAdv_S] +ip[iMEOH][iAdv_N] +ip[iMEOH][iDif_W] +ip[iMEOH][iDif_E] +ip[iMEOH][iDif_S] +ip[iMEOH][iDif_N]
	net_rxn_masses[i2j(kk,iETOH)] = +ip[iETOH][iAdv_W] +ip[iETOH][iAdv_E] +ip[iETOH][iAdv_S] +ip[iETOH][iAdv_N] +ip[iETOH][iDif_W] +ip[iETOH][iDif_E] +ip[iETOH][iDif_S] +ip[iETOH][iDif_N]


# VERTICAL TRANSPORTATION:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = +ip[iNO  ][iAdv_B] +ip[iNO  ][iAdv_T] +ip[iNO  ][iDif_B] +ip[iNO  ][iDif_T] +ip[iNO  ][iDil_V] +ip[iNO  ][iDilut] +ip[iNO  ][iTrain]
	net_rxn_masses[i2j(kk,iNO2 )] = +ip[iNO2 ][iAdv_B] +ip[iNO2 ][iAdv_T] +ip[iNO2 ][iDif_B] +ip[iNO2 ][iDif_T] +ip[iNO2 ][iDil_V] +ip[iNO2 ][iDilut] +ip[iNO2 ][iTrain]
	net_rxn_masses[i2j(kk,iO3  )] = +ip[iO3  ][iAdv_B] +ip[iO3  ][iAdv_T] +ip[iO3  ][iDif_B] +ip[iO3  ][iDif_T] +ip[iO3  ][iDil_V] +ip[iO3  ][iDilut] +ip[iO3  ][iTrain]
	net_rxn_masses[i2j(kk,iOLE )] = +ip[iOLE ][iAdv_B] +ip[iOLE ][iAdv_T] +ip[iOLE ][iDif_B] +ip[iOLE ][iDif_T] +ip[iOLE ][iDil_V] +ip[iOLE ][iDilut] +ip[iOLE ][iTrain]
	net_rxn_masses[i2j(kk,iPAN )] = +ip[iPAN ][iAdv_B] +ip[iPAN ][iAdv_T] +ip[iPAN ][iDif_B] +ip[iPAN ][iDif_T] +ip[iPAN ][iDil_V] +ip[iPAN ][iDilut] +ip[iPAN ][iTrain]
	net_rxn_masses[i2j(kk,iN2O5)] = +ip[iN2O5][iAdv_B] +ip[iN2O5][iAdv_T] +ip[iN2O5][iDif_B] +ip[iN2O5][iDif_T] +ip[iN2O5][iDil_V] +ip[iN2O5][iDilut] +ip[iN2O5][iTrain]
	net_rxn_masses[i2j(kk,iPAR )] = +ip[iPAR ][iAdv_B] +ip[iPAR ][iAdv_T] +ip[iPAR ][iDif_B] +ip[iPAR ][iDif_T] +ip[iPAR ][iDil_V] +ip[iPAR ][iDilut] +ip[iPAR ][iTrain]
	net_rxn_masses[i2j(kk,iTOL )] = +ip[iTOL ][iAdv_B] +ip[iTOL ][iAdv_T] +ip[iTOL ][iDif_B] +ip[iTOL ][iDif_T] +ip[iTOL ][iDil_V] +ip[iTOL ][iDilut] +ip[iTOL ][iTrain]
	net_rxn_masses[i2j(kk,iXYL )] = +ip[iXYL ][iAdv_B] +ip[iXYL ][iAdv_T] +ip[iXYL ][iDif_B] +ip[iXYL ][iDif_T] +ip[iXYL ][iDil_V] +ip[iXYL ][iDilut] +ip[iXYL ][iTrain]
	net_rxn_masses[i2j(kk,iFORM)] = +ip[iFORM][iAdv_B] +ip[iFORM][iAdv_T] +ip[iFORM][iDif_B] +ip[iFORM][iDif_T] +ip[iFORM][iDil_V] +ip[iFORM][iDilut] +ip[iFORM][iTrain]
	net_rxn_masses[i2j(kk,iALD2)] = +ip[iALD2][iAdv_B] +ip[iALD2][iAdv_T] +ip[iALD2][iDif_B] +ip[iALD2][iDif_T] +ip[iALD2][iDil_V] +ip[iALD2][iDilut] +ip[iALD2][iTrain]
	net_rxn_masses[i2j(kk,iETH )] = +ip[iETH ][iAdv_B] +ip[iETH ][iAdv_T] +ip[iETH ][iDif_B] +ip[iETH ][iDif_T] +ip[iETH ][iDil_V] +ip[iETH ][iDilut] +ip[iETH ][iTrain]
	net_rxn_masses[i2j(kk,iCRES)] = +ip[iCRES][iAdv_B] +ip[iCRES][iAdv_T] +ip[iCRES][iDif_B] +ip[iCRES][iDif_T] +ip[iCRES][iDil_V] +ip[iCRES][iDilut] +ip[iCRES][iTrain]
	net_rxn_masses[i2j(kk,iMGLY)] = +ip[iMGLY][iAdv_B] +ip[iMGLY][iAdv_T] +ip[iMGLY][iDif_B] +ip[iMGLY][iDif_T] +ip[iMGLY][iDil_V] +ip[iMGLY][iDilut] +ip[iMGLY][iTrain]
	net_rxn_masses[i2j(kk,iOPEN)] = +ip[iOPEN][iAdv_B] +ip[iOPEN][iAdv_T] +ip[iOPEN][iDif_B] +ip[iOPEN][iDif_T] +ip[iOPEN][iDil_V] +ip[iOPEN][iDilut] +ip[iOPEN][iTrain]
	net_rxn_masses[i2j(kk,iPNA )] = +ip[iPNA ][iAdv_B] +ip[iPNA ][iAdv_T] +ip[iPNA ][iDif_B] +ip[iPNA ][iDif_T] +ip[iPNA ][iDil_V] +ip[iPNA ][iDilut] +ip[iPNA ][iTrain]
	net_rxn_masses[i2j(kk,iCO  )] = +ip[iCO  ][iAdv_B] +ip[iCO  ][iAdv_T] +ip[iCO  ][iDif_B] +ip[iCO  ][iDif_T] +ip[iCO  ][iDil_V] +ip[iCO  ][iDilut] +ip[iCO  ][iTrain]
	net_rxn_masses[i2j(kk,iHONO)] = +ip[iHONO][iAdv_B] +ip[iHONO][iAdv_T] +ip[iHONO][iDif_B] +ip[iHONO][iDif_T] +ip[iHONO][iDil_V] +ip[iHONO][iDilut] +ip[iHONO][iTrain]
	net_rxn_masses[i2j(kk,iH2O2)] = +ip[iH2O2][iAdv_B] +ip[iH2O2][iAdv_T] +ip[iH2O2][iDif_B] +ip[iH2O2][iDif_T] +ip[iH2O2][iDil_V] +ip[iH2O2][iDilut] +ip[iH2O2][iTrain]
	net_rxn_masses[i2j(kk,iHNO3)] = +ip[iHNO3][iAdv_B] +ip[iHNO3][iAdv_T] +ip[iHNO3][iDif_B] +ip[iHNO3][iDif_T] +ip[iHNO3][iDil_V] +ip[iHNO3][iDilut] +ip[iHNO3][iTrain]
	net_rxn_masses[i2j(kk,iISOP)] = +ip[iISOP][iAdv_B] +ip[iISOP][iAdv_T] +ip[iISOP][iDif_B] +ip[iISOP][iDif_T] +ip[iISOP][iDil_V] +ip[iISOP][iDilut] +ip[iISOP][iTrain]
	net_rxn_masses[i2j(kk,iMEOH)] = +ip[iMEOH][iAdv_B] +ip[iMEOH][iAdv_T] +ip[iMEOH][iDif_B] +ip[iMEOH][iDif_T] +ip[iMEOH][iDil_V] +ip[iMEOH][iDilut] +ip[iMEOH][iTrain]
	net_rxn_masses[i2j(kk,iETOH)] = +ip[iETOH][iAdv_B] +ip[iETOH][iAdv_T] +ip[iETOH][iDif_B] +ip[iETOH][iDif_T] +ip[iETOH][iDil_V] +ip[iETOH][iDilut] +ip[iETOH][iTrain]
	
	
	# final concentration:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iFinal]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iFinal]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iFinal]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iFinal]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iFinal]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iFinal]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iFinal]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iFinal]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iFinal]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iFinal]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iFinal]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iFinal]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iFinal]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iFinal]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iFinal]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iFinal]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iFinal]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iFinal]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iFinal]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iFinal]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iFinal]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iFinal]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iFinal]
	
	# units conversion:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iUn_Con]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iUn_Con]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iUn_Con]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iUn_Con]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iUn_Con]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iUn_Con]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iUn_Con]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iUn_Con]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iUn_Con]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iUn_Con]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iUn_Con]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iUn_Con]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iUn_Con]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iUn_Con]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iUn_Con]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iUn_Con]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iUn_Con]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iUn_Con]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iUn_Con]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iUn_Con]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iUn_Con]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iUn_Con]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iUn_Con]
	
	# average cell volume:
	
	kk += 1
	
	net_rxn_masses[i2j(kk,iNO  )] = ip[iNO  ][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iNO2 )] = ip[iNO2 ][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iO3  )] = ip[iO3  ][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iOLE )] = ip[iOLE ][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iPAN )] = ip[iPAN ][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iN2O5)] = ip[iN2O5][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iPAR )] = ip[iPAR ][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iTOL )] = ip[iTOL ][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iXYL )] = ip[iXYL ][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iFORM)] = ip[iFORM][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iALD2)] = ip[iALD2][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iETH )] = ip[iETH ][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iCRES)] = ip[iCRES][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iMGLY)] = ip[iMGLY][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iOPEN)] = ip[iOPEN][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iPNA )] = ip[iPNA ][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iCO  )] = ip[iCO  ][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iHONO)] = ip[iHONO][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iH2O2)] = ip[iH2O2][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iHNO3)] = ip[iHNO3][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iISOP)] = ip[iISOP][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iMEOH)] = ip[iMEOH][iAv_Cel_Vol]
	net_rxn_masses[i2j(kk,iETOH)] = ip[iETOH][iAv_Cel_Vol]
	
	
	kk += 1
	
	
	
	
	
	# <<<add more net reactions assignments here....>>>
	
	
	# save these hourly results for later output....
	hour_number.append(time)
	hourly_net_rxn_masses.append(copy.deepcopy(net_rxn_masses))
	
	# accumulate this hour's net masses into a daily total..
	for i in range(0,len(net_rxn_masses)):
		daily_net_rxn_masses[i] += net_rxn_masses[i]
	
	# all reactions computed, read in a new set of ir and go back to top
	(time, ir, ip) = get_irr_data(fin);
	
	### end of the hour loop 


######### repeat to here for each set of net reactions ######

## >>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
##  P A R T   T H R E E  : C A L C U L A T E  T H E  S U M M A R I E S
## ===================================================================

# Set up data representation for the diagram_values

# Start with info about the diagram sections...
#   allocate vector of names of the net_rxns
diagram_sect_names = []
diagram_sect_start = []
#   both indexed by next_diag_section

# provide a counter to give integer index to sections
sec_id_num = 0

# the number for the next diagram section to be added...
next_diag_section = 0
num_diag_sections = 0


# allocate diagram values table storage

# allocate section labels vector (label of a row)
section_labels = []
#   adjacent elements in this string vector hold the row label

num_hrs = len(hour_number)

# allocate a list of lists of hourly values
hourly_diagram_values = []
#   adjacent elements in this vector hold hourly_values vectors for num_hrs for
#    whatever the property is in a given row.

daily_diagram_values = []
#   accumulate the hourly diagram values over all time.

# next location in net_rxn_masses vector to append another net rxn set
#   and in the net_rxn_spcname vector to insert the iSPC number
jstart_next_diagram_row = 0

# allocate a row accumulator vector for general use...
row_acc = [0.0]*num_hrs
#      
total_acc = 0
#      for the sum of the row_accumulator

def copyhours(hourlynet = [], elem = 0, row = []):
	"copy same element at each hour to an hourly vector of values"
	
	if len(row) == num_hrs :
		for t in range(0,num_hrs):
			row[t] = hourlynet[t][elem]
	else :
		raise IndexError
	
def accumulate_rxn(hourlynet = [], elem = 0, acc = []):
	"accumulate same element at each hour to an hourly vector of values"
	if len(cc) == num_hrs :
		for t in range(0,num_hrs):
			acc[t]      += hourlynet[t][elem]
	else :
		raise IndexError
	
def accumulate_row(row = [], acc = []):
	"accumulate same element at each hour to an hourly vector of values"
	if len(acc) == num_hrs :
		for t in range(0,num_hrs):
			acc[t]      += row[t]
	else :
		raise IndexError
	


# @@@@@@@@@ D I A G R A M   S E C T I O N S @@@@@@@@@@@@
#
#   1) "C2O3     new, prod, loss" 
#   2) "XO2/XO2N new, prod, loss"
#   3) "HO2      new, prod, loss"
#   4) "OH       new, prod, loss"



## begin Diagram Section storage allocation
a_row   = [0.0]*num_hrs

#-------------------------------------------------------
### The Diagram Section for ** C2O3 prod and loss **
this_diag_section = next_diag_section  
next_diag_section += 1

#    ... save the name and starting location of this net_rxn
diagram_sect_names.append("C2O3 New, Prod, and Loss")
n_C2O3_prodloss = sec_id_num  # provide an interger with a name for this section
sec_id_num += 1

this_jstart = jstart_next_diagram_row
diagram_sect_start.append(this_jstart)

jj = this_jstart

section_labels.append(" Ald+hv  C2O3") 
copyhours(hourly_net_rxn_masses, i2j(n_Aldhvrad, iC2O3 ), a_row)
hourly_diagram_values.append([e for e in a_row])
hourly_tot_c2o3p = [e for e in a_row]

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_Aldhvrad, iC2O3 )])

jj += 1

section_labels.append(" Ox+org  C2O3") 
copyhours(hourly_net_rxn_masses, i2j(n_OxOrgrad, iC2O3 ), a_row)
hourly_diagram_values.append([e for e in a_row])
accumulate_row(a_row,hourly_tot_c2o3p)

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_OxOrgrad, iC2O3 )])

jj += 1

section_labels.append(" NO3+org C2O3") 
copyhours(hourly_net_rxn_masses, i2j(n_NO3Orgrad, iC2O3 ), a_row)
hourly_diagram_values.append([e for e in a_row])
accumulate_row(a_row,hourly_tot_c2o3p)

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_NO3Orgrad, iC2O3 )])

jj += 1

section_labels.append("Total  C2O3 New") 
#save hourly and daily new C3O3 for later use with OH
hourly_newc2o3 = [e for e in hourly_tot_c2o3p]
daily_newc2o3  = sum(hourly_newc2o3)
hourly_diagram_values.append([e for e in hourly_newc2o3])

daily_diagram_values.append(daily_newc2o3)

jj += 1

section_labels.append("OH+Org C2O3 Prod") 
copyhours(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iC2O3 ), a_row)
hourly_diagram_values.append([e for e in a_row])
accumulate_row(a_row,hourly_tot_c2o3p)

daily_diagram_values.append(sum(hourly_tot_c2o3p))

jj += 1

section_labels.append("Total  C2O3 Prod") 
hourly_diagram_values.append([e for e in hourly_tot_c2o3p])
# remember the total C2O3 produced by sum of hourly_tot_c2o3p
daily_tot_c2o3p = sum(hourly_tot_c2o3p)
daily_diagram_values.append(daily_tot_c2o3p)


jj += 1

section_labels.append("C2O3+NO2 PAN Prod") 
copyhours(hourly_net_rxn_masses, i2j(n_PANProd, iPAN ), a_row)
hourly_diagram_values.append([e for e in a_row])
tot_c2o3l = [e for e in a_row]  # PAN plus C2O2 + NO2

daily_diagram_values.append(sum(tot_c2o3l))

jj += 1

section_labels.append("C2O3+NO Loss")
copyhours(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iC2O3 ), a_row)
# these are negative values, chg sign
b_row = [-e for e in a_row]
hourly_diagram_values.append([e for e in b_row])
accumulate_row(b_row,tot_c2o3l)

daily_diagram_values.append(sum(b_row))

jj += 1

section_labels.append("  NO2 prod") 
copyhours(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iNO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
daily_no2_prod = sum(a_row)
daily_diagram_values.append(daily_no2_prod)

# remember the NO2 produced by C2O3
hourly_no2_prod_c2o3 = [ n/c for (n,c) in zip(a_row,hourly_tot_c2o3p)]
daily_no2_prod_c2o3 = daily_no2_prod/daily_tot_c2o3p

jj += 1

section_labels.append("  XO2 prod") 
copyhours(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iXO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

# remember the XO2 produced by C2O3
hourly_xo2_prod_c2o3 = [ x/c for (x,c) in zip(a_row,hourly_tot_c2o3p)]
daily_xo2_prod_c2o3 = sum(a_row)/daily_tot_c2o3p

jj += 1

section_labels.append("  OOX prod") 
copyhours(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iOOX ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("Total C2O3 Loss") 
hourly_diagram_values.append([e for e in tot_c2o3l])
daily_diagram_values.append(sum(tot_c2o3l))

jj += 1

section_labels.append("NO2 Prod / C2O3 Prod") 
hourly_diagram_values.append([e for e in hourly_no2_prod_c2o3])
daily_diagram_values.append(daily_no2_prod_c2o3)

jj += 1

section_labels.append("XO2 Prod / C2O3 Prod") 
hourly_diagram_values.append([e for e in hourly_xo2_prod_c2o3])
daily_diagram_values.append(daily_xo2_prod_c2o3)

jj += 1

# how many elements in vector are needed.
num_rows_in_section = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_diagram_row  = jj



#-------------------------------------------------------
### The Diagram Section for ** XO2 Prod and Loss **
this_diag_section = next_diag_section  
next_diag_section += 1

#    ... save the name and starting location of this net_rxn
diagram_sect_names.append("XO2 New, Prod, and Loss")
n_XO2_prodloss = sec_id_num  # provide an interger with a name for this section
sec_id_num += 1

this_jstart = jstart_next_diagram_row
diagram_sect_start.append(this_jstart)

jj = this_jstart

section_labels.append(" Ald+hv  XO2") 
copyhours(hourly_net_rxn_masses, i2j(n_Aldhvrad, iXO2 ), a_row)
# start acculating the sources of XO2
hourly_tot_xo2p = [e for e in a_row]
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_Aldhvrad, iXO2 )])

jj += 1

section_labels.append(" Ox+org  XO2") 
copyhours(hourly_net_rxn_masses, i2j(n_OxOrgrad, iXO2 ), a_row)
accumulate_row(a_row,hourly_tot_xo2p)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_OxOrgrad, iXO2 )])

jj += 1

section_labels.append(" NO3+org XO2") 
copyhours(hourly_net_rxn_masses, i2j(n_NO3Orgrad, iXO2 ), a_row)
accumulate_row(a_row,hourly_tot_xo2p)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_NO3Orgrad, iXO2 )])

jj += 1

section_labels.append("Total  XO2  New") 
# save the new XO2
hourly_newxo2 = [e for e in hourly_tot_xo2p]
daily_newxo2  = sum(hourly_newxo2)
hourly_diagram_values.append([e for e in hourly_newxo2])

daily_diagram_values.append(daily_newxo2)

jj += 1

section_labels.append("C2O3   XO2  Prod") 
copyhours(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iXO2 ), a_row)
accumulate_row(a_row,hourly_tot_xo2p)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("OH+Org XO2  Prod") 
copyhours(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iXO2 ), a_row)
accumulate_row(a_row,hourly_tot_xo2p)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("Total  XO2  Prod") 
# remember the total XO2 produced by saving hourly_tot_xo2p
daily_tot_xo2p = sum(hourly_tot_xo2p)
hourly_diagram_values.append([e for e in hourly_tot_xo2p])

daily_diagram_values.append(daily_tot_xo2p)

jj += 1

section_labels.append(" Ox+org  XO2N") 
copyhours(hourly_net_rxn_masses, i2j(n_OxOrgrad, iXO2N ), a_row)
# accumulate the XO2N 
hourly_tot_xo2n = [e for e in a_row]
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_OxOrgrad, iXO2N )])

jj += 1

section_labels.append(" NO3+org XO2N") 
copyhours(hourly_net_rxn_masses, i2j(n_NO3Orgrad, iXO2N ), a_row)
accumulate_row(a_row,hourly_tot_xo2n)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_NO3Orgrad, iXO2N )])

jj += 1

section_labels.append("Total  XO2N New") 
hourly_diagram_values.append([e for e in hourly_tot_xo2n])

daily_diagram_values.append(sum(hourly_tot_xo2n))

jj += 1

section_labels.append("OH+Org XO2N Prod") 
copyhours(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iXO2N ), a_row)
accumulate_row(a_row,hourly_tot_xo2n)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("Total  XO2N Prod") 
hourly_diagram_values.append([e for e in hourly_tot_xo2n])
daily_diagram_values.append(sum(hourly_tot_xo2n))

jj += 1

section_labels.append("Total  XO2* Prod")
accumulate_row(hourly_tot_xo2p,hourly_tot_xo2n)
hourly_diagram_values.append([e for e in hourly_tot_xo2n])
daily_diagram_values.append(sum(hourly_tot_xo2n))

jj += 1

section_labels.append("XO2 +NO Loss")
copyhours(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iXO2 ), a_row)
# these are negative values, chg sign
b_row = [-e for e in a_row]
hourly_diagram_values.append([e for e in b_row])
tot_rows = [e for e in b_row]
daily_diagram_values.append(sum(b_row))

jj += 1

section_labels.append("XO2N+NO Loss")
copyhours(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iXO2N ), a_row)
# these are negative values, chg sign
b_row = [-e for e in a_row]
hourly_diagram_values.append([e for e in b_row])
accumulate_row(b_row,tot_rows)
daily_diagram_values.append(sum(b_row))

jj += 1

section_labels.append("Total XO2* Loss") 
hourly_diagram_values.append([e for e in tot_rows])
daily_diagram_values.append(sum(tot_rows))

jj += 1

section_labels.append("  NO  loss") 
copyhours(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iNO ), a_row)
# these are negative values, chg sign
b_row = [-e for e in a_row]
hourly_diagram_values.append([e for e in b_row])

daily_diagram_values.append(sum(b_row))

jj += 1

section_labels.append("  NO2 prod") 
copyhours(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iNO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
daily_no2_prod = sum(a_row)

daily_diagram_values.append(daily_no2_prod)

# remember NO2 produced per XO2 produced
no2_prod_xo2 = [ n/x for (n,x) in zip(a_row,hourly_tot_xo2p)]
daily_no2_prod_xo2 = daily_no2_prod/daily_tot_xo2p
jj += 1

section_labels.append("  HO2 prod") 
copyhours(hourly_net_rxn_masses, i2j(n_XO2NOOxid, ixHO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
daily_ho2_prod = sum(a_row)

daily_diagram_values.append(daily_ho2_prod)

# remember HO2 produced per XO2 produced
hourly_ho2_prod_xo2 = [ h/x for (h,x) in zip(a_row,hourly_tot_xo2p)]
daily_ho2_prod_xo2 = daily_ho2_prod/daily_tot_xo2p

jj += 1

section_labels.append("  NTR prod") 
copyhours(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iNTR ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append(" *OOX prod") 
copyhours(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iOOX ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("NO2 Prod / XO2 Prod") 
hourly_diagram_values.append([e for e in no2_prod_xo2])
daily_diagram_values.append(daily_no2_prod_xo2)

jj += 1

section_labels.append("HO2 Prod / XO2 Prod") 
hourly_diagram_values.append([e for e in hourly_ho2_prod_xo2])
daily_diagram_values.append(daily_ho2_prod_xo2)

jj += 1

# how many elements in vector are needed.
num_rows_in_section = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_diagram_row  = jj


#-------------------------------------------------------
### The Diagram Section for ** HO2 prod and loss **
this_diag_section = next_diag_section  
next_diag_section += 1

#    ... save the name and starting location of this net_rxn
diagram_sect_names.append("HO2 New, Prod, and Loss")
sec_id_num += 1


this_jstart = jstart_next_diagram_row
diagram_sect_start.append(this_jstart)

jj = this_jstart

section_labels.append(" Ald+hv  HO2") 
copyhours(hourly_net_rxn_masses, i2j(n_Aldhvrad, iHO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
# initialize total ho2 accumulator...
hourly_tot_ho2p = [e for e in a_row]

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_Aldhvrad, iHO2 )])

jj += 1

section_labels.append(" Ox+org  HO2") 
copyhours(hourly_net_rxn_masses, i2j(n_OxOrgrad, iHO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
# add in this HO2..
accumulate_row(a_row,hourly_tot_ho2p)

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_OxOrgrad, iHO2 )])

jj += 1

section_labels.append("Total  HO2 New") 
# save total new HO2
hourly_newho2 = [e for e in hourly_tot_ho2p]
daily_newho2  = sum(hourly_newho2)
hourly_diagram_values.append([e for e in hourly_newho2])

daily_diagram_values.append(daily_newho2)

jj += 1

section_labels.append("OH+Org HO2 Prod") 
copyhours(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iHO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
# add in this HO2..
accumulate_row(a_row,hourly_tot_ho2p)

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("XO2+NO HO2 Prod") 
copyhours(hourly_net_rxn_masses, i2j(n_XO2NOOxid, ixHO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
# add in this HO2..
accumulate_row(a_row,hourly_tot_ho2p)

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("Total  HO2 Prod") 
hourly_diagram_values.append([e for e in hourly_tot_ho2p])
# daily total HO2 produced 
daily_tot_ho2p = sum(hourly_tot_ho2p)

daily_diagram_values.append(daily_tot_ho2p)

jj += 1

# now start the losses of HO2...
section_labels.append("HO2+(O3,rads) Loss") 
copyhours(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iHO2 ), a_row)
# these are negative values, chg sign
b_row = [-e for e in a_row]
hourly_diagram_values.append([e for e in b_row])
#   and start a new accumulator for HO2 loss
hourly_tot_ho2l = [e for e in b_row]

daily_diagram_values.append(sum(b_row))

jj += 1

section_labels.append("  O3   Loss") 
copyhours(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iO3  ), a_row)
# these are negative values, chg sign
b_row = [-e for e in a_row]
hourly_diagram_values.append([e for e in b_row])

daily_diagram_values.append(sum(b_row))

jj += 1

section_labels.append("  XO2  Loss") 
copyhours(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iXO2  ), a_row)
# these are negative values, chg sign
b_row = [-e for e in a_row]
hourly_diagram_values.append([e for e in b_row])

daily_diagram_values.append(sum(b_row))

jj += 1

section_labels.append("  C2O3 Loss") 
copyhours(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iC2O3  ), a_row)
# these are negative values, chg sign
b_row = [-e for e in a_row]
hourly_diagram_values.append([e for e in b_row])

daily_diagram_values.append(sum(b_row))

jj += 1

section_labels.append("  H2O2 Prod") 
copyhours(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iH2O2 ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  -OOX Prod") 
copyhours(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iOOX ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  OH   Prod") 
copyhours(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iOH ), a_row)
hourly_diagram_values.append([e for e in a_row])
# start a new accumulator for OH produced 
hourly_tot_ohp = [e for e in a_row]

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("HO2+OH Loss")
copyhours(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iH2O ), a_row)
# these are negative values, chg sign
b_row = [-e for e in a_row]
hourly_diagram_values.append([e for e in b_row])
# add this HO2 loss to total HO2 loss so far..
accumulate_row(b_row,hourly_tot_ho2l)

daily_diagram_values.append(sum(b_row))

jj += 1

section_labels.append("HO2+NO Loss")
copyhours(hourly_net_rxn_masses, i2j(n_HO2NOOxid, iHO2 ), a_row)
# these are negative values, chg sign
b_row = [-e for e in a_row]
hourly_diagram_values.append([e for e in b_row])
# add this HO2 loss to total HO2 loss so far..
accumulate_row(b_row,hourly_tot_ho2l)

daily_diagram_values.append(sum(b_row))

jj += 1

section_labels.append("  NO2 prod") 
copyhours(hourly_net_rxn_masses, i2j(n_HO2NOOxid, iNO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
# save the only NO2 produced in this section...
daily_no2_prod = sum(a_row)

daily_diagram_values.append(daily_no2_prod)

# calculate the NO2 produced per HO2 produced
hourly_no2_prod_ho2 = [ n/c for (n,c) in zip(a_row,hourly_tot_ho2p)]
daily_no2_prod_ho2 = daily_no2_prod/daily_tot_ho2p

jj += 1

section_labels.append("  OH prod") 
copyhours(hourly_net_rxn_masses, i2j(n_HO2NOOxid, iOH ), a_row)
hourly_diagram_values.append([e for e in a_row])
# add the OH produced via HO2+NO
accumulate_row(a_row,hourly_tot_ohp)

daily_diagram_values.append(sum(a_row))

# calculate the OH produced per HO2 prod
hourly_oh_prod_ho2 = [ o/c for (o,c) in zip(hourly_tot_ohp,hourly_tot_ho2p)]
daily_oh_prod = sum(hourly_tot_ohp)
daily_oh_prod_ho2 = daily_oh_prod/daily_tot_ho2p 

jj += 1

section_labels.append("Total HO2 Loss") 
hourly_diagram_values.append([e for e in hourly_tot_ho2l])
daily_diagram_values.append(sum(hourly_tot_ho2l))

jj += 1

section_labels.append("Total OH Prod") 
hourly_diagram_values.append([e for e in hourly_tot_ohp])
# calculate the daily total OH produced
daily_tot_ohp = sum(hourly_tot_ohp)
daily_diagram_values.append(daily_tot_ohp)

jj += 1

section_labels.append("NO2 Prod / HO2 Prod") 
hourly_diagram_values.append([e for e in hourly_no2_prod_ho2])
daily_diagram_values.append(daily_no2_prod_ho2)

jj += 1

section_labels.append("OH  Prod / HO2 Prod") 
hourly_diagram_values.append([e for e in hourly_oh_prod_ho2])
daily_diagram_values.append(daily_oh_prod_ho2)

jj += 1


# how many elements in vector are needed.
num_rows_in_section = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_diagram_row  = jj


#-------------------------------------------------------
### The Diagram Section for ** OH prod and loss **
this_diag_section = next_diag_section  
next_diag_section += 1

#    ... save the name and starting location of this net_rxn
diagram_sect_names.append("OH New, Prod, and Loss")
sec_id_num += 1

# compute radical thruputs for OH from XO2
hourly_oh_prod_xo2 = [ oh * ho2 for (oh,ho2) \
                       in zip(hourly_oh_prod_ho2,hourly_ho2_prod_xo2) ]
daily_oh_prod_xo2  = daily_oh_prod_ho2 * daily_ho2_prod_xo2

# compute radical thruputs for OH from C2O3 via XO2 
hourly_oh_prod_c2o3 = [ ohxo2 * c for (ohxo2,c) \
                       in zip(hourly_oh_prod_xo2,hourly_xo2_prod_c2o3) ]
daily_oh_prod_c2o3  = daily_oh_prod_xo2 * daily_xo2_prod_c2o3

this_jstart = jstart_next_diagram_row
diagram_sect_start.append(this_jstart)

jj = this_jstart

# accumulate the secondary new OH 
section_labels.append(" new C2O3     OH")
a_row = [ p*t for (p,t) in zip(hourly_newc2o3,hourly_oh_prod_c2o3)]
hourly_diagram_values.append([e for e in a_row])
# initialize total OH accumulator...
hourly_tot_newoh = [e for e in a_row]

daily_diagram_values.append(daily_newc2o3*daily_oh_prod_c2o3)

jj += 1

section_labels.append(" new XO2      OH")
a_row = [ p*t for (p,t) in zip(hourly_newxo2,hourly_oh_prod_xo2)]
hourly_diagram_values.append([e for e in a_row])
# add in this OH..
accumulate_row(a_row,hourly_tot_newoh)

daily_diagram_values.append(daily_newxo2*daily_oh_prod_xo2)

jj += 1

section_labels.append(" new HO2      OH")
a_row = [ p*t for (p,t) in zip(hourly_newho2,hourly_oh_prod_ho2)]
hourly_diagram_values.append([e for e in a_row])
# add in this OH..
accumulate_row(a_row,hourly_tot_newoh)

daily_diagram_values.append(daily_newho2*daily_oh_prod_ho2)

jj += 1

section_labels.append("Total sec new OH") 
# save the total secondary new OH
hourly_tot_sec_newoh = [e for e in hourly_tot_newoh]
daily_tot_sec_newoh  = sum(hourly_tot_sec_newoh)
hourly_diagram_values.append([e for e in hourly_tot_sec_newoh])

daily_diagram_values.append(daily_tot_sec_newoh)

jj += 1

section_labels.append(" O3/H2O2+hv   OH") 
copyhours(hourly_net_rxn_masses, i2j(n_O3hvrad, iOH  ), a_row)
hourly_diagram_values.append([e for e in a_row])
# add in this OH..
accumulate_row(a_row,hourly_tot_newoh)

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_O3hvrad, iOH )])

jj += 1

section_labels.append(" Ox+org       OH") 
copyhours(hourly_net_rxn_masses, i2j(n_OxOrgrad, iOH ), a_row)
hourly_diagram_values.append([e for e in a_row])
# add in this OH..
accumulate_row(a_row,hourly_tot_newoh)

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_OxOrgrad, iOH )])

jj += 1

section_labels.append("Total new OH")

daily_tot_newoh = sum(hourly_tot_newoh)
hourly_diagram_values.append([e for e in hourly_tot_newoh])

daily_diagram_values.append(daily_tot_newoh)

jj += 1

section_labels.append("recreated OH")
# first fetch the total OH that reacted in the OHOrgOxid net reaction
copyhours(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iOH ), a_row)
hourly_tot_oh_reacted = [-e for e in a_row]
daily_tot_oh_reacted  = sum(hourly_tot_oh_reacted)

# then compute the recreated OH by difference...
hourly_recreated_oh = [ t-n for (t,n) \
                         in zip(hourly_tot_oh_reacted, hourly_tot_newoh)] 

hourly_diagram_values.append([e for e in hourly_recreated_oh])

daily_diagram_values.append(sum(hourly_recreated_oh))

jj += 1

section_labels.append("Total OH reacted") 
copyhours(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iOH ), a_row)
hourly_tot_oh_reacted = [-e for e in a_row]
daily_tot_oh_reacted  = sum(hourly_tot_oh_reacted)
hourly_diagram_values.append([e for e in hourly_tot_oh_reacted])

daily_diagram_values.append(daily_tot_oh_reacted)

jj += 1

section_labels.append("OH chain length")
# first fetch the total OH that reacted in the OHOrgOxid net reaction
hourly_oh_chain = [ Q / q for (Q,q) \
                      in zip(hourly_tot_oh_reacted, hourly_tot_newoh)] 
daily_oh_chain  = daily_tot_oh_reacted / daily_tot_newoh

hourly_diagram_values.append([e for e in hourly_oh_chain])

daily_diagram_values.append(daily_oh_chain)

jj += 1

section_labels.append("OH Pr")
hourly_oh_pr = [ 1- (q/Q) for (Q,q) \
                      in zip(hourly_tot_oh_reacted, hourly_tot_newoh)] 
daily_oh_chain  = daily_tot_oh_reacted / daily_tot_newoh

hourly_diagram_values.append([e for e in hourly_oh_pr])

daily_diagram_values.append(1-(1/daily_oh_chain))

jj += 1

section_labels.append("Total VOC reacted") 
copyhours(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iVOC ), a_row)
hourly_tot_voc_reacted = [-e for e in a_row]
daily_tot_voc_reacted  = sum(hourly_tot_voc_reacted)
hourly_diagram_values.append([e for e in hourly_tot_voc_reacted])

daily_diagram_values.append(daily_tot_voc_reacted)

jj += 1

section_labels.append("  NO2 by C2O3") 
copyhours(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iNO2 ), a_row)
hourly_tot_no2_oxid = [e for e in a_row]
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  NO2 by XO2") 
copyhours(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iNO2 ), a_row)
accumulate_row(a_row,hourly_tot_no2_oxid)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  NO2 by HO2") 
copyhours(hourly_net_rxn_masses, i2j(n_HO2NOOxid, iNO2 ), a_row)
accumulate_row(a_row,hourly_tot_no2_oxid)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("Total NO2 Oxidized") 

hourly_diagram_values.append([e for e in hourly_tot_no2_oxid])

daily_tot_no2_oxid = sum(hourly_tot_no2_oxid)

daily_diagram_values.append(daily_tot_no2_oxid)

jj += 1

section_labels.append("NO2 per VOC") 
hourly_no2_voc = [ n/v for (n,v) \
                 in zip(hourly_tot_no2_oxid,hourly_tot_voc_reacted)]
hourly_diagram_values.append([e for e in hourly_no2_voc])

daily_no2_voc = daily_tot_no2_oxid / daily_tot_voc_reacted

daily_diagram_values.append(daily_no2_voc)

jj += 1

# how many elements in vector are needed.
num_rows_in_section = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_diagram_row  = jj

if options.output_format=="xml":
	from jaxml import XML_document
	fout.close()
	delim="\t"
	
	#Create the generator object
	xdoc = XML_document()
	xdoc.netreactions(version="0.1")
	xdoc._push()
	
	xdoc.delimiter(delim)
	xdoc.starthour(hour_number[0])
	xdoc.endhour(hour_number[-1])
	xdoc.summaries(irrfile=doc1.replace('"','').replace('\n',''))
	xdoc._push()
	
	kk = 0
	for i in range(0,len(hourly_diagram_values)):
		if i == diagram_sect_start[kk] :
			if kk!=0:
				xdoc._pop()
				xdoc._pop()
			
			xdoc._push()
			
			xdoc.summary()
			xdoc._push()
			xdoc.title(diagram_sect_names[kk])
	
			if kk < len(diagram_sect_start)-1:
				kk += 1
		xdoc.source()	
		xdoc.name(section_labels[i])
		xdoc.values(delim.join(map(str,hourly_diagram_values[i][0:num_hrs])), uom="ppb")
		xdoc.total(daily_diagram_values[i], uom='ppb')
		xdoc._pop()
		xdoc._push()
		
	xdoc._pop()
	xdoc._pop()
	xdoc._pop()
	xdoc._pop()
	
	xdoc.reactions(irrfile=doc1.replace('"','').replace('\n',''))
	xdoc._push()
	
	kk = 0
	for i in range(0,len(net_rxn_masses)):
		if i == net_rxn_jindex[kk] :
			xdoc._pop()
			xdoc._pop()
			xdoc._push()
			xdoc.reaction()
			xdoc._push()
			xdoc.title(net_rxn_names[kk])
	
			if kk < len(net_rxn_names)-1:
				kk += 1
		xdoc._pop()
		xdoc._push()
		xdoc.species()
		xdoc.name(SPC_Names[net_rxn_spcname[i]])
		dia_temp = []
		for t in range(0,len(hour_number)):
			dia_temp.append(hourly_net_rxn_masses[t][i])
		
		xdoc.values(delim.join(map(str,dia_temp)), uom="ppb")
		xdoc.total(daily_net_rxn_masses[i], uom='ppb')
	
	#write out values to file
	xdoc._output(output_filename)
else:
	# Output the hourly and total net reactions to file...
	print >>fout, "IRR file doc line was"
	print >>fout, doc1
	
	kk = 0
	print >>fout, "%-21s" % "Hour",
	for t in hour_number:
		print >>fout, "      %02d" % t,
	print >>fout,  "    Daily"	
	print >>fout
	print >>fout, "Summary Section Name"
	print >>fout, "Item             ppb"
	for i in range(0,len(hourly_diagram_values)):
		if i == diagram_sect_start[kk] :
			print >>fout
			print >>fout
			print >>fout, "%-20s" % diagram_sect_names[kk]
			print >>fout
			if kk < len(diagram_sect_start)-1:
				kk += 1
		print  >>fout, "%-20s " % section_labels[i],
		for t in range(0,num_hrs):
			print >>fout, "%8.4f" % (hourly_diagram_values[i][t]),
		print >>fout, "%8.4f" % daily_diagram_values[i]
	
	
	
	print >>fout
	print >>fout
	print >>fout
	
	
	# Output the hourly and total net reactions to file...
	print >>fout, "IRR file doc line was"
	print >>fout, doc1
	
	kk = 0
	print >>fout, "%-21s" % "Hour",
	for t in hour_number:
		print >>fout, "      %02d" % t,
	print >>fout,  "     Daily"	
	print >>fout
	print >>fout, "Net Reaction Name"
	print >>fout, "Species     ppb"
	for i in range(0,len(net_rxn_masses)):
		if i == net_rxn_jindex[kk] :
			print >>fout
			print >>fout
			print >>fout, "%-20s" % net_rxn_names[kk]
			print >>fout
			if kk < len(net_rxn_names)-1:
				kk += 1
		print  >>fout, "%-20s " % SPC_Names[net_rxn_spcname[i]],
		for t in range(0,len(hour_number)):
			print >>fout, "%8.4f" % (hourly_net_rxn_masses[t][i]),
		print >>fout, "%8.4f" % daily_net_rxn_masses[i]
	# finished the whole set of net rxn masses for all hours

fout.close()


print "F I N I S H E D"
