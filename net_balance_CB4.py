#!/usr/bin/env python2.3
__doc__ = \
"""
Make Files of Net Reactions and Chemical Parameters from IRR/IPR Data

net_balance.py [options] EXTFILE.EXT OUTFILE_ROOT_NAME

This version is hard-coded for CB4 Mechanism in CAMx

use -h or --help for details on options

"""

"""
This file official source is
$HeadURL$
"""
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"



import os.path
import sys
import re
	
from optparse import OptionParser
	

SCRIPT_ID_STRING = "Net_Balance_CB4.py " + ChangeDate[18:-2]

__version__ = "Version 1.5, (R%s), last changed by %s" % (RevisionNum[22:-2], ChangedBy[16:-2])



# create a cmd line parser...
usage  = "usage: %prog [options] EXTFILE OUTFILE"
version = "%%prog %s" % __version__
parser = OptionParser(usage=usage, version=version)

parser.add_option("-i", "--indir",
	dest="indir",
	default="",
	help="prepend input dir to inputfilename")

parser.add_option("-o", "--outdir",
	dest="outdir",
	default="",
	help="prepend output dir to outputfilename")

# extent of output
parser.add_option("-q", "--quiet",
	action="store_false", 
	dest="verbose", 
	default=True,
	help="don't print status messages to stdout")

# get options and arguments
(options, args) = parser.parse_args()

# requires filenames as argument
if len(args) != 2:
	parser.error("Need two arguments")

# process given file names to include optional in/out paths
inpath  = args[0]
outpath = args[1]
if options.indir :
	inpath = os.path.join(options.indir,inpath)
if options.outdir :
	outpath = os.path.join(options.outdir,outpath)
input_filename = os.path.abspath(inpath)
given_filename = os.path.abspath(outpath)
(root_output_filename, oldext) = os.path.splitext(given_filename)

# requires a valid input file
if not os.path.exists(input_filename):
	err_msg = "Input file '%s' does not exist" % input_filename
	parser.error(err_msg)


# document all these options...
if options.verbose:
	print
	print "======================================================"
	print "Script ",  SCRIPT_ID_STRING
	print "This is version ",  __version__
	print ""
	print "Reading from:"
	print "%s" % input_filename
	print ""
	print "Writing to  :"
	print "%s" % root_output_filename
	print ""

# ----- options are now processed


## ==========================================================
# create regular expressions to read the input file..
time_re  = re.compile('Time =[0-9]{6}', re.IGNORECASE)
irr_re   = re.compile('\{\s*\d+\}\s+\d+', re.IGNORECASE)
ipr_re   = re.compile('"\w+\s*"\s*\d+\s*', re.IGNORECASE)
split_re = re.compile('[ ]+')



## ==========================================================
# function to get an hour's worth of data..
def get_irr_data(f):
	"""
	f -- opened irr&ipr file
	return (time (as integer), ir_rates (as list), ip_rates (as list of lists)
	if time < 0, end of file
	"""
	# initialize ir and ip storage...
	ir_time  = -1
	ir_rates = [0,]  # account for 1-origin of rates
	ip_rates = []
	
	line = f.readline()
	if line[0] == '|':
		return (ir_time, ir_rates, ip_rates)
		
	# next line is first Time =
	if time_re.match(line) != None:   # if this is a 'Time =" line...
		time = line.split('=')[1][0:2] # pick 'HH' out of 'HH0000'
		ir_time = int(time)           # store hh
	else :
		print "ERROR:: in get_data read, did not find a time."
		sys.exit(1)
	
	# read in the !"Rxn No"    "Int Rate" header
	line = f.readline()
	
	# read the irr values
	while line[0] != ';' :
		line = f.readline()
		if irr_re.match(line) != None:    # if this is a { n} n.nnnn line ...
			ir_value = float(split_re.split(line.replace('{','').replace('}','').strip())[1])
			ir_rates.append(ir_value)
	
	# read in the ! Species      Initial conc. ... header
	line = f.readline()
	# 
	# ip_rates = [ [ip_values], [ip_values], [ip_values], ... ]
	#   order of sub-lists 'ip_values' is order of iProcess id Number
	#   for this hour.
	line = f.readline()
	while line[0] != ';' :
		if ipr_re.match(line) != None:
			ip_set = split_re.split(line[11:-1].strip())
			ip_values = [ float(v) for v in ip_set ]
			ip_rates.append(ip_values)
		line = f.readline()
	
	return (ir_time, ir_rates, ip_rates)






## ==========================================================
#  GLOBAL DATA FOR CAMX CB4 Mechanism 3 ...

# ------------ Constants for Species ------------------------

Net_spc_names = [
  'NO  ', 'NO2 ', 'O3  ', 'OLE ', 'PAN ', 'N2O5', 'PAR ', 'TOL ', 'XYL ',\
  'FORM', 'ALD2', 'ETH ', 'CRES', 'MGLY', 'OPEN', 'PNA ', 'CO  ', 'HONO',\
  'H2O2', 'HNO3', 'ISOP', 'MEOH', 'ETOH', 'CH4 ', 'O   ', 'O1D' , 'OH  ', 'HO2 ',\
  'NO3 ', 'C2O3', 'XO2 ', 'XO2N', 'NTR ', 'CRO ', 'ISPD', 'TO2 ', 'ROR ',\
  'SO2 ', 'xHO2', 'H2O ', '-OOX', 'VOC ', 'xFORM', 'xALD', 'xMGLY',\
  'O3_o', 'O__o', 'NO_o', 'NO2o', 'OH_o', 'NO3o', 'NA_o' ]

# Use xHO2 to track 'prompt HO2' and separate it from direct HO2.
# Use H2O  as a product in some reactions to track a term pathway
# Use -OOX to track the XO2+XO2 type reaction products.
# Use VOC  as a summing of individual HC species that reacted.
# Use xFORM and xALD to track the secondary FORM and ALD produced
# Use O3_o, O__o, OH_o, NO_o, NO2_o, NOzo to track the consumption of O3 
#    by early and late titration of NO to NO2 in NO cycle net reaction.


# integer indexing for CB4 species
#   ORDER IS IMPORTANT -- Must agree with CAMx output order
z = 0
iNO   =  z ; z += 1
iNO2  =  z ; z += 1
iO3   =  z ; z += 1
iOLE  =  z ; z += 1
iPAN  =  z ; z += 1
iN2O5 =  z ; z += 1
iPAR  =  z ; z += 1
iTOL  =  z ; z += 1
iXYL  =  z ; z += 1
iFORM =  z ; z += 1
iALD2 =  z ; z += 1 
iETH  =  z ; z += 1
iCRES =  z ; z += 1
iMGLY =  z ; z += 1
iOPEN =  z ; z += 1
iPNA  =  z ; z += 1
iCO   =  z ; z += 1
iHONO =  z ; z += 1
iH2O2 =  z ; z += 1
iHNO3 =  z ; z += 1
iISOP =  z ; z += 1
iMEOH =  z ; z += 1
iETOH =  z ; z += 1;  zp = z # zp is number of last species in physical process
iCH4  =  z ; z += 1
iO    =  z ; z += 1
iO1D  =  z ; z += 1
iOH   =  z ; z += 1
iHO2  =  z ; z += 1
iNO3  =  z ; z += 1
iC2O3 =  z ; z += 1
iXO2  =  z ; z += 1
iXO2N =  z ; z += 1
iNTR  =  z ; z += 1
iCRO  =  z ; z += 1
iISPD =  z ; z += 1
iTO2  =  z ; z += 1
iROR  =  z ; z += 1
iSO2  =  z ; z += 1
ixHO2 =  z ; z += 1
iH2O  =  z ; z += 1
iOOX  =  z ; z += 1
iVOC  =  z ; z += 1
ixFOR =  z ; z += 1
ixALD =  z ; z += 1
ixMGL =  z ; z += 1
iO3_o =  z ; z += 1
iO__o =  z ; z += 1
iNO_o =  z ; z += 1
iNO2o =  z ; z += 1
iOH_o =  z ; z += 1
iNO3o =  z ; z += 1
iNA_o =  z ; z += 1

num_pa_spc = z

# only the first 23 species (zp) are included in the physical processes...
num_phy_proc_spc = zp

# add some 'aggregate' species to processes using same slots...
iNOx  = iCH4
iNOy  = iO
iOx   = iO1D
iVOCm = iOH
iVOCC = iHO2
ikOHa = iNO3

#  NOx  = NO + NO2
#  NOy  = NOx + PAN + N2O5 + PNA + HONO + HNO3 
#  Ox   = O3 + NO2 + PAN + NO3 + 2*N2O5 + PNA + HNO3
#  VOCm = OLE + PAR + TOL + XYL + FORM + ALD2 + ETH + CRES + MGLY
#          + OPEN + ISOP + MEOH + ETOH
#  VOCC = 2*OLE + PAR + 7*TOL + 8*XYL + FORM + 2*ALD2 + 2*ETH + 7*CRES + 3*MGLY
#          + 7*OPEN + 5*ISOP + MEOH + 2*ETOH
#  kOHa = (VOC_i * k_OH_i)/VOCm

tot_num_phy_proc_spc = num_phy_proc_spc + 6

# provide labels for these physical process and aggregate species...
Phy_spc_names = [
  'NO  ', 'NO2 ', 'O3  ', 'OLE ', 'PAN ', 'N2O5', 'PAR ', 'TOL ', 'XYL ',\
  'FORM', 'ALD2', 'ETH ', 'CRES', 'MGLY', 'OPEN', 'PNA ', 'CO  ', 'HONO',\
  'H2O2', 'HNO3', 'ISOP', 'MEOH', 'ETOH', 'NOx ', 'NOy ', 'Ox  ', 'VOCm',\
  'VOCC', 'kOHa' ]


# provide constants needed to compute some of the aggregate variables

nox_spc  = [ iNO, iNO2 ]
noy_spc  = [ iNO, iNO2, iPAN, iN2O5, iPNA, iHONO, iHNO3 ]
ox_spc   = [ iO3, iNO2, iPAN, iN2O5, iPNA, iHNO3 ]

voc_comp_spc = [ iALD2, iFORM, iXYL, iTOL,
                 iISOP, iOLE, iETH, iETOH, iMEOH, iPAR ] 

num_composition_spc = len(voc_comp_spc)+1

splt_xport_spc = [ iNO, iNO2, iO3, iNOx, iNOy, iVOCm, iVOCC ]

num_splt_xport_spc = len(splt_xport_spc)+1

k_OH_spc = [ iOLE, iPAR, iTOL, iXYL, iFORM, iALD2, 
             iETH, iCRES, iMGLY, iISOP, iMEOH, iETOH ]
             
k_OH = [ 0.0 for p in range(num_phy_proc_spc) ]
scale_koh =  1.0E-4

# the k_OH at 298 from CB4 mech in per_ppm_min
k_OH[iPAR ] = scale_koh * 1.203E+3
k_OH[iMEOH] = scale_koh * 1.600E+3
k_OH[iETOH] = scale_koh * 4.300E+3
k_OH[iTOL ] = scale_koh * 9.150E+3
k_OH[iFORM] = scale_koh * 1.500E+4
k_OH[iETH ] = scale_koh * 1.192E+4
k_OH[iNO2 ] = scale_koh * 1.682E+4
k_OH[iALD2] = scale_koh * 2.400E+4
k_OH[iMGLY] = scale_koh * 2.600E+4
k_OH[iXYL ] = scale_koh * 3.620E+4
k_OH[iOLE ] = scale_koh * 4.200E+4
k_OH[iOPEN] = scale_koh * 4.400E+4
k_OH[iCRES] = scale_koh * 6.100E+4
k_OH[iISOP] = scale_koh * 1.476E+5

VOC_carb_num = [ 0.0 for p in range(num_phy_proc_spc) ]

VOC_carb_num[iOLE ] = 2
VOC_carb_num[iPAR ] = 1
VOC_carb_num[iTOL ] = 7
VOC_carb_num[iXYL ] = 8
VOC_carb_num[iFORM] = 1
VOC_carb_num[iALD2] = 2
VOC_carb_num[iETH ] = 2
VOC_carb_num[iCRES] = 7
VOC_carb_num[iMGLY] = 3
VOC_carb_num[iOPEN] = 5
VOC_carb_num[iISOP] = 5
VOC_carb_num[iMEOH] = 1
VOC_carb_num[iETOH] = 2


# ------------ Constants for Reactions ----------------------

# the number of irr lines in the input file
num_cb4_rxn = 97


# ------------ Constants for Processes ------------------------
# integer indexing for input file processes
#   ORDER IS IMPORTANT -- must agree with CAMx Output order
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

num_i_processes = iFinal + 1

# integer indexing for aggregated species processes
jInitial    =  0
jCarryOver  =  1
jEmission   =  2
jChemistry  =  3
jHor_Trans  =  4
jVer_Trans  =  5
jEntrain    =  6
jDepo       =  7
jFinal      =  8
jTempAdj    =  9
jFinalTadj  = 10

num_j_processes = jFinalTadj + 1

Proc_Names = [
 'Initial', 'CarryOver', 'Emissions', 'Chemistry', 'H_Trans', 'V_Trans', 
 'Entrain', 'Deposit', 'Final', 'T adj', 'Final, T adj']


split_gainloss_processes = [ jChemistry, jHor_Trans, jVer_Trans, jEntrain ]

hourly_sum_processes = [jInitial, jCarryOver, jEmission, jChemistry, 
                        jHor_Trans, jVer_Trans, jEntrain, jDepo, jTempAdj ]

daily_sum_processes  = [jInitial, jEmission, jChemistry, 
                        jHor_Trans, jVer_Trans, jEntrain, jDepo,
                        jTempAdj ]




## ==========================================================
# DATA STORAGE for Net Reactions, Summaries, Physical Processes

# Start with info about the net_rxns sets...
#   allocate vector of names of the net_rxns
net_rxn_names = []
#   this is indexed by next_net_rxn_set

# allocate a look up table between global species ID num (called iSPC)
#    and net_rxn_masses and net_rxn_spcname vector element 
#   [where SPC is any species name in the above global list].
net_rxn_species = []
#   this is indexed by next_net_rxn_set
#   each row of net_rxn_species will be a vector of len(num_pa_spc)
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

# allocate species processes masses array
# this is a list of lists of lists
#   outer list is by species, e.g.  [ [spc1],  [spc2],  [spc3], ...]
#   a species list is list of list by time, e.g. [ [ [proc], [proc], [proc]..]
#   a proc list is list of masses by process, e.g., [ [ [I,C,T,D,F], [I,C,T,D,F] ] ]
#   Index order [spc][time][process]
#  
species_process_masses = []
for s in range(tot_num_phy_proc_spc): 
	species_process_masses.append([])   # set up list for process species

#   accumulate the species_process_masses over all time.
total_species_process_masses = []
#   create a 0-filled initial vector of j-processes..
j_proc_null = [ 0.0 for p in range(num_j_processes) ]
for s in range(tot_num_phy_proc_spc): 
	# for each species, append a 0-valued j-processes vector
	total_species_process_masses.append([z for z in j_proc_null])
	# use incremented addition into each element to accumulate, 
	#  eg, total_species_process_masses[s][jInitial] += this_jInitial




## >>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
##  P A R T   O N E  : D E F I N E  T H E  N E T  R E A C T I O N S 


# @@@@@@@@ N E T   R E A C T I O N  G R O U P S  @@@@@@@@@
#
#  1) ** HONO+hv Radical Source  **   IdNum = n_HONOhvrad
#  2) ** Ald+hv  Radical Source  **   IdNum = n_Aldhvrad
#  3) ** Ox+Org  Radical Source  **   IdNum = n_OxOrgrad
#  4) ** NO3+Org Radical Source  **   IdNum = n_NO3Orgrad
#  5) ** OH  + (organic+NO2)     **   IdNum = n_OHOrgOxid
#  6) ** C2O3     + NO Oxidation **   IdNum = n_C2O3NOOxid
#  7) ** XO2/XO2N + NO Oxidation **   IdNum = n_XO2NOOxid
#  8) ** PAN Production          **   IdNum = n_PANProd
#  9) ** HO2 to OH via Radical   **   IdNum = n_HO2toOHrad
# 10) ** HO2      + NO Oxidation **   IdNum = n_HO2NOOxid
# 11) ** NO2 + hv O3 production  **   IdNum = n_NO2hvO3prod
# 12) ** NO2 Termination         **   IdNum = n_NO2term
# 13) ** O1D (hv) Radical Source **   IdNum = n_O3hvrad


## ==========================================================
# define the species storage for each set of net reactions

# init a counter for number of net reaction sets so far...
#    used to assign a name that refers to a given net_reaction.
nr_num = 0


## begin net_reaction storage allocation

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
jndx_net_rxn = [jNONE]*num_pa_spc
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

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

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
jndx_net_rxn = [jNONE]*num_pa_spc
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
jndx_net_rxn[ixFOR] = jj; jj += 1; indx_net_rxn.append(ixFOR)
jndx_net_rxn[ixALD] = jj; jj += 1; indx_net_rxn.append(ixALD)

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

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
jndx_net_rxn = [jNONE]*num_pa_spc
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
jndx_net_rxn[ixFOR] = jj; jj += 1; indx_net_rxn.append(ixFOR)
jndx_net_rxn[ixALD] = jj; jj += 1; indx_net_rxn.append(ixALD)
jndx_net_rxn[ixMGL] = jj; jj += 1; indx_net_rxn.append(ixMGL)

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

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
jndx_net_rxn = [jNONE]*num_pa_spc
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
jndx_net_rxn[ixFOR] = jj; jj += 1; indx_net_rxn.append(ixFOR)
jndx_net_rxn[ixALD] = jj; jj += 1; indx_net_rxn.append(ixALD)

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

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
jndx_net_rxn = [jNONE]*num_pa_spc
indx_net_rxn = []


# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
# use NA_o to catch mag of OH + HNO3 = NO3 reaction
jndx_net_rxn[iNA_o] = jj; jj += 1; indx_net_rxn.append(iNA_o)
# use H2O to catch mag of OH + HO2 = H2O + O2 reaction..
jndx_net_rxn[iH2O ] = jj; jj += 1; indx_net_rxn.append(iH2O )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iCH4 ] = jj; jj += 1; indx_net_rxn.append(iCH4 )
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH)
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)
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
jndx_net_rxn[ixFOR] = jj; jj += 1; indx_net_rxn.append(ixFOR)
jndx_net_rxn[ixALD] = jj; jj += 1; indx_net_rxn.append(ixALD)
jndx_net_rxn[ixMGL] = jj; jj += 1; indx_net_rxn.append(ixMGL)

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

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
jndx_net_rxn = [jNONE]*num_pa_spc
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
jndx_net_rxn[ixFOR] = jj; jj += 1; indx_net_rxn.append(ixFOR)

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

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
jndx_net_rxn = [jNONE]*num_pa_spc
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

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

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
jndx_net_rxn = [jNONE]*num_pa_spc
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

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

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
jndx_net_rxn = [jNONE]*num_pa_spc
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
jndx_net_rxn[ixFOR] = jj; jj += 1; indx_net_rxn.append(ixFOR)

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

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
jndx_net_rxn = [jNONE]*num_pa_spc
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

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

# add and init new elements to the net_rxn_masses vector the net rxn
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------
### The Net Rxns for ** NO2 +hv O3 Production **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('NO2+hv O3 Production')
n_NO2hvO3prod = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)


# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*num_pa_spc
indx_net_rxn = []


# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iNO  ] = jj; jj += 1; indx_net_rxn.append(iNO  )
jndx_net_rxn[iO3  ] = jj; jj += 1; indx_net_rxn.append(iO3  )
jndx_net_rxn[iO1D ] = jj; jj += 1; indx_net_rxn.append(iO1D )
jndx_net_rxn[iO   ] = jj; jj += 1; indx_net_rxn.append(iO   )
jndx_net_rxn[iNO3 ] = jj; jj += 1; indx_net_rxn.append(iNO3 )
jndx_net_rxn[iN2O5] = jj; jj += 1; indx_net_rxn.append(iN2O5)
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)

jndx_net_rxn[iO3_o] = jj; jj += 1; indx_net_rxn.append(iO3_o)
jndx_net_rxn[iNO_o] = jj; jj += 1; indx_net_rxn.append(iNO_o)
jndx_net_rxn[iNO2o] = jj; jj += 1; indx_net_rxn.append(iNO2o)
jndx_net_rxn[iO__o] = jj; jj += 1; indx_net_rxn.append(iO__o)
jndx_net_rxn[iOH_o] = jj; jj += 1; indx_net_rxn.append(iOH_o)
jndx_net_rxn[iNO3o] = jj; jj += 1; indx_net_rxn.append(iNO3o)
jndx_net_rxn[iNA_o] = jj; jj += 1; indx_net_rxn.append(iNA_o)

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

# add and init new elements to the net_rxn_masses vector for net rxn
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------
### The Net Rxns for ** NO2 Termination **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('NO2 Termination')
n_NO2term = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)


# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*num_pa_spc
indx_net_rxn = []


# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart
jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iPAN ] = jj; jj += 1; indx_net_rxn.append(iPAN )
jndx_net_rxn[iNTR ] = jj; jj += 1; indx_net_rxn.append(iNTR )


# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

# add and init new elements to the net_rxn_masses vector for species 
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


#-------------------------------------------------------
### The Net Rxns for ** O1D (hv) Radical Source **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1

#    ... save the name and starting location of this net_rxn
net_rxn_names.append('O1D (hv) radical source')
n_O3hvrad = nr_num  # provide an interger with a name for this nr
nr_num += 1

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)

# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*num_pa_spc
indx_net_rxn = []

# fill in the jindex at the iSPC locations for the net rxn
jj = this_jstart

jndx_net_rxn[iO1D ] = jj; jj += 1; indx_net_rxn.append(iO1D )
jndx_net_rxn[iH2O2] = jj; jj += 1; indx_net_rxn.append(iH2O2)
jndx_net_rxn[iOH  ] = jj; jj += 1; indx_net_rxn.append(iOH  )


# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# copy a new jndx vector for net_rxn to the net_rxn_species array
net_rxn_species.append([e for e in jndx_net_rxn])

# add and init new elements to the net_rxn_masses vector for species
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	daily_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)


## >>>>>> finished setting up all net reaction storage <<<<<<
num_net_rxn_sets   = this_net_rxn_set
max_j_net_rxn_set  = jstart_next_nr_set




## ==========================================================

## >>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
##  P A R T   T W O  : C A L C U L A T E  T H E  N E T  M A S S E S
##                     R E O R D E R  &  S U M  N E T  P R O D E S S E S

## ==========================================================
if options.verbose:
	print "Opening input file.."

fin = open(input_filename, 'r')

# first two lines are 'doc' lines giving data source
#   skip the quote error in doc line from file...
irrfile_docline = fin.readline()[1:-1]
# skip the second doc line...
line = fin.readline()

if options.verbose:
	print "IRR file doc line was"
	print irrfile_docline
	print


## read first hours' data... 
(time, ir, ip) = get_irr_data(fin);  


## ==========================================================
## big loop repeated for each time in file...
##
first_time = -1
while ( time > 0 ) :
	if options.verbose:
		print "Time %2d  hours" % time
	if len(ir) != num_cb4_rxn:
		print "ERROR: incorrect number of ir values!"
		sys.exit(1)
		
	if len(ip) != num_phy_proc_spc:
		print "ERROR: incorrect number of ip values!"
		sys.exit(1)
	
	hour_number.append(time)
	
	if first_time < 0:
		first_time = time
	
	#-------------------------------------------------------
	# ... process the ip data for this time...
	#-------------------------------------------------------
	# fill out accumulator vectors for aggregates...
	nox_j_proc  = [ 0.0 for p in range(num_j_processes)]
	noy_j_proc  = [ p   for p in nox_j_proc]
	ox_j_proc   = [ p   for p in nox_j_proc]
	vocm_j_proc = [ p   for p in nox_j_proc]
	vocc_j_proc = [ p   for p in nox_j_proc]
	koha_j_proc = [ p   for p in nox_j_proc]
	
	for s in range(num_phy_proc_spc):
		# use 'i' to refer to the full list of processes in file
		# use 'j' to refer to the lumped processes we are saving
		
		# get the i-processes read in for this species at this time...
		sp = ip[s]
		# create a j-process vector for this species at this time
		j_proc = []
		# start by appending the Initial concs
		#   this is true Initial the first time, CarryOver rest of time
		if time == first_time:
			j_proc.append(sp[iInitial])  # this is real initial
			j_proc.append(0.0)           # this is first carry over
		else:
			j_proc.append(0.0)           # second plus hrs are not initial
			j_proc.append(sp[iInitial])  #  but are carry over 
		
		# sum all of the emissions into a single value
		t_emiss = sp[iEmiss_Area] + sp[iEmiss_Pnt] + sp[iEmiss_PiG]
		j_proc.append(t_emiss)
		
		# append the chemistry
		j_proc.append(sp[iChemistry])
		
		# sum all the horizontal transport into a single value
		t_htrans = sp[iAdv_W] + sp[iAdv_E] + sp[iAdv_S] + sp[iAdv_N] \
		          +sp[iDif_W] + sp[iDif_E] + sp[iDif_S] + sp[iDif_N]
		j_proc.append(t_htrans)
		
		# sum the vertical transport into a single value
		t_vtrans = sp[iAdv_B] + sp[iAdv_T] + sp[iDif_B] + sp[iDif_T]
		j_proc.append(t_vtrans)
		
		# sum the entrainment and dilution into a single value
		j_proc.append(sp[iTrain]+sp[iDilut])
		
		# sum the wet and dry deposition into a single value
		j_proc.append(sp[iDep_D]+sp[iDep_W])
		
		# store the file's iFinal
		j_proc.append(sp[iFinal])
		
		# save a copy of iFinal in jTempAdj for temp_adjustment calculation
		j_proc.append(sp[iFinal])
		
		# save a copy of iFinal in jFinalTadj for temp adj calculation
		j_proc.append(sp[iFinal])
		
		# this is species s at this time by j-processes
		species_process_masses[s].append(j_proc)
		
		if s in nox_spc:
			for j in range(num_j_processes):
				 nox_j_proc[j] += j_proc[j]
		if s in noy_spc:
			for j in range(num_j_processes):
				 noy_j_proc[j] += j_proc[j]
		if s in ox_spc:
			for j in range(num_j_processes):
				 ox_j_proc[j] += j_proc[j]
		if s in k_OH_spc:
			for j in range(num_j_processes):
				vocm_j_proc[j] += j_proc[j]
				vocc_j_proc[j] += VOC_carb_num[s] * j_proc[j]
				koha_j_proc[j] += k_OH[s] * abs(j_proc[j])
		
		# now accumulate the values over time
		for jp in hourly_sum_processes :
			total_species_process_masses[s][jp] += j_proc[jp]
		# the Final will be done after read loop is over.
	
	for j in range(num_j_processes):
		if abs(vocm_j_proc[j]) > 0.0 :
			koha_j_proc[j] = koha_j_proc[j] / abs(vocm_j_proc[j])
		else:
			koha_j_proc[j] = k_OH[iNO2]
	
	species_process_masses[iNOx ].append( nox_j_proc)
	species_process_masses[iNOy ].append( nox_j_proc)
	species_process_masses[iOx  ].append(  ox_j_proc)
	species_process_masses[iVOCm].append(vocm_j_proc)
	species_process_masses[iVOCC].append(vocc_j_proc)
	species_process_masses[ikOHa].append(koha_j_proc)
	
	# accumulate relevant processes for the aggregates for each time.
	for jp in hourly_sum_processes  :
		total_species_process_masses[iNOx ][jp] +=  nox_j_proc[jp]
		total_species_process_masses[iNOy ][jp] +=  nox_j_proc[jp]
		total_species_process_masses[iOx  ][jp] +=   ox_j_proc[jp]
		total_species_process_masses[iVOCm][jp] += vocm_j_proc[jp]
		total_species_process_masses[iVOCC][jp] += vocc_j_proc[jp]
		total_species_process_masses[ikOHa][jp] += koha_j_proc[jp]
	# the aggregated Final will be done after read loop is over.
		
		
	#-------------------------------------------------------
	# ... process the ir data for this time...
	#-------------------------------------------------------
	#      ... zero out the net reaction masses vector
	for i in range(0,len(net_rxn_masses)):
		net_rxn_masses[i] = 0.0
	
	# initialize the counter for net reaction sets so
	# that it can be incremented in each reaction...
	kk = -1
	
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
	
	#   ... secondary aldehydes
	net_rxn_masses[i2j(kk,ixFOR)] =  ir[45]+0.9*ir[95]
	net_rxn_masses[i2j(kk,ixALD)] =  0.067*ir[95]
	
	
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
	
	#   ... secondary aldehydes
	net_rxn_masses[i2j(kk,ixFOR)] = 0.20*ir[56]+0.56*ir[75]\
										+0.74*ir[58]+ir[62]+0.70*ir[71]+0.60*ir[77]\
										+0.15*ir[93]
	net_rxn_masses[i2j(kk,ixALD)] = 0.63*ir[56]\
										+0.50*ir[58]       +0.03*ir[71]+0.15*ir[77]\
										+0.02*ir[93]
	net_rxn_masses[i2j(kk,ixMGL)] = 0.20*ir[71]+0.85*ir[93]
										
	
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
	
		#   ... secondary aldehydes
	net_rxn_masses[i2j(kk,ixFOR)] = ir[59]+0.282*ir[94]
	net_rxn_masses[i2j(kk,ixALD)] = ir[59]+0.80 *ir[78]+0.357*ir[94]
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** OH + (organic+NO2) **
	kk += 1  
	
	# { 26} OH + NO2  = HNO3    # use HNO3 for this HNO3
	# { 27} OH + HNO3 = NO3     # use NA_o for this HNO3
	
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
	#                    [ code 0.1*NTR as 0.1*NTR in XO2/XO2N set]
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
	#                    but tracked in 'NO2 Term' net reaction set.
	#   4) for ir[66], CRO, NTR is from NO2+rad not NO+RO2, omitted the CRO production.
	#                   but ir[68] production is tracked in 'NO2 Term' net reaction set.
	#
	#   5) for ir[64], TO2, the 0.9*NO2 is coded as 0.9*XO2 with 0.9*xHO2
	#                       and the 0.1*NTR is coded as 0.1*NTR in XO2N+NO net rxn
	#                       The NO is omitted as a a loss here but coded as
	#                       loss in the XO2 + NO net reaction below.
	
	
	#   ... the other reactant losses
	net_rxn_masses[i2j(kk,iNO2 )] = -ir[26]
	net_rxn_masses[i2j(kk,iNA_o)] = -ir[27]
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
	net_rxn_masses[i2j(kk,iOH  )] = -ir[26]-ir[27]-ir[90]-ir[12]-ir[35]-ir[36]\
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
	
	#   ... secondary aldehydes
	net_rxn_masses[i2j(kk,ixFOR)] = ir[84]+ir[70]+ir[51]+1.56*ir[61]+ir[57]
	net_rxn_masses[i2j(kk,ixALD)] = ir[85]+0.11*ir[52]+1.10*ir[53]+0.22*ir[61]\
										+ir[57]+0.273*ir[92]
	net_rxn_masses[i2j(kk,ixMGL)] = 0.80*ir[72]+0.168*ir[92]
	
	
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
	#   ... secondary aldehydes
	net_rxn_masses[i2j(kk,ixFOR)] = +ir[46] +2*ir[49] +0.79*ir[50]
	
		
	
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
	net_rxn_masses[i2j(kk,iXO2N)] = -ir[81] -2*ir[88] -ir[89]
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
	
	net_rxn_masses[i2j(kk,iNTR )] = +ir[81]+0.1*ir[64]
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
	#   ... secondary aldehydes
	net_rxn_masses[i2j(kk,ixFOR)] = +0.79*ir[50]
	
	
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
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** NO2+hv O3 Production **
	kk += 1  
	
	# {  1} NO2     = NO+O
	# {  2} O       = O3
	# {  3} O3+NO   = NO2
	# {  4} O+NO2   = NO
	# {  5} O+NO2   = NO3
	# {  6} O+NO    = NO2
	# {  7} O3+NO2  = NO3
	# {  8} O3      = O
	# {  9} O3      = O1D   # add back ir[11] cause in  ** O3+hv Rad Source **
	# { 10} O1D     = O
	#                       { 11} O1D+H2O = 2*OH   # in ** O3+hv Rad Source **
	#                       { 12} O3+OH   = HO2    # in ** OH + (Organic&NO2) **
	#                       { 13} O3+HO2  = OH     # in ** HO2 + OH via radical **
	# { 14} NO3     = 0.89*NO2 + 0.89*O + 0.11*NO
	# { 15} NO3+NO  = 2*NO2
	# { 16} NO3+NO2 = NO+NO2
	# { 17} NO3+NO2 = N2O5
	# { 18} N2O5+H2O= 2*HNO3
	# { 19} N2O5    = NO3+NO2
	# { 20} NO+NO   = 2*NO2
	
	#                       { 27} OH+HNO3 = NO3    # in ** OH + (organic+NO2) **

	
	## Notes:
	#  two conditions are expected:
	#    no net production of NO2 by other processes, and
	#      O3+NO == NO2 here       (this is 'old_O3')
	#    or
	#    HO2/RO2 + NO is producing NO2 and
	#      NO2+hv == O3 + NO here  (this is 'new O3')
	#    these cases can be distinguished by the sign of O3 in net reaction
	#
	#  In the former case, 
	#    save the mag of 'old O3' in net_reaction_mass[iO3_o] 
	#    zero out net_reaction_mass[iO3]
	#
	#  In the latter case,
	#    zero out the mag of net_reaction_mass[iO3_o]
	#  and expect to be exporting (ie producing)
	#    O3  - THE net production of ozone (still have O3+org losses, etc)
	#    O   - for reaction with organics in ** Ox+org new radicals **
	#    NO3 - for reaction with organics in ** NO3+org new radicals **
	#    OH  - for reaction with anything
	#
	
	
	
	net_rxn_masses[i2j(kk,iNO2 )] = +ir[3] +ir[6] +0.89*ir[14] +2*ir[15] +ir[16]\
	                                +ir[19] +2*ir[20]\
	                                -ir[1] -ir[4] -ir[5] -ir[7] -ir[16] -ir[17]
	                                
	net_rxn_masses[i2j(kk,iNO  )] = +ir[1] +ir[4] +0.11*ir[14] +ir[16]\
	                                -ir[3] -ir[6] -ir[15] -2*ir[20]
	
	net_rxn_masses[i2j(kk,iHNO3)] = +2*ir[18]
	
	net_rxn_masses[i2j(kk,iO   )] = +ir[1] +ir[8] +ir[10] +0.89*ir[14]\
	                                -ir[2] -ir[4] -ir[5] -ir[6]
	                                
	net_rxn_masses[i2j(kk,iO3  )] = +ir[2] -ir[3] -ir[7] -ir[8] -ir[9] +ir[11]
	
	net_rxn_masses[i2j(kk,iO1D )] = +ir[9] -ir[10]
	
	net_rxn_masses[i2j(kk,iNO3 )] = +ir[5] +ir[7] +ir[19] \
	                                -ir[14] -ir[15] -ir[16] -ir[17] 
	                                
	net_rxn_masses[i2j(kk,iN2O5)] = +ir[17] -ir[18] -ir[19]
	
	
	# prepare for the 'no net production of O3' case...
	net_rxn_masses[i2j(kk,iO3_o)] = 0.0 
	net_rxn_masses[i2j(kk,iNO_o)] = 0.0 
	net_rxn_masses[i2j(kk,iNO2o)] = 0.0 
	net_rxn_masses[i2j(kk,iO__o)] = 0.0 
	net_rxn_masses[i2j(kk,iOH_o)] = 0.0 
	net_rxn_masses[i2j(kk,iNO3o)] = 0.0 
	net_rxn_masses[i2j(kk,iNA_o)] = 0.0 
	
	if net_rxn_masses[i2j(kk,iO3  )] < 0 :
		# this is old O3 titrating NO and old O3 is oxidizing NO2
		#   save in special variables and zero out regular variables
		net_rxn_masses[i2j(kk,iO3_o)] = net_rxn_masses[i2j(kk,iO3  )]
		net_rxn_masses[i2j(kk,iNO_o)] = net_rxn_masses[i2j(kk,iNO  )]
		net_rxn_masses[i2j(kk,iNO2o)] = net_rxn_masses[i2j(kk,iNO2 )]
		net_rxn_masses[i2j(kk,iO__o)] = net_rxn_masses[i2j(kk,iO   )]
		net_rxn_masses[i2j(kk,iOH_o)] = net_rxn_masses[i2j(kk,iOH  )]
		net_rxn_masses[i2j(kk,iNO3o)] = net_rxn_masses[i2j(kk,iNO3 )]
		net_rxn_masses[i2j(kk,iNA_o)] = net_rxn_masses[i2j(kk,iHNO3)]
		
		net_rxn_masses[i2j(kk,iO3  )] = 0.0
		net_rxn_masses[i2j(kk,iNO  )] = 0.0
		net_rxn_masses[i2j(kk,iNO2 )] = 0.0
		net_rxn_masses[i2j(kk,iO   )] = 0.0
		net_rxn_masses[i2j(kk,iOH  )] = 0.0
		net_rxn_masses[i2j(kk,iNO3 )] = 0.0
		net_rxn_masses[i2j(kk,iN2O5)] = 0.0
		net_rxn_masses[i2j(kk,iHNO3)] = 0.0
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** NO2 Termination **  
	kk += 1  
	
	# { 26} NO2+OH=HNO3
	# { 47} C2O3+NO2=PAN
	# { 48} PAN=C2O3+NO2
	# { 55} ROR+NO2=NTR
	# { 68} CRO+NO2=NTR
	# 
	
	net_rxn_masses[i2j(kk,iNO2 )] = -ir[26]-ir[47]+ir[48]\
										-ir[55]-ir[68]
	net_rxn_masses[i2j(kk,iHNO3)] = +ir[26]
	net_rxn_masses[i2j(kk,iPAN )] = +ir[47] -ir[48]
	net_rxn_masses[i2j(kk,iNTR )] = +ir[55] +ir[68] 
	# { 64} TO2+NO=0.9*NO2+0.9*HO2+0.9*OPEN+0.1*NTR
	#  Omit +0.1*ir[64] cause TO2+NO is in prop cycle not RO2+NO2 loss
	# { 96} NO2+ISOP=0.2*ISPD+0.8*NTR+XO2+0.8*HO2+0.2*NO+0.8*ALD2+2.4*PAR
	#  -0.8*ir[96]  +0.8*ir[96]
	
	
	#-------------------------------------------------------
	### The Net Rxns for ** O1D (hv) Radical Source **
	kk += 1  
	
	# new OH from O3+hv
	# { 11} O1D+H2O =2*OH
	# { 34} H2O2    =2*OH
	
	# { 27} OH+HNO3 = NO3

	net_rxn_masses[i2j(kk,iO1D )] = -ir[11]
	net_rxn_masses[i2j(kk,iH2O2)] = -ir[34]
	net_rxn_masses[i2j(kk,iOH  )] = 2.0*ir[11] -ir[27] + 2.0*ir[34]
	
	
	# <<<add more net reactions assignments here....>>>
	
	
	
	
	# save these hourly results for later output....
	hourly_net_rxn_masses.append([e for e in net_rxn_masses])
	
	# accumulate this hour's net masses into a daily total..
	for i in range(0,len(net_rxn_masses)):
		daily_net_rxn_masses[i] += net_rxn_masses[i]
	
	# all reactions computed, read in a new set of ir and go back to top
	(time, ir, ip) = get_irr_data(fin);
	
	### end of the hour loop 


######### repeat to here for each hour's worth of net reactions ######


# Complete the 'end of loop' calculations....
# compute T adj between hours and correct jFinal 
# insert the last time jFinal into the total_species_process_masses
for s in range(tot_num_phy_proc_spc):
	# adjust each final values for temperature changes to next initial value
	for t in range(1,len(hour_number)):
		species_process_masses[s][t-1][jTempAdj] -= \
			species_process_masses[s][t][jCarryOver]
		species_process_masses[s][t-1][jTempAdj] *= -1.0
		species_process_masses[s][t-1][jFinalTadj] += \
			species_process_masses[s][t-1][jTempAdj]
	# do last hour ...
	species_process_masses[s][-1][jTempAdj] = 0.0
	
	# update the Final and FinalTadj for each species
	total_species_process_masses[s][jFinal] = \
	      species_process_masses[s][-1][jFinal]
	total_species_process_masses[s][jFinalTadj] = \
	      species_process_masses[s][-1][jFinalTadj]
	# compute the sum of the TempAdj
	total_species_process_masses[s][jTempAdj] = 0.0
	for t in range(len(hour_number)):
		total_species_process_masses[s][jTempAdj] += \
			species_process_masses[s][t][jTempAdj]
	
	# compute the averages for jCarryOver
	total_species_process_masses[s][jCarryOver] /= len(hour_number)
	      
total_species_process_masses[iNOx ][jFinal] = \
              species_process_masses[iNOx ][-1][jFinal]
total_species_process_masses[iNOy ][jFinal] = \
              species_process_masses[iNOy ][-1][jFinal]
total_species_process_masses[iOx  ][jFinal] = \
              species_process_masses[iOx  ][-1][jFinal]
total_species_process_masses[iVOCm][jFinal] = \
              species_process_masses[iVOCm][-1][jFinal]
total_species_process_masses[iVOCC][jFinal] = \
              species_process_masses[iVOCC][-1][jFinal]
total_species_process_masses[ikOHa][jFinal] = \
              species_process_masses[ikOHa][-1][jFinal]




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

def copyhours_net(hourlynet = [], elem = 0, row = []):
	"copy same element at each hour to an hourly vector of values"
	
	if len(row) == num_hrs :
		for t in range(0,num_hrs):
			row[t] = hourlynet[t][elem]
	else :
		raise IndexError
	
def copyhours_phy(hourlyphy = [], spc = 0, proc = 0, row = []):
	"copy same element at each hour to an hourly vector of values"
	
	if len(row) == num_hrs :
		for t in range(0,num_hrs):
			row[t] = hourlyphy[spc][t][proc]
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
#   5) "New NO2, physical NO2 sources"
#   6) "Available NO2"
#   7) "Physical NO2 losses"
#   7) "Chemical NO2 losses"
#   8) "new NO, Physical NO sources"
#   9) "new O3 balance"




## begin Diagram Section storage allocation
# intermediate working vector, must be of correct len
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

section_labels.append("  Ald+hv -> n_C2O3") 
copyhours_net(hourly_net_rxn_masses, i2j(n_Aldhvrad, iC2O3 ), a_row)
hourly_diagram_values.append([e for e in a_row])
hourly_tot_c2o3p = [e for e in a_row]

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_Aldhvrad, iC2O3 )])

jj += 1

section_labels.append("  Ox+org -> n_C2O3") 
copyhours_net(hourly_net_rxn_masses, i2j(n_OxOrgrad, iC2O3 ), a_row)
hourly_diagram_values.append([e for e in a_row])
accumulate_row(a_row,hourly_tot_c2o3p)

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_OxOrgrad, iC2O3 )])

jj += 1

section_labels.append("  NO3+org-> n_C2O3") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO3Orgrad, iC2O3 ), a_row)
hourly_diagram_values.append([e for e in a_row])
accumulate_row(a_row,hourly_tot_c2o3p)

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_NO3Orgrad, iC2O3 )])

jj += 1

section_labels.append("Total  New C2O3") 
#save hourly and daily new C3O3 for later use with OH
hourly_newc2o3 = [e for e in hourly_tot_c2o3p]
daily_newc2o3  = sum(hourly_newc2o3)
hourly_diagram_values.append([e for e in hourly_newc2o3])

daily_diagram_values.append(daily_newc2o3)

jj += 1

section_labels.append("OH+Org C2O3 Prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iC2O3 ), a_row)
hourly_diagram_values.append([e for e in a_row])
accumulate_row(a_row,hourly_tot_c2o3p)

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_OHOrgOxid, iC2O3 )])

jj += 1

section_labels.append("Total  C2O3 Prod") 
hourly_diagram_values.append([e for e in hourly_tot_c2o3p])
# remember the total C2O3 produced by sum of hourly_tot_c2o3p
daily_tot_c2o3p = sum(hourly_tot_c2o3p)
daily_diagram_values.append(daily_tot_c2o3p)

jj += 1

section_labels.append("C2O3+NO2 PAN Prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_PANProd, iPAN ), a_row)
hourly_diagram_values.append([e for e in a_row])
# count PAN production as a C2O3 loss
tot_c2o3l = [-e for e in a_row]  # PAN plus C2O2 + NO2

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("C2O3+NO Loss")
copyhours_net(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iC2O3 ), a_row)
hourly_diagram_values.append([e for e in a_row])
# add rxn with NO as C2O3 loss
accumulate_row(a_row,tot_c2o3l)

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  NO2 prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iNO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
daily_no2_prod = sum(a_row)
daily_diagram_values.append(daily_no2_prod)

# remember the NO2 produced by C2O3
hourly_no2_prod_c2o3 = [ n/c for (n,c) in zip(a_row,hourly_tot_c2o3p)]
daily_no2_prod_c2o3 = daily_no2_prod/daily_tot_c2o3p

jj += 1

section_labels.append("  XO2 prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iXO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

# remember the XO2 produced by C2O3
hourly_xo2_prod_c2o3 = [ x/c for (x,c) in zip(a_row,hourly_tot_c2o3p)]
daily_xo2_prod_c2o3 = sum(a_row)/daily_tot_c2o3p

jj += 1

section_labels.append("  OOX prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iOOX ), a_row)
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

section_labels.append("  Ald+hv -> n_XO2") 
copyhours_net(hourly_net_rxn_masses, i2j(n_Aldhvrad, iXO2 ), a_row)
# start acculating the sources of XO2
hourly_tot_xo2p = [e for e in a_row]
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_Aldhvrad, iXO2 )])

jj += 1

section_labels.append("  Ox+org -> n_XO2") 
copyhours_net(hourly_net_rxn_masses, i2j(n_OxOrgrad, iXO2 ), a_row)
accumulate_row(a_row,hourly_tot_xo2p)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_OxOrgrad, iXO2 )])

jj += 1

section_labels.append("  NO3+org-> n_XO2") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO3Orgrad, iXO2 ), a_row)
accumulate_row(a_row,hourly_tot_xo2p)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_NO3Orgrad, iXO2 )])

jj += 1

section_labels.append("Total  New XO2") 
# save the new XO2
hourly_newxo2 = [e for e in hourly_tot_xo2p]
daily_newxo2  = sum(hourly_newxo2)
hourly_diagram_values.append([e for e in hourly_newxo2])

daily_diagram_values.append(daily_newxo2)

jj += 1

section_labels.append("  C2O3   XO2 Prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iXO2 ), a_row)
accumulate_row(a_row,hourly_tot_xo2p)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  OH+Org XO2 Prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iXO2 ), a_row)
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

section_labels.append("  Ox+org -> n_XO2N") 
copyhours_net(hourly_net_rxn_masses, i2j(n_OxOrgrad, iXO2N ), a_row)
# accumulate the XO2N 
hourly_tot_xo2n = [e for e in a_row]
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_OxOrgrad, iXO2N )])

jj += 1

section_labels.append("  NO3+org-> n_XO2N") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO3Orgrad, iXO2N ), a_row)
accumulate_row(a_row,hourly_tot_xo2n)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_NO3Orgrad, iXO2N )])

jj += 1

section_labels.append("Total  New XO2N") 
hourly_diagram_values.append([e for e in hourly_tot_xo2n])

daily_diagram_values.append(sum(hourly_tot_xo2n))

jj += 1

section_labels.append("  OH+Org XO2N Prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iXO2N ), a_row)
accumulate_row(a_row,hourly_tot_xo2n)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("Total  XO2N Prod") 
hourly_diagram_values.append([e for e in hourly_tot_xo2n])
daily_diagram_values.append(sum(hourly_tot_xo2n))

jj += 1

section_labels.append("Total  XO2+XO2N Prod")
accumulate_row(hourly_tot_xo2p,hourly_tot_xo2n)
hourly_diagram_values.append([e for e in hourly_tot_xo2n])
daily_diagram_values.append(sum(hourly_tot_xo2n))

jj += 1

section_labels.append("  XO2  + NO Loss")
copyhours_net(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iXO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
tot_rows = [e for e in a_row]
daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  XO2N + NO Loss")
copyhours_net(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iXO2N ), a_row)
hourly_diagram_values.append([e for e in a_row])
accumulate_row(a_row,tot_rows)
daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("Total XO2+XO2N Loss") 
hourly_diagram_values.append([e for e in tot_rows])
daily_diagram_values.append(sum(tot_rows))

jj += 1

section_labels.append("  NO  loss") 
copyhours_net(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iNO ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  NO2 prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iNO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
daily_no2_prod = sum(a_row)

daily_diagram_values.append(daily_no2_prod)

# remember NO2 produced per XO2 produced
no2_prod_xo2 = [ n/x for (n,x) in zip(a_row,hourly_tot_xo2p)]
daily_no2_prod_xo2 = daily_no2_prod/daily_tot_xo2p
jj += 1

section_labels.append("  HO2 prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_XO2NOOxid, ixHO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
daily_ho2_prod = sum(a_row)

daily_diagram_values.append(daily_ho2_prod)

# remember HO2 produced per XO2 produced
hourly_ho2_prod_xo2 = [ h/x for (h,x) in zip(a_row,hourly_tot_xo2p)]
daily_ho2_prod_xo2 = daily_ho2_prod/daily_tot_xo2p

jj += 1

section_labels.append("  NTR prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iNTR ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append(" *OOX prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iOOX ), a_row)
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

section_labels.append("  Ald+hv -> n_HO2") 
copyhours_net(hourly_net_rxn_masses, i2j(n_Aldhvrad, iHO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
# initialize total ho2 accumulator...
hourly_tot_ho2p = [e for e in a_row]

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_Aldhvrad, iHO2 )])

jj += 1

section_labels.append("  Ox+org -> n_HO2") 
copyhours_net(hourly_net_rxn_masses, i2j(n_OxOrgrad, iHO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
# add in this HO2..
accumulate_row(a_row,hourly_tot_ho2p)

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_OxOrgrad, iHO2 )])

jj += 1

section_labels.append("Total New HO2 Prod") 
# save total new HO2
hourly_newho2 = [e for e in hourly_tot_ho2p]
daily_newho2  = sum(hourly_newho2)
hourly_diagram_values.append([e for e in hourly_newho2])

daily_diagram_values.append(daily_newho2)

jj += 1

section_labels.append("  OH+Org HO2 Prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iHO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
# add in this HO2..
accumulate_row(a_row,hourly_tot_ho2p)

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  XO2+NO HO2 Prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_XO2NOOxid, ixHO2 ), a_row)
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
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iHO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
#   and start a new accumulator for HO2 loss
hourly_tot_ho2l = [e for e in a_row]

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  O3   Loss") 
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iO3  ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  XO2  Loss") 
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iXO2  ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  C2O3 Loss") 
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iC2O3  ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  H2O2 Prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iH2O2 ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  -OOX Prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iOOX ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  OH   Prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iOH ), a_row)
hourly_diagram_values.append([e for e in a_row])
# start a new accumulator for OH produced 
hourly_tot_ohp = [e for e in a_row]

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("HO2+OH Loss")
copyhours_net(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iH2O ), a_row)
hourly_diagram_values.append([e for e in a_row])
# add this HO2 loss to total HO2 loss so far..
accumulate_row(a_row,hourly_tot_ho2l)

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("HO2+NO Loss")
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2NOOxid, iHO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
# add this HO2 loss to total HO2 loss so far..
accumulate_row(a_row,hourly_tot_ho2l)

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  NO2  Prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2NOOxid, iNO2 ), a_row)
hourly_diagram_values.append([e for e in a_row])
# save the only NO2 produced in this section...
daily_no2_prod = sum(a_row)

daily_diagram_values.append(daily_no2_prod)

# calculate the NO2 produced per HO2 produced
hourly_no2_prod_ho2 = [ n/c for (n,c) in zip(a_row,hourly_tot_ho2p)]
daily_no2_prod_ho2 = daily_no2_prod/daily_tot_ho2p

jj += 1

section_labels.append("  OH   Prod") 
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2NOOxid, iOH ), a_row)
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
section_labels.append("  new C2O3  -> n_OH")
a_row = [ p*t for (p,t) in zip(hourly_newc2o3,hourly_oh_prod_c2o3)]
hourly_diagram_values.append([e for e in a_row])
# initialize total OH accumulator...
hourly_tot_newoh = [e for e in a_row]

daily_diagram_values.append(daily_newc2o3*daily_oh_prod_c2o3)

jj += 1

section_labels.append("  new XO2   -> n_OH")
a_row = [ p*t for (p,t) in zip(hourly_newxo2,hourly_oh_prod_xo2)]
hourly_diagram_values.append([e for e in a_row])
# add in this OH..
accumulate_row(a_row,hourly_tot_newoh)

daily_diagram_values.append(daily_newxo2*daily_oh_prod_xo2)

jj += 1

section_labels.append("  new HO2   -> n_OH")
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

section_labels.append("  (O3&H2O2)+hv->n_OH") 
# this is directly from O3+hv-->2 OH, so it includes
#   O3_o + hv already.
copyhours_net(hourly_net_rxn_masses, i2j(n_O3hvrad, iOH  ), a_row)
hourly_diagram_values.append([e for e in a_row])
# add in this OH..
accumulate_row(a_row,hourly_tot_newoh)

daily_diagram_values.append(daily_net_rxn_masses[i2j(n_O3hvrad, iOH )])

jj += 1

section_labels.append("  Ox+org      ->n_OH") 
copyhours_net(hourly_net_rxn_masses, i2j(n_OxOrgrad, iOH ), a_row)
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

section_labels.append("recreated OH (diff)")
# first fetch the total OH that reacted in the OHOrgOxid net reaction
copyhours_net(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iOH ), a_row)
hourly_tot_oh_reacted = [e for e in a_row]
daily_tot_oh_reacted  = sum(hourly_tot_oh_reacted)

# then compute the recreated OH by difference...
hourly_recreated_oh = [ -(t)-n for (t,n) \
                         in zip(hourly_tot_oh_reacted, hourly_tot_newoh)] 

hourly_diagram_values.append([e for e in hourly_recreated_oh])

daily_diagram_values.append(sum(hourly_recreated_oh))

jj += 1

section_labels.append("Total OH reacted") 
hourly_diagram_values.append([e for e in hourly_tot_oh_reacted])

daily_diagram_values.append(daily_tot_oh_reacted)

jj += 1

# Chain Lenght is Q/q = 1 / (1 - Pn)
section_labels.append("OH chain length")
# first fetch the total OH that reacted in the OHOrgOxid net reaction
hourly_oh_chain = [ -Q / q for (Q,q) \
                      in zip(hourly_tot_oh_reacted, hourly_tot_newoh)] 
daily_oh_chain  = -daily_tot_oh_reacted / daily_tot_newoh

hourly_diagram_values.append([e for e in hourly_oh_chain])

daily_diagram_values.append(daily_oh_chain)

jj += 1

# Propagation factor is recreated OH over total OH reacted
#  Pr = 1 - (q/Q)
# 
section_labels.append("OH P_r")
hourly_oh_pr = [ 1- (q/-Q) for (Q,q) \
                      in zip(hourly_tot_oh_reacted, hourly_tot_newoh)] 
daily_oh_chain  = 1- ( daily_tot_newoh / -daily_tot_oh_reacted)

hourly_diagram_values.append([e for e in hourly_oh_pr])

daily_diagram_values.append(daily_oh_chain)

jj += 1

section_labels.append("Total VOC reacted") 
copyhours_net(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iVOC ), a_row)
hourly_tot_voc_reacted = [e for e in a_row]
daily_tot_voc_reacted  = sum(hourly_tot_voc_reacted)
hourly_diagram_values.append([e for e in hourly_tot_voc_reacted])

daily_diagram_values.append(daily_tot_voc_reacted)

jj += 1

section_labels.append("  NO  by C2O3") 
copyhours_net(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iNO ), a_row)
hourly_tot_no_oxid = [e for e in a_row]
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  NO  by XO2") 
copyhours_net(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iNO ), a_row)
accumulate_row(a_row,hourly_tot_no_oxid)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  NO  by HO2") 
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2NOOxid, iNO ), a_row)
accumulate_row(a_row,hourly_tot_no_oxid)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("Total NO  Oxidized") 

hourly_diagram_values.append([e for e in hourly_tot_no_oxid])

daily_tot_no_oxid = sum(hourly_tot_no_oxid)

daily_diagram_values.append(daily_tot_no_oxid)

jj += 1

section_labels.append("  NTR by XO2") 
copyhours_net(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iNTR ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  NO2 by C2O3") 
copyhours_net(hourly_net_rxn_masses, i2j(n_C2O3NOOxid, iNO2 ), a_row)
hourly_tot_no2_oxid = [e for e in a_row]
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  NO2 by XO2") 
copyhours_net(hourly_net_rxn_masses, i2j(n_XO2NOOxid, iNO2 ), a_row)
accumulate_row(a_row,hourly_tot_no2_oxid)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("  NO2 by HO2") 
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2NOOxid, iNO2 ), a_row)
accumulate_row(a_row,hourly_tot_no2_oxid)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("Total NO2 Produced") 

hourly_diagram_values.append([e for e in hourly_tot_no2_oxid])

daily_tot_no2_oxid = sum(hourly_tot_no2_oxid)

daily_diagram_values.append(daily_tot_no2_oxid)

jj += 1

section_labels.append("NO2 per VOC reacted") 
hourly_no2_voc = [ n/-v for (n,v) \
                 in zip(hourly_tot_no2_oxid,hourly_tot_voc_reacted)]
hourly_diagram_values.append([e for e in hourly_no2_voc])

daily_no2_voc = daily_tot_no2_oxid / -daily_tot_voc_reacted

daily_diagram_values.append(daily_no2_voc)

jj += 1

# how many elements in vector are needed.
num_rows_in_section = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_diagram_row  = jj


#-------------------------------------------------------
### The Diagram Section for ** new NO2, Src Outside Box **
this_diag_section = next_diag_section  
next_diag_section += 1

#    ... save the name and starting location of this net_rxn
diagram_sect_names.append("New NO2, Srcs Outside Box")
sec_id_num += 1

this_jstart = jstart_next_diagram_row
diagram_sect_start.append(this_jstart)

jj = this_jstart

#  Start with emissions
section_labels.append("  NO2 emissions")
copyhours_phy(species_process_masses, iNO2, jEmission, a_row)
hourly_total_new_no2 = [e for e in a_row ]
hourly_diagram_values.append([e for e in a_row])

t_value = total_species_process_masses[iNO2][jEmission]
daily_total_new_no2 = t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... add horiz trans as source
section_labels.append("  NO2 horiz trans")
copyhours_phy(species_process_masses, iNO2, jHor_Trans, a_row)
# only count horiz trans as new NO2 if positive value
for h in range(num_hrs):
	if a_row[h] < 0.0:
		a_row[h] = 0.0
accumulate_row(a_row,hourly_total_new_no2)
hourly_diagram_values.append([e for e in a_row ])

t_value = sum(a_row)
daily_total_new_no2 += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... add vert trans
section_labels.append("  NO2 vert  trans")
copyhours_phy(species_process_masses, iNO2, jVer_Trans, a_row)
# only count vert trans as new NO2 if positive value
for h in range(num_hrs):
	if a_row[h] < 0.0:
		a_row[h] = 0.0
accumulate_row(a_row,hourly_total_new_no2)
hourly_diagram_values.append([e for e in a_row ])

t_value = sum(a_row)
daily_total_new_no2 += abs(t_value)
daily_diagram_values.append(t_value)

jj += 1

#  ... add entrain/dil
section_labels.append("  NO2 entrain/dil")
copyhours_phy(species_process_masses, iNO2, jEntrain, a_row)
# only count entrain as new NO2 if positive value
for h in range(num_hrs):
	if a_row[h] < 0.0:
		a_row[h] = 0.0
accumulate_row(a_row,hourly_total_new_no2)
hourly_diagram_values.append([e for e in a_row ])

t_value = sum(a_row)
daily_total_new_no2 += t_value
daily_diagram_values.append(t_value)
jj += 1

#  ...sum
section_labels.append("Total new NO2")
hourly_diagram_values.append([e for e in hourly_total_new_no2])

daily_diagram_values.append(daily_total_new_no2)

jj += 1

# how many elements in vector are needed.
num_rows_in_section = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_diagram_row  = jj


#-------------------------------------------------------
### The Diagram Section for ** Total Available NO2 Inside Box **
this_diag_section = next_diag_section  
next_diag_section += 1

#    ... save the name and starting location of this net_rxn
diagram_sect_names.append("Total Available NO2 Inside Box")
sec_id_num += 1

this_jstart = jstart_next_diagram_row
diagram_sect_start.append(this_jstart)

jj = this_jstart

# accumulate the physical NO2 processes 
#  ... initial concentration
section_labels.append("  NO2 initial")
copyhours_phy(species_process_masses, iNO2, jInitial, a_row)
hourly_total_avail_no2 = [e for e in a_row ]
hourly_diagram_values.append([e for e in a_row ])

t_value = total_species_process_masses[iNO2][jInitial]
daily_total_avail_no2 = t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... carry over concentration
section_labels.append("  NO2 carry over")
copyhours_phy(species_process_masses, iNO2, jCarryOver, a_row)
hourly_no2_carryover = [e for e in a_row]
accumulate_row(a_row, hourly_total_avail_no2)
hourly_diagram_values.append([e for e in a_row ])

t_value = total_species_process_masses[iNO2][jCarryOver]
daily_total_avail_no2 += 0.0
daily_diagram_values.append(t_value)

jj += 1

#  ... NO+O3_o titration
section_labels.append("  NO2 by O3o+NO") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO2hvO3prod, iNO2o ), a_row)
accumulate_row(a_row, hourly_total_avail_no2)
hourly_diagram_values.append([e for e in a_row])

daily_total_avail_no2 = sum(a_row)
daily_diagram_values.append(sum(a_row))

jj += 1

#  ... add total new NO2 ...
section_labels.append("  new NO2")
accumulate_row(hourly_total_new_no2, hourly_total_avail_no2)
hourly_diagram_values.append([e for e in hourly_total_new_no2])

daily_total_avail_no2 += daily_total_new_no2
daily_diagram_values.append(daily_total_new_no2)

jj += 1

#   ... add the chemically produced NO2...
section_labels.append("  NO2 from NO oxid")
accumulate_row(hourly_tot_no2_oxid, hourly_total_avail_no2)
hourly_diagram_values.append([e for e in hourly_tot_no2_oxid])

daily_total_avail_no2 += daily_tot_no2_oxid
daily_diagram_values.append(daily_tot_no2_oxid)

jj += 1

#   ... show Total Available NO2 and its fraction from chemistry
section_labels.append("Total Available NO2")
hourly_diagram_values.append([e for e in hourly_total_avail_no2])

daily_diagram_values.append(daily_total_avail_no2)

jj += 1

#   ... fraction avail NO2 from chemistry
section_labels.append("Frac NO2 from chem")

hourly_frac_no2_chem = [0.0]*num_hrs
hourly_frac_no2_chem[0] = hourly_tot_no2_oxid[0]/hourly_total_avail_no2[0]
for t in range(1,num_hrs):
	hourly_frac_no2_chem[t] = \
		(hourly_frac_no2_chem[t-1]*hourly_no2_carryover[t] + \
		 hourly_tot_no2_oxid[t]) / hourly_total_avail_no2[t]
hourly_diagram_values.append([e for e in hourly_frac_no2_chem])

daily_frac_no2_chem = daily_tot_no2_oxid / daily_total_avail_no2

daily_diagram_values.append(daily_frac_no2_chem)

jj += 1

# how many elements in vector are needed.
num_rows_in_section = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_diagram_row  = jj


#-------------------------------------------------------
### The Diagram Section for ** physical NO2 losses **
this_diag_section = next_diag_section  
next_diag_section += 1

#    ... save the name and starting location of this net_rxn
diagram_sect_names.append("Physical NO2 Losses")
sec_id_num += 1

this_jstart = jstart_next_diagram_row
diagram_sect_start.append(this_jstart)

jj = this_jstart

# accumulate the physical NO2 processes 
#  ... final concentration (as a loss of available NO2)
section_labels.append("  NO2 Final")
copyhours_phy(species_process_masses, iNO2, jFinal, a_row)
hourly_tot_phy_no2_losses  = [-e for e in a_row ]
hourly_diagram_values.append([-e for e in a_row ])

t_value = total_species_process_masses[iNO2][jFinal]
daily_tot_phy_no2_losses =  -t_value
daily_diagram_values.append(-t_value)

jj += 1

#  ... deposition
section_labels.append("  NO2 deposition")
copyhours_phy(species_process_masses, iNO2, jDepo, a_row)
accumulate_row(a_row,hourly_tot_phy_no2_losses)

hourly_diagram_values.append([e for e in a_row ])

t_value = total_species_process_masses[iNO2][jDepo]
daily_tot_phy_no2_losses += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... horiz trans
section_labels.append("  NO2 horiz trans")
copyhours_phy(species_process_masses, iNO2, jHor_Trans, a_row)
# only count horiz trans a loss if neg value
for h in range(num_hrs):
	if a_row[h] > 0.0:
		a_row[h] = 0.0
accumulate_row(a_row,hourly_tot_phy_no2_losses)
hourly_diagram_values.append([e for e in a_row ])

t_value = sum(a_row)
daily_tot_phy_no2_losses += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... vert trans
section_labels.append("  NO2 vert  trans")
copyhours_phy(species_process_masses, iNO2, jVer_Trans, a_row)
# only count vert trans a loss if neg value
for h in range(num_hrs):
	if a_row[h] > 0.0:
		a_row[h] = 0.0
accumulate_row(a_row,hourly_tot_phy_no2_losses)
hourly_diagram_values.append([e for e in a_row ])

t_value = sum(a_row)
daily_tot_phy_no2_losses += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... entrain/dil
section_labels.append("  NO2 entrain/dil")
copyhours_phy(species_process_masses, iNO2, jEntrain, a_row)
# only count entrain a loss if neg value
for h in range(num_hrs):
	if a_row[h] > 0.0:
		a_row[h] = 0.0
accumulate_row(a_row,hourly_tot_phy_no2_losses)
hourly_diagram_values.append([e for e in a_row ])

t_value = sum(a_row)
daily_tot_phy_no2_losses += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ...sum
section_labels.append("Total NO2 phy losses")
hourly_diagram_values.append([e for e in hourly_tot_phy_no2_losses])

daily_diagram_values.append(daily_tot_phy_no2_losses)

jj += 1

#   ... compute the NO2 avail for chemistry as a diff...
section_labels.append("Total NO2 for Chem")

hourly_chem_sink_no2 = [ a+l for (a,l) in \
						zip(hourly_total_avail_no2, hourly_tot_phy_no2_losses)]
hourly_diagram_values.append([e for e in hourly_chem_sink_no2])

daily_chem_sink_no2 =sum(hourly_chem_sink_no2)
daily_diagram_values.append(daily_chem_sink_no2)

jj += 1

# how many elements in vector are needed.
num_rows_in_section = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_diagram_row  = jj


#-------------------------------------------------------
### The Diagram Section for ** chemical NO2 losses **
this_diag_section = next_diag_section  
next_diag_section += 1

#    ... save the name and starting location of this net_rxn
diagram_sect_names.append("Chemical NO2 Losses")
sec_id_num += 1

this_jstart = jstart_next_diagram_row
diagram_sect_start.append(this_jstart)

jj = this_jstart

section_labels.append("   HNO3 formation") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO2term, iHNO3 ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("   PAN  formation") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO2term, iPAN ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("   NTR  formation") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO2term, iNTR ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("Total NO2 Term Loss") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO2term, iNO2 ), a_row)
hourly_chem_no2_loss = [e for e in a_row ]
hourly_diagram_values.append([e for e in a_row])

daily_chem_no2_loss = sum(a_row)
daily_diagram_values.append(daily_chem_no2_loss)

jj += 1

section_labels.append("   NO2 Photolyzed") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO2hvO3prod, iNO2 ), a_row)
hourly_no2_phot = [e for e in a_row]
hourly_diagram_values.append(hourly_no2_phot)

daily_no2_phot = sum(a_row)
daily_diagram_values.append(daily_no2_phot)

jj += 1

section_labels.append("   NO  Produced") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO2hvO3prod, iNO  ), a_row)
hourly_no_prod = [e for e in a_row]
hourly_diagram_values.append(hourly_no_prod)

daily_no_prod = sum(a_row)
daily_diagram_values.append(daily_no_prod)

jj += 1

section_labels.append("   O3  Produced") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO2hvO3prod, iO3  ), a_row)
hourly_new_o3_prod = [e for e in a_row]
hourly_diagram_values.append(hourly_new_o3_prod)

daily_new_o3_prod = sum(a_row)
daily_diagram_values.append(daily_new_o3_prod)

jj += 1

section_labels.append("   NO3 Produced") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO2hvO3prod, iNO3  ), a_row)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

section_labels.append("Total NO2 Chem Loss") 
a_row = [t+p for (t,p) in \
						zip(hourly_chem_no2_loss, hourly_no2_phot) ]
hourly_tot_chem_no2_loss = [e for e in a_row ]
hourly_diagram_values.append(hourly_tot_chem_no2_loss)

daily_tot_chem_no2_loss = sum(hourly_tot_chem_no2_loss)
daily_diagram_values.append(daily_tot_chem_no2_loss)

jj += 1

hourly_o3_yield = [0.0] * num_hrs
section_labels.append("O3 yld per NO2 phot")
for h in range(num_hrs):
	if hourly_no2_phot[h] < 0.0 :
		hourly_o3_yield[h] = hourly_new_o3_prod[h]/abs(hourly_no2_phot[h])
hourly_diagram_values.append(hourly_o3_yield)

daily_o3_yield = daily_new_o3_prod/abs(daily_no2_phot)
daily_diagram_values.append(daily_o3_yield)

jj += 1

section_labels.append("NO yld per NO2 avail") 
hourly_no_yield = [ n/a for (n,a) in \
				zip(hourly_no_prod, hourly_total_avail_no2) ]
hourly_diagram_values.append(hourly_no_yield)

daily_no_yield = daily_no_prod/daily_total_avail_no2
daily_diagram_values.append(daily_no_yield)

jj += 1

section_labels.append("new NO frm new NO2") 
hourly_new_no_frm_no2 = [ n * y for (n,y) in \
								zip(hourly_total_new_no2, hourly_no_yield) ]
hourly_diagram_values.append([e for e in hourly_new_no_frm_no2])

daily_no_frm_no2 = daily_total_new_no2 * daily_no_yield
daily_diagram_values.append(daily_no_frm_no2)

jj += 1

section_labels.append("recreated NO") 
hourly_recreated_no = [ p - n for (p,n) in \
								zip(hourly_no_prod, hourly_new_no_frm_no2) ]
hourly_diagram_values.append([e for e in hourly_recreated_no])

daily_recreated_no = daily_no_prod - daily_no_frm_no2
daily_diagram_values.append(daily_recreated_no)

jj += 1

# how many elements in vector are needed.
num_rows_in_section = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_diagram_row  = jj


#-------------------------------------------------------
### The Diagram Section for ** new NO, Physical Sources **
this_diag_section = next_diag_section  
next_diag_section += 1

#    ... save the name and starting location of this net_rxn
diagram_sect_names.append("New NO, Physical NO sources")
sec_id_num += 1

this_jstart = jstart_next_diagram_row
diagram_sect_start.append(this_jstart)

jj = this_jstart

# accumulate the physical NO processes new 
#  ... initial concentration
section_labels.append("  NO initial")
copyhours_phy(species_process_masses, iNO, jInitial, a_row)
hourly_total_new_no = [e for e in a_row ]
hourly_diagram_values.append([e for e in a_row ])

t_value = total_species_process_masses[iNO][jInitial]
daily_total_new_no = t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... hourly carryover concentration
section_labels.append("  NO carryover")
copyhours_phy(species_process_masses, iNO, jCarryOver, a_row)
accumulate_row(a_row,hourly_total_new_no)
hourly_diagram_values.append([e for e in a_row ])

t_value = total_species_process_masses[iNO][jCarryOver]
daily_total_new_no += 0.0
daily_diagram_values.append(t_value/num_hrs)

jj += 1

#  ... NO loss via O3_o+NO
section_labels.append("  NO + O3o titration") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO2hvO3prod, iNO_o ), a_row)
accumulate_row(a_row,hourly_total_new_no)
hourly_diagram_values.append([e for e in a_row])

daily_diagram_values.append(sum(a_row))

jj += 1

#  ... emissions
section_labels.append("  NO emissions")
copyhours_phy(species_process_masses, iNO, jEmission, a_row)
accumulate_row(a_row,hourly_total_new_no)
hourly_diagram_values.append([e for e in a_row ])

t_value = total_species_process_masses[iNO][jEmission]
daily_total_new_no += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... horiz trans as source
section_labels.append("  NO horiz trans")
copyhours_phy(species_process_masses, iNO, jHor_Trans, a_row)
accumulate_row(a_row,hourly_total_new_no)
hourly_diagram_values.append([e for e in a_row ])

t_value = sum(a_row)
daily_total_new_no += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... vert trans
section_labels.append("  NO vert  trans")
copyhours_phy(species_process_masses, iNO, jVer_Trans, a_row)
accumulate_row(a_row,hourly_total_new_no)
hourly_diagram_values.append([e for e in a_row ])

t_value = sum(a_row)
daily_total_new_no += abs(t_value)
daily_diagram_values.append(t_value)

jj += 1

#  ... entrain/dil
section_labels.append("  NO entrain/dil")
copyhours_phy(species_process_masses, iNO, jEntrain, a_row)
accumulate_row(a_row,hourly_total_new_no)
hourly_diagram_values.append([e for e in a_row ])

t_value = sum(a_row)
daily_total_new_no += t_value
daily_diagram_values.append(t_value)

jj += 1

#   ... secondary new NO from new NO2 photolysis...
section_labels.append("  sec new NO")
hourly_diagram_values.append([e for e in hourly_new_no_frm_no2])
accumulate_row(hourly_new_no_frm_no2,hourly_total_new_no)

daily_total_new_no += daily_no_frm_no2
daily_diagram_values.append(daily_no_frm_no2)

jj += 1

#  ...sum
section_labels.append("Total new NO")
hourly_diagram_values.append([e for e in hourly_total_new_no])

daily_diagram_values.append(daily_total_new_no)

jj += 1

# the NO propagation factor is the fraction of total NO oxidized that
#   gets recreated after NO2 photolysis.  
#   That is, ( recreated NO / tot NO oxidized )
section_labels.append("NO P_n")
hourly_no_pr = [ (r/abs(t)) for (r,t) \
                      in zip(hourly_recreated_no, hourly_tot_no_oxid)] 

hourly_diagram_values.append([e for e in hourly_no_pr])

daily_diagram_values.append(daily_recreated_no/abs(daily_tot_no_oxid)) 

jj += 1

# the NO chain length is 
#   number of times new NO is recreated as NO:
#    if recreated NO = 0 and new NO is 20, chain length = 0.0
#    if recreated NO = new NO, chain length is 1.0
#    if recreated NO = 2 x new NO, chain lenghth is 2.0.
#
section_labels.append("NO chain length")
hourly_no_chain = [ r / q for (r,q) \
                      in zip(hourly_recreated_no, hourly_total_new_no)]

daily_no_chain  = daily_recreated_no / daily_total_new_no

hourly_diagram_values.append([e for e in hourly_no_chain])

daily_diagram_values.append(daily_no_chain)

jj += 1

# how many elements in vector are needed.
num_rows_in_section = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_diagram_row  = jj


#-------------------------------------------------------
### The Diagram Section for ** new O3, O3 srcs/losses **
this_diag_section = next_diag_section  
next_diag_section += 1

#    ... save the name and starting location of this net_rxn
diagram_sect_names.append("New O3, O3 Sources/Losses")
sec_id_num += 1

this_jstart = jstart_next_diagram_row
diagram_sect_start.append(this_jstart)

jj = this_jstart

# accumulate the physical and chemical ozone processes 
#  ... initial concentration
section_labels.append("  O3 initial")
copyhours_phy(species_process_masses, iO3, jInitial, a_row)
hourly_total_o3 = [e for e in a_row ]
hourly_diagram_values.append([e for e in a_row ])

t_value = total_species_process_masses[iO3][jInitial]
daily_total_o3 = t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... carryover concentration
section_labels.append("  O3 carryover")
copyhours_phy(species_process_masses, iO3, jCarryOver, a_row)
accumulate_row(a_row,hourly_total_o3)
hourly_diagram_values.append([e for e in a_row ])

t_value = total_species_process_masses[iO3][jCarryOver]
daily_total_o3 += 0.0
daily_diagram_values.append(t_value)

jj += 1

#  ... emissions
section_labels.append("  O3 emissions")
copyhours_phy(species_process_masses, iO3, jEmission, a_row)
accumulate_row(a_row,hourly_total_o3)
hourly_diagram_values.append([e for e in a_row ])

t_value = total_species_process_masses[iO3][jEmission]
daily_total_o3 += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... horiz trans as source
section_labels.append("  O3 horiz trans")
copyhours_phy(species_process_masses, iO3, jHor_Trans, a_row)
accumulate_row(a_row,hourly_total_o3)
hourly_diagram_values.append([e for e in a_row ])

t_value = sum(a_row)
daily_total_o3 += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... vert trans
section_labels.append("  O3 vert  trans")
copyhours_phy(species_process_masses, iO3, jVer_Trans, a_row)
accumulate_row(a_row,hourly_total_o3)
hourly_diagram_values.append([e for e in a_row ])

t_value = sum(a_row)
daily_total_o3 += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... entrain/dil
section_labels.append("  O3 entrain/dil")
copyhours_phy(species_process_masses, iO3, jEntrain, a_row)
accumulate_row(a_row,hourly_total_o3)
hourly_diagram_values.append([e for e in a_row ])

t_value = sum(a_row)
daily_total_o3 += t_value
daily_diagram_values.append(t_value)

jj += 1

#   ...  new o3 from NO2 photolysis...
section_labels.append("  O3 from NO2+hv")
hourly_diagram_values.append([e for e in hourly_new_o3_prod])
accumulate_row(hourly_new_o3_prod,hourly_total_o3)

daily_total_o3 += daily_new_o3_prod
daily_diagram_values.append(daily_new_o3_prod)

jj += 1

#   ...  o3 organic losses
section_labels.append("  O3 + org loss")
copyhours_net(hourly_net_rxn_masses, i2j(n_OxOrgrad, iO3 ), a_row)
accumulate_row(a_row,hourly_total_o3)
hourly_diagram_values.append([e for e in a_row])

t_value = sum(a_row)
daily_total_o3 += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... O3 loss via O3_o+NO
section_labels.append("  O3o + NO titration") 
copyhours_net(hourly_net_rxn_masses, i2j(n_NO2hvO3prod, iO3_o ), a_row)
accumulate_row(a_row,hourly_total_o3)
hourly_diagram_values.append([e for e in a_row])

t_value = sum(a_row)
daily_total_o3 += t_value
daily_diagram_values.append(t_value)

jj += 1

#   ...  o3 inorganic losses
section_labels.append("  O3 + OH loss")
copyhours_net(hourly_net_rxn_masses, i2j(n_OHOrgOxid, iO3 ), a_row)
accumulate_row(a_row,hourly_total_o3)
hourly_diagram_values.append([e for e in a_row])

t_value = sum(a_row)
daily_total_o3 += t_value
daily_diagram_values.append(t_value)

jj += 1

#   ...  o3 inorganic losses
section_labels.append("  O3 + HO2 loss")
copyhours_net(hourly_net_rxn_masses, i2j(n_HO2toOHrad, iO3 ), a_row)
accumulate_row(a_row,hourly_total_o3)
hourly_diagram_values.append([e for e in a_row])

t_value = sum(a_row)
daily_total_o3 += t_value
daily_diagram_values.append(t_value)

jj += 1

#   ...  o3 photolysis losses
section_labels.append("  O3 + hv loss")
copyhours_net(hourly_net_rxn_masses, i2j(n_O3hvrad, iO1D ), a_row)
accumulate_row(a_row,hourly_total_o3)
hourly_diagram_values.append([e for e in a_row])

t_value = sum(a_row)
daily_total_o3 += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... depositon
section_labels.append("  O3 deposition")
copyhours_phy(species_process_masses, iO3, jDepo, a_row)
accumulate_row(a_row,hourly_total_o3)
hourly_diagram_values.append([e for e in a_row ])

t_value = total_species_process_masses[iO3][jDepo]
daily_total_o3 += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... Temperature adjust
section_labels.append("  O3 Temp adjust")
copyhours_phy(species_process_masses, iO3, jTempAdj, a_row)
accumulate_row(a_row,hourly_total_o3)
hourly_diagram_values.append([e for e in a_row ])

t_value = total_species_process_masses[iO3][jTempAdj]
daily_total_o3 += t_value
daily_diagram_values.append(t_value)

jj += 1

#  ... Final Conc, temp adjusted
section_labels.append("  O3, final(Tadj)")
copyhours_phy(species_process_masses, iO3, jFinalTadj, a_row)
hourly_diagram_values.append([e for e in a_row ])

t_value = total_species_process_masses[iO3][jFinalTadj]
daily_total_o3 += 0.0
daily_diagram_values.append(t_value)

jj += 1

#  ...sum of processes
section_labels.append("Sum of Process O3")
hourly_diagram_values.append([e for e in hourly_total_o3])

daily_diagram_values.append(daily_total_o3)

jj += 1

# how many elements in vector are needed.
num_rows_in_section = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_diagram_row  = jj




## >>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
##  P A R T   F O U R  : P R I N T  T H E  S U M M A R I E S   
##                       P R I N T  T H E  P H Y  P R O C E S S E S
##                       P R I N T  T H E  V O C  C O M P O S I T I O N S
##                       P R I N T  T H E  N E T  R E A C T I O N S
## ===================================================================

# Output the diagram sections summaries to file...
# open the output file...
sum_output_filename = root_output_filename+".sum"
fout = open(sum_output_filename, 'w')

print >>fout, SCRIPT_ID_STRING
print >>fout, __version__
print >>fout, "Diagram Section Summaries"
print >>fout
print >>fout, "IRR file doc line was"
print >>fout, irrfile_docline
print >>fout

kk = 0
print >>fout, "%-21s" % "Ending Hour",
for t in hour_number:
	print >>fout, "        %02d" % t,
print >>fout,  "    Daily"	
print >>fout
print >>fout, "Summary Section Name"
print >>fout, "Item             ppb"
for i in range(len(hourly_diagram_values)):
	if i == diagram_sect_start[kk] :
		if ((kk+1) % 3) == 0:
			print >>fout
			print >>fout, "%-21s" % "Ending Hour",
			for t in hour_number:
				print >>fout, "        %02d" % t,
			print >>fout,  "    Daily"	
		print >>fout
		print >>fout
		print >>fout, "%-20s" % diagram_sect_names[kk]
		print >>fout
		if kk < len(diagram_sect_start)-1:
			kk += 1
	print  >>fout, "%-20s " % section_labels[i],
	for t in range(num_hrs):
		print >>fout, "%10.4f" % (hourly_diagram_values[i][t]),
	print >>fout, "%10.4f" % daily_diagram_values[i]
print >>fout
fout.close()

# Output the physical processes to file...
# open the output file...
phy_output_filename = root_output_filename+".phy"
fout = open(phy_output_filename, 'w')

print >>fout, SCRIPT_ID_STRING
print >>fout, __version__
print >>fout, "Physical Processes By Species"
print >>fout
print >>fout, "IRR file doc line was"
print >>fout, irrfile_docline
print >>fout

print >>fout, "%-21s" % "Ending Hour",
for t in hour_number:
	print >>fout, "        %02d" % t,
print >>fout,  "    Daily"	
print >>fout
print >>fout, "Species"
print >>fout, "  Process, ppb or ppb/h"
for s in range(tot_num_phy_proc_spc):
	if ((s+1) % 8) == 0:
		print >>fout, "%-21s" % "Ending Hour",
		for t in hour_number:
			print >>fout, "        %02d" % t,
		print >>fout,  "    Daily"	
		print >>fout
	# if current species is in the list of spc that need split xport, then...
	if s in splt_xport_spc:
		spc_mark = Phy_spc_names[s] + "--XT"
		print >>fout, "%-20s" % spc_mark
		hour_sum = [0.0]*(num_hrs+1)
		for p in range(num_j_processes):
			if p in split_gainloss_processes:
				for sgn in [1,-1]:
					if sgn == 1:
						proc_name = Proc_Names[p]+", gain"
					else:
						proc_name = Proc_Names[p]+", loss"
					print  >>fout, "  %-18s " % proc_name,
					a_row = [0.0]*num_hrs
					for t in range(num_hrs):
						a_row[t] = species_process_masses[s][t][p]
						if ((sgn > 0) & (a_row[t] < 0.0)):
							a_row[t] = 0.0
						elif ((sgn < 0) & (a_row[t] > 0.0)):
							a_row[t] = 0.0
					for t in range(num_hrs):
						print >>fout, "%10.4f" % (a_row[t]),
						if p in hourly_sum_processes:
							hour_sum[t] += a_row[t]
					total_mass = sum(a_row)
					print >>fout, "%10.4f" % total_mass
					if p in daily_sum_processes:
						hour_sum[-1] += total_mass
			else:
				print  >>fout, "  %-18s " % Proc_Names[p],
				for t in range(num_hrs):
					print >>fout, "%10.4f" % (species_process_masses[s][t][p]),
					if p in hourly_sum_processes:
						hour_sum[t] += species_process_masses[s][t][p]
				print >>fout, "%10.4f" % total_species_process_masses[s][p]
				if p in daily_sum_processes:
					hour_sum[-1] += total_species_process_masses[s][p]
		print  >>fout, "  %-18s " % "Sum Processes",
		for t in range(num_hrs):
			print >>fout, "%10.4f" % (hour_sum[t]),
		print >>fout, "%10.4f" % hour_sum[-1]
		print >>fout
	# also print split xport species the normal way
	print >>fout, "%-20s" % Phy_spc_names[s]
	hour_sum = [0.0]*(num_hrs+1)
	for p in range(num_j_processes):
		print  >>fout, "  %-18s " % Proc_Names[p],
		for t in range(num_hrs):
			print >>fout, "%10.4f" % (species_process_masses[s][t][p]),
			if p in hourly_sum_processes:
				hour_sum[t] += species_process_masses[s][t][p]
		print >>fout, "%10.4f" % total_species_process_masses[s][p]
		if p in daily_sum_processes:
			hour_sum[-1] += total_species_process_masses[s][p]
	print  >>fout, "  %-18s " % "Sum Processes",
	for t in range(num_hrs):
		print >>fout, "%10.4f" % (hour_sum[t]),
	print >>fout, "%10.4f" % hour_sum[-1]
	print >>fout
print >>fout
fout.close()


# Output the VOC composition by processes to file...
# open the output file...
voc_comp_output_filename = root_output_filename+".voc"
fout = open(voc_comp_output_filename, 'w')

print >>fout, SCRIPT_ID_STRING
print >>fout, __version__
print >>fout, "VOC Composition By Processes"
print >>fout
print >>fout, "IRR file doc line was"
print >>fout, irrfile_docline
print >>fout

print >>fout, "%-21s" % "Ending Hour",
for t in hour_number:
	print >>fout, "        %02d" % t,
print >>fout,  "    Daily"	
print >>fout
print >>fout, "Process, ppb or ppb/h"
print >>fout, "  Species"
print >>fout
for p in range(num_j_processes):
	print  >>fout, "%-21s " % Proc_Names[p]
	if p in [jInitial, jEmission, jChemistry, jHor_Trans, jVer_Trans, jEntrain, jFinal]:
		hour_sum = [0.0]*(num_hrs+1)
		if p in [jInitial, jFinal]:
			print >>fout, " %-20s " % "Concs"
		else:
			print >>fout, " %-20s " % "Gains"
		for s in voc_comp_spc:
			print >>fout, "  %-19s" % Phy_spc_names[s],
			daily_total = 0.0
			for t in range(num_hrs):
				c = species_process_masses[s][t][p]
				if c >= 0.0 :
					print >>fout, "%10.4f" % (c),
					hour_sum[t] += c
					daily_total += c
				else:
					print >>fout, "%10.4f" % (0.0),
			if p in [jInitial, jFinal]:
				print >>fout, "%10.4f" % total_species_process_masses[s][p]
				if p == jInitial :
					hour_sum[-1] = hour_sum[0]
				else:
					hour_sum[-1] = hour_sum[-2]
			else:
				print >>fout, "%10.4f" % daily_total
				hour_sum[-1] += daily_total
		print  >>fout, "  %-19s" % "VOCm",
		for t in range(num_hrs):
			print >>fout, "%10.4f" % (hour_sum[t]),
		print >>fout, "%10.4f" % hour_sum[-1]
		print >>fout
	if p in [ jChemistry, jHor_Trans, jVer_Trans, jEntrain, jDepo ]:
		print >>fout, " %-20s " % "Losses"
		hour_sum = [0.0]*(num_hrs+1)
		for s in voc_comp_spc:
			print >>fout, "  %-19s" % Phy_spc_names[s],
			daily_total = 0.0
			for t in range(num_hrs):
				c = species_process_masses[s][t][p]
				if c < 0.0 :
					print >>fout, "%10.4f" % (c),
					hour_sum[t] += c
					daily_total += c
				else:
					print >>fout, "%10.4f" % (0.0),
			print >>fout, "%10.4f" % daily_total
			hour_sum[-1] += daily_total
		print  >>fout, "  %-19s" % "VOCm",
		for t in range(num_hrs):
			print >>fout, "%10.4f" % (hour_sum[t]),
		print >>fout, "%10.4f" % hour_sum[-1]
		print >>fout
print >>fout
print >>fout
print >>fout, "%-21s" % "Ending Hour",
for t in hour_number:
	print >>fout, "        %02d" % t,
print >>fout,  "    Daily"	
print >>fout
print >>fout, "Process, ppbC or ppbC/h"
print >>fout, "  Species"
print >>fout
for p in range(num_j_processes):
	print  >>fout, "%-21s " % Proc_Names[p]
	if p in [jInitial, jEmission, jChemistry, jHor_Trans, jVer_Trans, jEntrain, jFinal]:
		hour_sum = [0.0]*(num_hrs+1)
		if p in [jInitial, jFinal]:
			print >>fout, " %-20s " % "Concs"
		else:
			print >>fout, " %-20s " % "Gains"
		for s in voc_comp_spc:
			print >>fout, "  %-19s" % Phy_spc_names[s],
			daily_total = 0.0
			for t in range(num_hrs):
				c = VOC_carb_num[s]*species_process_masses[s][t][p]
				if c >= 0.0 :
					print >>fout, "%10.4f" % (c),
					hour_sum[t] += c
					daily_total += c
				else:
					print >>fout, "%10.4f" % (0.0),
			if p in [jInitial, jFinal]:
				c = VOC_carb_num[s]*total_species_process_masses[s][p]
				print >>fout, "%10.4f" % (c)
				if p == jInitial :
					hour_sum[-1] = hour_sum[0]
				else:
					hour_sum[-1] = hour_sum[-2]
			else:
				print >>fout, "%10.4f" % daily_total
				hour_sum[-1] += daily_total
		print  >>fout, "  %-19s" % "VOCm",
		for t in range(num_hrs):
			print >>fout, "%10.4f" % (hour_sum[t]),
		print >>fout, "%10.4f" % hour_sum[-1]
		print >>fout
	if p in [ jChemistry, jHor_Trans, jVer_Trans, jEntrain, jDepo ]:
		print >>fout, " %-20s " % "Losses"
		hour_sum = [0.0]*(num_hrs+1)
		for s in voc_comp_spc:
			print >>fout, "  %-19s" % Phy_spc_names[s],
			daily_total = 0.0
			for t in range(num_hrs):
				c = VOC_carb_num[s]*species_process_masses[s][t][p]
				if c < 0.0 :
					print >>fout, "%10.4f" % (c),
					hour_sum[t] += c
					daily_total += c
				else:
					print >>fout, "%10.4f" % (0.0),
			print >>fout, "%10.4f" % daily_total
			hour_sum[-1] += daily_total
		print  >>fout, "  %-19s" % "VOCm",
		for t in range(num_hrs):
			print >>fout, "%10.4f" % (hour_sum[t]),
		print >>fout, "%10.4f" % hour_sum[-1]
		print >>fout
print >>fout
fout.close()


# Output the hourly and total net reactions to file...
# open the output file...
net_output_filename = root_output_filename+".net"
fout = open(net_output_filename, 'w')

print >>fout, SCRIPT_ID_STRING
print >>fout, __version__
print >>fout, "Net Reactions"
print >>fout
print >>fout, "IRR file doc line was"
print >>fout, irrfile_docline
print >>fout

kk = 0
print >>fout, "%-21s" % "Ending Hour",
for t in hour_number:
	print >>fout, "        %02d" % t,
print >>fout,  "     Daily"	
print >>fout
print >>fout, "Net Reaction Name"
print >>fout, "Species     ppb"
for i in range(len(net_rxn_masses)):
	if i == net_rxn_jindex[kk] :
		if ((kk+1) % 4) == 0:
			print >>fout
			print >>fout, "%-21s" % "Ending Hour",
			for t in hour_number:
				print >>fout, "        %02d" % t,
			print >>fout,  "    Daily"	
		print >>fout
		print >>fout
		print >>fout, "%-20s" % net_rxn_names[kk]
		print >>fout
		if kk < len(net_rxn_names)-1:
			kk += 1
	print  >>fout, "%-20s " % Net_spc_names[net_rxn_spcname[i]],
	for t in range(0,len(hour_number)):
		print >>fout, "%10.4f" % (hourly_net_rxn_masses[t][i]),
	print >>fout, "%10.4f" % daily_net_rxn_masses[i]
# finished the whole set of net rxn masses for all hours
print >>fout
fout.close()


print "F I N I S H E D"
