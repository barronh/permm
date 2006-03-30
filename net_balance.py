#!/usr/bin/env python2.3
#

import os.path
import sys

import re
import copy   # for deepcopy
	
	

### beginning of MAIN BODY #######
"""
  From the command line this script takes a file name and two integers for begin and end hour of interest
  It currently uses a hardcoded CB4 mechanism, nets the reactions and returns them to the screen.
"""
from optparse import OptionParser

# create a cmd line parser...
parser = OptionParser("usage: %prog filename.ext ")
parser.add_option("-s", "--show",
	dest="show",
	default=False,
	action="store_true",
	help="Display graphs on screen")


# get options and arguments
(options, args) = parser.parse_args()

# requires filename as argument
if len(args) != 1:
	parser.error("Invalid number of arguments")
# requires a valid ext file
elif not os.path.exists(args[0]):
	parser.error("Input file does not exist")


	
input_filename = args[0]

time_re  = re.compile('Time =[0-9]{6}', re.IGNORECASE)
irr_re   = re.compile('\{\s*\d+\}\s+\d+', re.IGNORECASE)
ipr_re   = re.compile('"\w+\s*"\s*', re.IGNORECASE)
split_re = re.compile('[ ]+')


f = open(input_filename, 'r')
# first two lines are 'doc' lines giving data source
doc1 = f.readline()
line = f.readline()

print "IRR file doc line was"
print doc1
print


def get_irr_data(f):
	"""
	f -- opened file
	return (time (as integer), ir_rates (as list), ip_rates (as list of lists)
	if time < 0, end of file
	"""
	# initialize ir and ip storage...
	ir_time  = -1
	ir_rates = [0,]
	ip_rates = []
	for i in range(0,23):  # make space for process rates for CB4's 24 species
		ip_rates.append([])
	
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
	
	# for now skip reading the process rates..
	# read in the ! Species      Initial conc. ... header
	line = f.readline()
	while line[0] != ';' :
		line = f.readline()
	
	return (ir_time, ir_rates, ip_rates)





#  GLOBAL DATA FOR CAMX CB4 Mechanism 3 ...

SPC_Names = [
  'NO  ', 'NO2 ', 'O3  ', 'OLE ', 'PAN ', 'N2O5', 'PAR ', 'TOL ', 'XYL ',\
  'FORM', 'ALD2', 'ETH ', 'CRES', 'MGLY', 'OPEN', 'PNA ', 'CO  ', 'HONO',\
  'H2O2', 'HNO3', 'ISOP', 'MEOH', 'ETOH', 'CH4 ', 'O   ', 'OH  ', 'HO2 ', 'NO3 ',\
  'C2O3', 'XO2 ', 'XO2N', 'NTR ', 'CRO ', 'ISPD', 'TO2 ', 'ROR ', 'SO2 ' ]

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

max_i_spc = iSO2 + 1

# Start with info about the net_rxns sets...
#   allocate vector of names of the net_rxns
net_rxn_names = []
#   this is indexed by next_net_rxn_set

# allocate a look up table between global species ID num (called iSPC)
#    and net_rxn ID num (called jSPC) [where SPC is any species name in 
#    the above global list].
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
num_net_rxn_sets = next_net_rxn_set



# allocate net rxn masses vector
net_rxn_masses = []
#   adjacent elements in this numeric vector hold and accumulate
#    the net masses for species iSPC in the set of rxns being lumped

total_net_rxn_masses = []
#   accumulate the net_rxn_masses over all time.

# allocate net rxn species names index vector
net_rxn_spcname = []
#   adjacent elements in this integer vector hold the iSPC value
#    for the net masses in net_rxn_masses

# next location in net_rxn_masses vector to append another net rxn set
#   and in the net_rxn_spcname vector to insert the iSPC number
jstart_next_nr_set = 0

# create a function to return the j-index associated with species i in net 
#    reaction set k
def i2j(k,i):
	j = net_rxn_species[k][i]
	# print "i2j  k i  j", (k, i, j)
	return j
	

######### repeat section below for each set of net reactions ######

# The Net Rxns for ** OH + organic **
this_net_rxn_set = next_net_rxn_set  
next_net_rxn_set += 1


#    ... save the name and starting location of this net_rxn
net_rxn_names.append('OH+organic')

this_jstart = jstart_next_nr_set

# save the starting location of species for this net_reaction
net_rxn_jindex.append(this_jstart)


# create a cross reference vector for iSPC to j and fill it with jNONE
jndx_net_rxn = [jNONE]*max_i_spc
indx_net_rxn = []

#  fill in only the species numbers for species in 
#    the net_rxn for the 'OH+organic' reactions

# fill in the jindex at the iSPC locations for the 'OH+organic' net rxn
jj = this_jstart

jndx_net_rxn[iNO2 ] = jj; jj += 1; indx_net_rxn.append(iNO2 )
jndx_net_rxn[iOLE ] = jj; jj += 1; indx_net_rxn.append(iOLE )
jndx_net_rxn[iPAR ] = jj; jj += 1; indx_net_rxn.append(iPAR )
jndx_net_rxn[iTOL ] = jj; jj += 1; indx_net_rxn.append(iTOL )
jndx_net_rxn[iXYL ] = jj; jj += 1; indx_net_rxn.append(iXYL )
jndx_net_rxn[iFORM] = jj; jj += 1; indx_net_rxn.append(iFORM)
jndx_net_rxn[iALD2] = jj; jj += 1; indx_net_rxn.append(iALD2)
jndx_net_rxn[iETH ] = jj; jj += 1; indx_net_rxn.append(iETH )
jndx_net_rxn[iCRES] = jj; jj += 1; indx_net_rxn.append(iCRES)
jndx_net_rxn[iMGLY] = jj; jj += 1; indx_net_rxn.append(iMGLY)
jndx_net_rxn[iOPEN] = jj; jj += 1; indx_net_rxn.append(iOPEN)
jndx_net_rxn[iCO  ] = jj; jj += 1; indx_net_rxn.append(iCO  )
jndx_net_rxn[iHNO3] = jj; jj += 1; indx_net_rxn.append(iHNO3)
jndx_net_rxn[iISOP] = jj; jj += 1; indx_net_rxn.append(iISOP)
jndx_net_rxn[iMEOH] = jj; jj += 1; indx_net_rxn.append(iMEOH)
jndx_net_rxn[iETOH] = jj; jj += 1; indx_net_rxn.append(iETOH)
jndx_net_rxn[iCH4 ] = jj; jj += 1; indx_net_rxn.append(iCH4 )
jndx_net_rxn[iOH  ] = jj; jj += 1; indx_net_rxn.append(iOH  )
jndx_net_rxn[iHO2 ] = jj; jj += 1; indx_net_rxn.append(iHO2 )
jndx_net_rxn[iC2O3] = jj; jj += 1; indx_net_rxn.append(iC2O3)
jndx_net_rxn[iXO2 ] = jj; jj += 1; indx_net_rxn.append(iXO2 )
jndx_net_rxn[iXO2N] = jj; jj += 1; indx_net_rxn.append(iXO2N)

# how many elements in vector are needed.
num_spc_in_net_rxn = jj - this_jstart

# save the starting value of the vector index offset for the next cycle...
jstart_next_nr_set  = jj

# DEEP copy a new jndx vector for 'OH+organic' net_rxn to the net_rxn_species array
net_rxn_species.append(copy.deepcopy(jndx_net_rxn))

# add and init new elements to the net_rxn_masses vector for species in the 
#   OH+organic net rxn
for z in [0]*num_spc_in_net_rxn:
	net_rxn_masses.append(z)
	total_net_rxn_masses.append(z)
	
# and add the names of the species to the net_rxn_spcname vector
for n in indx_net_rxn:
	net_rxn_spcname.append(n)
	
# kk is this net_rxn set number
kk = this_net_rxn_set   

	
num_net_rxn_sets  = this_net_rxn_set
max_j  = jstart_next_nr_set


(time, ir, ip) = get_irr_data(f);  ## read new data and...
while time >= 0:
	print "-"*20
	print "Time ", time
	print
	print "For Net Reaction ", net_rxn_names[kk]
	
	# ... process the ir data...
	#      ... zero out the net reaction masses
	for i in range(0,len(net_rxn_masses)):
		net_rxn_masses[i] = 0.0
		
	# compute the losses and gains in this net_rxn set...
	
	# { 26} NO2+OH=HNO3
	# { 36} OH+CO=HO2
	# { 37} FORM+OH=HO2+CO
	# { 43} ALD2+OH=C2O3
	# { 51} OH+CH4=FORM+XO2+HO2
	# { 52} PAR+OH=0.87*XO2+0.13*XO2N+0.11*HO2+0.11*ALD2+-0.11*PAR+0.76*ROR
	# 	  { 53} ROR=0.96*XO2+1.1*ALD2+0.94*HO2+-2.1*PAR+0.04*XO2N
	# 	  { 54} ROR=HO2
	# 	  { 55} ROR+NO2=NTR    ! OMITTED from this net reaction set
	# { 57} OH+OLE=FORM+ALD2+-1*PAR+XO2+HO2
	# { 61} OH+ETH=XO2+1.56*FORM+0.22*ALD2+HO2
	# { 63} TOL+OH=0.44*HO2+0.08*XO2+0.36*CRES+0.56*TO2
	# 	  { 64} TO2+NO=0.9*NO2+0.9*HO2+0.9*OPEN+0.1*NTR
	# 	  { 65} TO2=CRES+HO2
	# { 66} OH+CRES=0.4*CRO+0.6*XO2+0.6*HO2+0.3*OPEN
	# 	  { 68} CRO+NO2=NTR
	# { 70} OPEN+OH=XO2+2*CO+2*HO2+C2O3+FORM
	# { 72} OH+XYL=0.7*HO2+0.5*XO2+0.2*CRES+0.8*MGLY+1.1*PAR+0.3*TO2
	# { 73} OH+MGLY=XO2+C2O3
	# { 76} OH+ISOP=0.912*ISPD+0.629*FORM+0.991*XO2+0.912*HO2+0.088*XO2N
	# { 84} MEOH+OH=FORM+HO2
	# { 85} ETOH+OH=HO2+ALD2
	# { 92} OH+ISPD=1.565*PAR+0.167*FORM+0.713*XO2+0.503*HO2+0.334*CO+0.168*MGLY+
	# 			  0.273*ALD2+0.498*C2O3
	
	# notes:: 
	#   1) for ir[55], ROR, NTR is from NO2+rad not NO+RO2, omitted here
	#   2) for ir[64], TO2, the 0.9*NO2 is coded as 0.9*XO2
	#                       and the NO is omitted as a reactant
	#                       and the 0.1*NTR is coded as 0.1*XO2N
	#          ir[65], TO2, included as a HO2 source
	#   3) for ir[66], CRO, we have omitted the CRO production.
	#                   the ir[68] production of NTR is from NO2, not NO
	#                   thus it is omitted in this reaction..

	
	# reactant losses in OH + org...
	#   ... the organics losses
	net_rxn_masses[i2j(kk,iNO2 )] = -ir[26]
	net_rxn_masses[i2j(kk,iOLE )] = -ir[57]
	net_rxn_masses[i2j(kk,iPAR )] = -ir[52]
	net_rxn_masses[i2j(kk,iTOL )] = -ir[63]
	net_rxn_masses[i2j(kk,iXYL )] = -ir[72]
	net_rxn_masses[i2j(kk,iFORM)] = -ir[37]
	net_rxn_masses[i2j(kk,iALD2)] = -ir[43]
	net_rxn_masses[i2j(kk,iETH )] = -ir[61]
	net_rxn_masses[i2j(kk,iCRES)] = -ir[66]
	net_rxn_masses[i2j(kk,iMGLY)] = -ir[73]
	net_rxn_masses[i2j(kk,iOPEN)] = -ir[70]
	net_rxn_masses[i2j(kk,iCO  )] = -ir[36]
	net_rxn_masses[i2j(kk,iISOP)] = -ir[76]
	net_rxn_masses[i2j(kk,iMEOH)] = -ir[84]
	net_rxn_masses[i2j(kk,iETOH)] = -ir[85]
	net_rxn_masses[i2j(kk,iCH4 )] = -ir[51]
	net_rxn_masses[i2j(kk,iISPD)] = -ir[92]
	#    ... the OH losses
	net_rxn_masses[i2j(kk,iOH  )] = -ir[26]-ir[57]-ir[52]-ir[63]-ir[72]-ir[37]\
									-ir[43]-ir[61]-ir[66]-ir[73]-ir[70]-ir[36]\
									-ir[76]-ir[84]-ir[85]-ir[51]-ir[92]
	#    ... the radical products...
	net_rxn_masses[i2j(kk,iHO2 )] =  ir[36]+ir[37]+ir[51]+0.11*ir[52]\
									+ir[57]+ir[61]+0.44*ir[63]+0.6*ir[66]\
									+2*ir[70]+0.7*ir[72]+0.912*ir[76]+ir[84]\
									+ir[85]+0.503*ir[92]+0.94*ir[53]+ir[54]\
									+0.9*ir[64]+ir[65]
	
	net_rxn_masses[i2j(kk,iC2O3)] =  ir[43]+ir[70]+ir[73]+0.498*ir[92]
	
	net_rxn_masses[i2j(kk,iXO2 )] =  ir[51]+0.87*ir[52]+ir[57]+ir[61]\
									+0.08*ir[63]+0.6*ir[66]\
									+ir[70]+0.5*ir[72]+ir[73]+0.991*ir[76]\
									+0.713*ir[92]+0.96*ir[53]+0.9*ir[64]
									
	#     ... radical losses via NO2 or NO->NO2 reactions ...
	net_rxn_masses[i2j(kk,iHNO3)] =  ir[26]
	net_rxn_masses[i2j(kk,iXO2N)] =  0.13*ir[52]+0.088*ir[76]+0.04*ir[53]\
										+0.1*ir[64]
	
	
	
	(time, ir, ip) = get_irr_data(f);  ## read new data and...
	
	
	for i in range(0,len(net_rxn_masses)):
		print '%s  %10.6f' % (SPC_Names[net_rxn_spcname[i]], net_rxn_masses[i])
		total_net_rxn_masses[i] += net_rxn_masses[i]
	print
	print



######### repeat to here for each set of net reactions ######

print "-"*20
print "-"*20
print "Total of whole period:"
print
print "For Net Reaction ", net_rxn_names[kk]
print

for i in range(0,len(net_rxn_masses)):
	print '%s  %10.6f' % (SPC_Names[net_rxn_spcname[i]], total_net_rxn_masses[i])
print
print

print "F I N I S H E D"
