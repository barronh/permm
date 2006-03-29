#!/usr/bin/env python2.3

import os.path
import sys

from Numeric import array, zeros, reshape, resize, transpose, Float as nFloat, add as nAdd

def parseext(filename):
	"""
	parseext parses ext file into irr rates and ipr quantities
	"""
	import re
	
	irr_rates = {}
	ipr_qtys = {}
	
	irr_re = re.compile('\{\s*\d+\}\s+\d+', re.IGNORECASE)
	ipr_re = re.compile('"\w+\s*"\s*', re.IGNORECASE)
	time_re = re.compile('Time =[0-9]{6}', re.IGNORECASE)
	splitter = re.compile('[ ]+')
	time=''
	f = open(filename, 'r')
	for line in f:
		line = line.replace('\n','')
		if time_re.match(line) != None:
			time = line.split('=')[1]
			irr_rates[time] = []
			ipr_qtys[time] = {}
		if irr_re.match(line) != None:
			irr_rates[time].append(float(splitter.split(line.replace('"','').replace('{','').replace('}','').strip())[1]))
		elif ipr_re.match(line) != None:
			ipr = splitter.split(line.replace('"','').replace('{','').replace('}','').strip())
			ipr_qtys[time][ipr[0]] = ipr[1:]
	return irr_rates, ipr_qtys
	
def NetReactions(species, reactions, rates, included):
	"""
	NetReactions shapes one dimensional the 1 dimensional rates and inclusion arrays
	to match the three dimensional reactions array
	"""
	#extend rates vertically to match species
	included = array([[included]*len(species)]*2)
	rates = array([[rates]*len(species)]*2)

	#transpose to match reactions layout
	rates = transpose(rates,(2,1,0))
	included = transpose(included,(2,1,0))

	net = nAdd.reduce(rates * reactions * included)
	
	return net

if __name__ == '__main__':
	"""
	  From the command line this script takes a file name and two integers for begin and end hour of interest
	  It currently uses a hardcoded CB4 mechanism, nets the reactions and returns them to the screen.
	"""
	from optparse import OptionParser
	
	#initialize output
	output = []
	parser = OptionParser("usage: %prog [-oOUTDIR] filename.ext iHrStart iHrEnd\niHrStart and iHrEnd should be between 1 and 24")
	parser.add_option("-s", "--show",
		dest="show",
		default=False,
		action="store_true",
		help="Display graphs on screen")
 
	parser.add_option("-o", "--output",
		dest="outdir",
		metavar="OUTDIR",
		help="Destination directory for output graphs")
	
	#get options and arguments
	(options, args) = parser.parse_args()
	
	#requires 3 arguments
	if len(args) != 3:
		parser.error("Invalid number of arguments")
	#requires a valid ext file
	elif not os.path.exists(args[0]):
		parser.error("Input file does not exist")
	#requires integer inputs
	elif float(args[1]) != int(args[1]) or float(args[2]) != int(args[2]):
		parser.error("iHrStart and End must be integers")
	#requires hours between 1 and 24
	elif int(args[1]) > 24 or int(args[1]) < 1 or int(args[2]) > 24 or int(args[2]) < 1:
		parser.error("iHrStart and End must be between 1 and 24")
	#the end must be equal to or greater than the beginning
	elif int(args[1]) > int(args[2]):
		parser.error('IHrEnd must be atleast as big as iHrStart')
	
	#if an output directory is supplied, ensure that it exists
	if options.outdir and not os.path.exists(options.outdir):
		#create it if it does not
		os.makedirs(options.outdir)

    #integer indexing for CB4 species
	iNULL = -1
	iNO   = 0 
	iNO2  = 1 
	iO3   = 2 
	iOLE  = 3 
	iPAN  = 4 
	iN2O5 = 5 
	iPAR  = 6 
	iTOL  = 7 
	iXYL  = 8 
	iFORM = 9 
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
	iC2O3 = 27
	iXO2  = 28
	iXO2N = 29
	iNTR  = 30
	iCRO  = 31
	iISPD = 32
	iTO2  = 33
	iROR  = 34
	iSO2  = 35
	
	#creating the look up table between global species ID num. and cycle ID num.
	cycle_ref = []

	next_cycle = 0


	#parse rates from ext file
	parsedrates = parseext(args[0])
	
	##### Some type of loop here over the hours....
	ir=parsedrates[0]['150000']
	
	ir.insert(0,0.00) # account for 1-origin indexing of rxns
	                  # eg, rxn 1..96

	ii = 0
	print "The current IR vector"
	for r in ir :
		print ii, r
		ii+=1
	print

	#allocate net mass vector
	net_mass = []
    
    
	#creating a temp NULL cross reference vector
	temp_xref_vec = [iNULL]*36
    
	next_j = 0
    
 ######### repeat this section for each set of net reactions ######
 
 
 
 # net mass accumulation for OH organic loss

    # number this cycle..
	cycle_OH_loss = next_cycle
	next_cycle += 1


    #  fill in correct cross reference for OH organic reactions


	jOLE  = 0  + next_j 
	jPAR  = 1  + next_j 
	jTOL  = 2  + next_j
	jXYL  = 3  + next_j 
	jFORM = 4  + next_j 
	jALD2 = 5  + next_j 
	jETH  = 6  + next_j 
	jCRES = 7  + next_j 
	jMGLY = 8  + next_j 
	jOPEN = 9  + next_j 
	jCO   = 10 + next_j
	jISOP = 11 + next_j 
	jMEOH = 12 + next_j 
	jETOH = 13 + next_j 
	jCH4  = 14 + next_j
	jOH   = 15 + next_j
	jHO2  = 16 + next_j
	jC2O3 = 17 + next_j
	jXO2  = 18 + next_j
	jXO2N = 19 + next_j
	jNTR  = 20 + next_j
	jCRO  = 21 + next_j
	jISPD = 22 + next_j
	jTO2  = 23 + next_j
	jROR  = 24 + next_j
	
	# save the starting value of the array index offset for the next cycle...
	next_j = jROR + 1
    
    # create the contents of the vector for the OH organic loss cycle
	temp_xref_vec[iOLE ] =  jOLE 
	temp_xref_vec[iPAR ] =  jPAR 
	temp_xref_vec[iTOL ] =  jTOL 
	temp_xref_vec[iXYL ] =  jXYL 
	temp_xref_vec[iFORM] =  jFORM
	temp_xref_vec[iALD2] =  jALD2
	temp_xref_vec[iETH ] =  jETH 
	temp_xref_vec[iCRES] =  jCRES
	temp_xref_vec[iMGLY] =  jMGLY
	temp_xref_vec[iOPEN] =  jOPEN
	temp_xref_vec[iCO  ] =  jCO  
	temp_xref_vec[iISOP] =  jISOP
	temp_xref_vec[iMEOH] =  jMEOH
	temp_xref_vec[iETOH] =  jETOH
	temp_xref_vec[iCH4 ] =  jCH4 
	temp_xref_vec[iOH  ] =  jOH  
	temp_xref_vec[iHO2 ] =  jHO2 
	temp_xref_vec[iC2O3] =  jC2O3
	temp_xref_vec[iXO2 ] =  jXO2 
	temp_xref_vec[iXO2N] =  jXO2N
	temp_xref_vec[iNTR ] =  jNTR
	temp_xref_vec[iCRO ] =  jCRO 
	temp_xref_vec[iISPD] =  jISPD
	temp_xref_vec[iTO2 ] =  jTO2 
	temp_xref_vec[iROR ] =  jROR 

	# add the xref vector for OH Cycle to the cycle_ref array
	cycle_ref.append(temp_xref_vec)
    
	print temp_xref_vec
	print
    
    #Build the net mass for the OH loss reactions

	for z in [0]*25:
		net_mass.append(z)


	print "The contents of net_mass.."
	ii = 0
	for n in net_mass:
		print ii, n
		ii+=1
	print
	

	net_mass[jOLE ] = -ir[57]
	net_mass[jPAR ] = -ir[52]
	net_mass[jTOL ] = -ir[63]
	net_mass[jXYL ] = -ir[72]
	net_mass[jFORM] = -ir[37]
	net_mass[jALD2] = -ir[43]
	net_mass[jETH ] = -ir[61]
	net_mass[jCRES] = -ir[66]
	net_mass[jMGLY] = -ir[73]
	net_mass[jOPEN] = -ir[70]
	net_mass[jCO  ] = -ir[36]
	net_mass[jISOP] = -ir[76]
	net_mass[jMEOH] = -ir[84]
	net_mass[jETOH] = -ir[85]
	net_mass[jCH4 ] = -ir[51]
	net_mass[jOH  ] = -ir[52]-ir[57]-ir[63]-ir[72]-ir[37]-ir[43]-ir[61]-ir[66]-\
						ir[73]-ir[70]-ir[36]-ir[76]-ir[84]-ir[85]-ir[51]
	net_mass[jHO2 ] =  ir[36]+ir[37]+ir[51]+0.11*ir[52]+ir[57]+ir[61]+0.44*ir[63]+\
						0.6*ir[66]+2*ir[70]+0.7*ir[72]+0.912*ir[76]+ir[84]+\
						ir[85]+.503*ir[92]+0.94*ir[53]+ir[54]+0.9*ir[64]+ir[65]
	net_mass[jC2O3] =  ir[43]+ir[70]+ir[73]+0.498*ir[92]
	net_mass[jXO2 ] =  ir[51]+0.87*ir[52]+ir[57]+ir[61]+0.08*ir[63]+0.6*ir[66]+\
						ir[70]+0.5*ir[72]+ir[73]+0.991*ir[76]+0.713*ir[92]+\
						0.96*ir[53]+0.9*ir[64]
	net_mass[jXO2N] =  0.13*ir[52]+0.088*ir[76]+0.04*ir[53]
	net_mass[jNTR ] =  ir[55]+0.1*ir[64]
	net_mass[jCRO ] =  0.4*ir[66]
	net_mass[jISPD] =  -ir[92]+0.912*ir[76]
	net_mass[jTO2 ] =  0.56*ir[63]+0.30*ir[72]-ir[64]-ir[65]
	net_mass[jROR ] =  0.76*ir[52]

    # notes:: 
    #   1) for ir[56], ROR, we have omitted the NO2 loss
    #   2) for ir[64], TO2, the 0.9*NO2 is coded as 0.9*XO2
    #                       and the NO is omitted as a reactant
    #   3) for ir[68], CRO, we have omitted the NO2 loss
    


  
	for i in range(0,next_j):
		print i, net_mass[i]

 #test2#		
 ######### repeat to here for each set of net reactions ######

