#!/usr/bin/python
import os.path
import sys
import operator
import re

from xml.dom import minidom
from optparse import OptionParser

SCRIPT_ID_STRING = "whoknows.py, 2006-04-09 (c) 2006-04-09 (m)"

__version__ = "R8"


# create a cmd line parser...
usage  = "usage: %prog [options] XMLFILE"
version = "%%prog %s" % __version__
parser = OptionParser(usage=usage, version=version)

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

# extent of output
parser.add_option("-s", "--starthr",
	dest="starthr", 
	default=-1,
	help="starting hour")

parser.add_option("-e", "--endhr",
	dest="endhr", 
	default=-1,
	help="end hour")

	
# get options and arguments
(options, args) = parser.parse_args()

if len(args)==0:
	pass
	#for testing
	#parser.error("Invalid number of arguments")
else:
	filename = args[0]
	
starthr=-1
endhr=-1


def SingleTextNode(parent, tag):
	return SingleNode(parent,tag).childNodes[0].nodeValue

def SingleNode(parent, tag):
	nodes = parent.getElementsByTagName(tag)
	if nodes!=None:
		if len(nodes)>1:
			raise ValueError, "Source has more than one %s node" % tag
		elif len(nodes)==0:
			raise ValueError, "Source has fewer than one %s node" % tag
		else:
			return nodes[0]


def ParseDiagramDictionary(filename):
	doc = minidom.parse(filename)
	
	delim = re.compile(SingleTextNode(doc,'delimiter'))
	starthr = int(SingleTextNode(doc,'starthour'))
	endhr = int(SingleTextNode(doc,'endhour'))
	Summaries = []
	Reactions = []
	Processes = []
	TimeSeries = {}
	#child nodes of summaries are summary
	for summary in SingleNode(doc,'summaries').getElementsByTagName('summary'):
		stitle = SingleTextNode(summary, 'title').strip()
		Summaries.append({'title': stitle})
		#sources are the components of the diagram
		for source in summary.getElementsByTagName('source'):
			#attributes should include name, else raise error
			sname = SingleTextNode(source,'name').strip()
			
			#there should only be one value node per source; values has one text value
			svals = [float(val) for val in delim.split(SingleTextNode(source, 'values'))]
			
			#there should only be one value node per source; values has one text value
			stot = [float(val) for val in delim.split(SingleTextNode(source, 'total'))]
			if len(stot) != 1:
				raise ValueError, "Total should be a single value"
			else:
				if Summaries[-1].has_key(sname):
					raise ValueError, "%s summary has more than one %s instance" % (stitle, sname)
				else:
					Summaries[-1][sname] = {'values': svals, 'total' : stot[0]}
		
	# Each children of netreactions are reactions
	for reaction in SingleNode(doc,'reactions').getElementsByTagName('reaction'):
		rtitle = SingleTextNode(reaction,'title').strip()
		Reactions.append({'title': rtitle})
		#reactions are made up of products and reactants
		for species in reaction.getElementsByTagName('species'):
			#nodes should include name, else raise error
			rname = SingleTextNode(species,'name').strip()

			#there should only be one value node per source; values has one text value
			rvals = [float(val) for val in delim.split(SingleTextNode(species, 'values'))]
			
			#there should only be one value node per source; values has one text value
			rtot = [float(tot) for tot in delim.split(SingleTextNode(species, 'total'))]
			if len(rtot) != 1:
				raise ValueError, "Total should be a single value"
			else:
				#I think this can be refacotred to look like processes
				#having the species as the first key would make the ts
				#process simpler without using the timeseries seq
				#i.e.  if Reactions = {'O3',{'Rxn1':[1,2,3], 'Rxn2':[6,5,4]}, 
				#then Reactions.iteritems() yields 'O3' + RxnDict
				#and RxnDict.iteritems() = TimeSeries
				#checking in before refactoring
				if Reactions[-1].has_key(rname):
					raise ValueError, "%s reaction has more than one %s instance" % (rtitle, rname)
				else:
					Reactions[-1][rname] = {'values': rvals, 'total' : rtot[0]}
				if not TimeSeries.has_key(rname):
					TimeSeries[rname] =[]
				
				TimeSeries[rname].append(rvals)
	
	# Each child of processes is a process
	for process in []: #SingleNode(doc,'processes').getElementsByTagName('process'):
		ptitle = SingleTextNode(process,'title').strip()
		#processes affect species
		for species in process.getElementsByTagName('species'):
			#species should include a name node
			pname = SingleTextNode(species,'name').strip()
			
			#there should only be one value node per species
			pvals = [float(val) for val in delim.split(SigleTextNode(species, 'values'))]

			#there should only be one value node per species
			ptot = [float(tot) for tot in delim.split(SigleTextNode(species, 'total'))]
			
			if len(ptot) != 1:
				raise ValueError, "Total should be a single value"
			else:
				if not Processes.has_key(pname):
					Processes[pname]={}
				if Processes[pname].has_key(ptitle):
					raise ValueError, "%s process has more than one %s instance" % (pname, ptitle)
				else:
					Processes[pname][ptitle] = {'values': pvals, 'total' : ptot[0]}

	return {'summaries': Summaries, 'reactions': Reactions, 'processes': Processes, 'starthour': starthr, 'endhour': endhr, 'timeseries': TimeSeries}


#For testing purposes, the file name is assumed.
filename = '/Users/barronh/Documents/Development/net_balance/inout/Bayland_4timeCO_0825_818_420.xml'

#Parse the xml file into the dictionary including 
#reactions, summaries, timeseries, delim, start and end
MasterDict = ParseDiagramDictionary(filename)

#if start and end were undefined; set them to the bounds
if starthr == -1: starthr = MasterDict['starthour']
if endhr == -1: endhr = MasterDict['endhour']


#convert values to range position
istarthr = int(starthr - MasterDict['starthour'])
iendhr = int(endhr - MasterDict['starthour'])+1

def ValueLookup(diagram_list,master):
	"""
	Takes a dictionary of "mappings" for boxes on diagrams 
	(eg reacted NO: summaries, NO summary, reacted NO);
	returns a dictionary of hour vectors for those boxes
	"""
	DiaValueDict = None
	def XMLFindIt(domainstr, grouptitle, itemname):
		"""
		XMLFindIt uses the 3 pieces of info to identify the source 
		or species	of interest in the net reactions output
		
		*assumes the existence of master from parent function
		"""
		try:
			retval = [group[itemname]['values'] for group in master[domainstr] if group['title']==grouptitle]
		except KeyError:
			raise KeyError, "Couldn't find %s in the %s from %s" % (itemname, grouptitle, domainstr)

		if len(retval)!=1: 
			raise ValueError, "%s in the %s from %s had more than one answer" % (itemname, grouptitle, domainstr)
		else:
			return retval[0]

	DiaValueDict = dict([(diagram, [reduce(operator.add, hrs[istarthr:iendhr]) for hrs in [XMLFindIt(domain, group, item) for domain, group, item in mappings]]) for diagram, mappings in diagram_list.iteritems()])

	return DiaValueDict

def CreateTSChemArrays(tsdict):
	"""
	Create time series arrays by taking a dictionary of lists of lists
	each key is a species.  Each value is a list of values from
	different reactions in which the species is a reactant or product.
	The function takes slice on the time axis and sums across the reaction
	axis.	
	"""
	#action!
	return dict([(spc, [reduce(operator.add, s[istarthr:iendhr]) for s in zip(*[s[istarthr:iendhr] for s in seq])]) for spc, seq in tsdict.iteritems()])

def CreateTSProcArrays(procdict):
	#Returns dict by species of dicts by proc
	#The only thing missing for a TS would be chemistry.
	return dict([(spc, dict([(proc, vals[istarthr:iendhr]) for proc, vals in procs.iteritems()])) for spc, procs in procdict])

#DiaMapDict will contain keys to the values from
#summaries and net reactions.  The keys will identify
#the values for use in each diagram
DiaMapDict ={}
DiaMapDict['OH_VOC_DIAGRAM']  = [ \
			('summaries','C2O3 New, Prod, and Loss','Ald+hv  C2O3'), \
			('reactions', 'O3+hv radical source','O3')
			]

#ALL HOURS	
DIAGRAM_VALUES = ValueLookup(DiaMapDict,MasterDict)
testspc='O3'
print 'OH VOC List (not real right now)'
print DIAGRAM_VALUES['OH_VOC_DIAGRAM']

print '%s timeseries' % testspc
print 'just chem right now'
print CreateTSChemArrays(MasterDict['timeseries'])[testspc]
#When proc are done, create a ts master by running process,
#since most lines are there.  Then append chemistry from the
#chem arrays
#
#TSMaster = CreateTSProcArrays(MasterDict['processes'])['H_Adv'][testspc]
#for spc, vals in CreateTSChemArrays(MasterDict['timeseries'])
#	TSMaster[spc]['chemistry']=vals
#print TSMaster[testspc]