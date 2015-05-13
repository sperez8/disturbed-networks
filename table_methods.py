
#library imports
import sys
import os
import argparse
import numpy as np
from tabulate import tabulate

from network_simulation import get_network_fullnames

BEGINNING = '''\\begin{table}
\caption[]{}
\label{tab:label}
\centering\n'''

END = '''\n\end{table}'''

def sample_sequence(net_path, networkNames, inputFolder, inputFileEnd):
	'''makes an OTU table with avg depth and othe features per OTU'''

	networks,treatments = get_network_fullnames(networkNames)

	otuTable = {}
	for n in networks:
		otuTable[n] = np.loadtxt(os.path.join(inputFolder,n.replace('BAC_','')+inputFileEnd), dtype='S1000')

	samples = []

	final = ''
	for location,treatments in networkNames.iteritems():
		for t in treatments:
			table = otuTable[location+'_'+t]
			sampleNames = table[0,1:-1]
			sampleCounts = table[1:-1,1:-1].astype(np.int).sum(axis=0)
			samplesAdd = np.concatenate((sampleNames[np.newaxis],sampleCounts[np.newaxis]),axis=0)
			if samples == []:
				samples = samplesAdd
			else:
				samples = np.concatenate((samples,samplesAdd),axis=1)
			text = convert(samples.T,header=True, add_headers=['Sample Id','Number of sequences'])
			caption = "[Number of sequences recovered for samples in ecozone {0} with treatment {1}]".format(location.split('_')[1],t)
			caption = caption + caption.replace("[","{").replace("]","}")
			newtext = text.replace('caption','caption'+caption)
			final += newtext
			final += '\n\n'

	f = open('dump.txt','w')
	f.write(final)
	f.close()
	return None

def convert_file(fileName, header=False, rows=False):
	table = np.loadtxt(fileName, delimiter='\t', dtype='S1000')
	convert(table,header,rows)
	return None

def convert(table, header=False, rows=False, add_headers= []):
	if header:
		headerNames = list(table[0,:])
		table = table[1:,:]
	if add_headers and header:
		headerNames = add_headers
	tableText = tabulate(table, headers=headerNames, tablefmt="latex")
	finalText = collect_pieces(tableText)
	# f = open('dump.txt','w')
	# f.write(finalText)
	# f.close()
	return finalText

def collect_pieces(table):
	text = BEGINNING
	text += table
	text += END
	return text