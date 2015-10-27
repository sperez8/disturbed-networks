'''
created  01/17/2014

by sperez

Runs attack/ex tinction simulations on networks
'''

#library imports
import sys
import os
import argparse
import numpy as np
import prettyplotlib as ppl
import math
import powerlaw 
from decimal import Decimal

# prettyplotlib imports 
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib_venn import venn3, venn3_circles, venn2
from prettyplotlib import brewer2mpl

# _cur_dir = os.path.dirname(os.path.realpath(__file__))
# _root_dir = os.path.dirname(_cur_dir)
# sys.path.insert(0, _root_dir)

import networkx as nx
from make_network import import_graph
import network_measures as nm

RANDSEED = 2
np.random.seed(RANDSEED)
DPI = 200 #resolution of plot #low for testing
RAND_NAME = 'random_network_size_of_'
SCALE_NAME = 'scalefree_network_size_of_'
FILTER_NON_OTUS = True
MARKER_SIZE = 200
NUM_BINS = 30.0

NOT_A_NODE_VALUE = 'NA'

TAXONOMY = ["kingdom","phylum","class","order","family","genus","species","subspecies","subsubspecies"]

INDVAL_P_CUTOFF = 0.01
INDVAL_CUTOFF = 0.6

BIG_FIG_WIDTH = 12
BIG_FIG_HEIGHT = 16

SMALL_FIG_WIDTH = 5
SMALL_FIG_HEIGHT = 5

TITLE_FONT = 20

OM_COLORS = {"OM0":"#238b45",
			"OM1":"#fd8d3c",
			"OM2":"#e31a1c",
			"OM3":"#800026"}

# OM_COLORS = {"All":"black",
# 			"Mod1":"#01665e",
# 			"Mod2":"#35978f",
# 			"Mod3":"#80cdc1"}

STRUCTURE_METRICS = [nm.number_of_nodes, 
					nm.number_of_edges,
					nm.number_of_components,
					nm.size_of_big_components,
					nm.number_of_nodes_of_largest_connected_component, 
					nm.number_of_edges_of_largest_connected_component,
					nm.diameter_of_largest_connected_component,
					nm.average_degree, 
					nm.connectance, 
					nm.global_clustering_coefficient,
					nm.fraction_of_possible_triangles,
					nm.size_of_largest_clique,
					nm.average_path_on_largest_connected_component,
					nm.degree_assortativity,
					nm.correlation_of_degree_and_betweenness_centrality,
					]

INPUT_METRICS = []
				# nm.richness,
				# nm.shannon_diversity,
				# ]
OTU_METRICS = [#nm.correlation_of_degree_and_depth,
			#	nm.correlation_of_edge_depth,
			#	nm.compute_modularity_horizon,
				]
MEASURES = [nm.node_degrees,
			nx.betweenness_centrality, 
			nx.clustering, 
			nm.in_largest_connected_component,
			nm.node_modularity,
			]

MODULE_METRICS = [nm.number_of_nodes, 
					nm.number_of_edges,
					nm.diameter_of_largest_connected_component,
					nm.average_degree, 
					nm.connectance, 
					nm.global_clustering_coefficient,
					nm.fraction_of_possible_triangles,
					nm.average_path_on_largest_connected_component,
					nm.degree_assortativity,
					nm.correlation_of_degree_and_betweenness_centrality,
					]

MODULE_OTU_METRICS = []
					# nm.correlation_of_degree_and_depth,
					# nm.correlation_of_edge_depth,
					#FIX MEEE
					#nm.avg_depth,
					#]



def get_graph(nodeFile, edgeFile,edgetype):
	'''imports the node and edge file and makes the graph'''
	G = import_graph(nodeFile,edgeFile,edgetype,FILTER_NON_OTUS)
	return G

def get_multiple_graphs(networks, path, edgetype, add_random, add_scalefree, LCC=False):
	'''makes multiple graphs from names of networks and a file path'''
	graphs = {}
	for netName in networks:
		nodeFile = os.path.join(path,netName+'_nodes.txt')
		edgeFile = os.path.join(path,netName+'_edges.txt')
		G = get_graph(nodeFile,edgeFile,edgetype)
		if LCC:
			G = nm.get_LCC(G)
			print "keeping only connected component"
		graphs[netName] = G
		print 'Made the networkx graph {0} with N = {1}, E = {2}.'.format(netName,G.number_of_nodes(),G.number_of_edges())
		
		##adding random graph for comparaison
		if add_random:
			M = nx.number_of_edges(G)
			N = nx.number_of_nodes(G)
			H = nx.gnm_random_graph(N,M,seed=RANDSEED)
			if LCC:
				H = nm.get_LCC(H)
			graphs[RAND_NAME+netName] = H
		if add_scalefree:
			N = nx.number_of_nodes(G)
			H = nx.scale_free_graph(N,seed=RANDSEED)
			UH = H.to_undirected()
			UH = nx.Graph(UH)
			if LCC:
				UH = nm.get_LCC(UH)		
			graphs[SCALE_NAME+netName] = UH
	return graphs

def get_network_fullnames(networkNames):
	networks = []
	key = networkNames.keys()[0]
	if networkNames[key] == []:
		return networkNames.keys(),None
	for location,treatments in networkNames.iteritems():
		location, treatments
		for t in treatments:
			networks.append(location+'_'+t)
	return networks,treatments

def load_samples_info(samplesFile):
	samplesTable = np.loadtxt(samplesFile, comments=None, delimiter='\t', dtype='S1000')
	return samplesTable

def get_info_per_samples(samplesFile, samples, feature):
	'''gets a list of samples and an ecological feature,
	 and returns a dicttionary of sample values for the desired feature'''
	sampleInfo = {}
	table = load_samples_info(samplesFile)
	for s in samples:
		row = np.where(table[:,0]==s)[0][0]
		column = np.where(table[0,:]==feature)[0][0]
		sampleInfo[s] = table[row,column]
	return sampleInfo

def get_ind(table, otu):
	ind = 0
	indrowindex = np.where(table[:,0]==otu)
	if indrowindex[0]: #check that otu is in fact in that table, otherwise ignore and return 0
		indrow = table[indrowindex[0][0],:]
		pcol = np.where(table[0,:]=="output$pval[rownames(output$indval)]")[0][0]
		valcol = np.where(table[0,:]=="output$indcls[rownames(output$indval)]")[0][0]
		clustercol = np.where(table[0,:]=="output$maxcls[rownames(output$indval)]")[0][0]
		if float(indrow[pcol]) < INDVAL_P_CUTOFF:
			if float(indrow[valcol])>INDVAL_CUTOFF:
				return indrow[clustercol]
	else:
		print 'not found:',otu
	return ind

def make_OTU_feature_table(net_path, networkNames, inputFolder, inputFileEnd, indvalFolder, indvalFileEnd, samplesFile, features, path, featureFile,edgetype, factor):
	'''makes an OTU table with avg depth and othe features per OTU'''

	networks,treatments = get_network_fullnames(networkNames)
	if len(MEASURES)>=1:
		graphs = get_multiple_graphs(networks, net_path, edgetype, False, False)

	otuTable = {}
	for n in networks:
		otuTable[n] = np.loadtxt(os.path.join(inputFolder,n.replace('BAC_','')+inputFileEnd), dtype='S1000')

	indTable = {}
	for n in networkNames.keys():
		indTable[n] = np.loadtxt(os.path.join(indvalFolder,n.replace('BAC_','')+indvalFileEnd), dtype='S1000')

	header = ['OTUs','Abundance','IndforCluster']
	headerStart = len(header)
	for f in features:
		header.append(f+ ' avg')
		header.append(f+ ' std')
	header.extend([m.__name__.replace('_',' ').capitalize() for m in MEASURES])
	header.extend(TAXONOMY)
	tax_index = len(TAXONOMY)

	for location,treatments in networkNames.iteritems():
		for t in treatments:
			abundances = otuTable[location+'_'+t]
			sampleNames = abundances[0,1:-1]
			sampleCounts = abundances[1:-1,1:-1].astype(np.float).sum(axis=0)
			featureTable = np.zeros(shape=(abundances.shape[0]-1,headerStart+len(features)*2+len(MEASURES)+tax_index), dtype='S1000')
			featureTable[0,:] = np.array(header)
			for r,row in enumerate(abundances[1:-1,]):
				otu = row[0]
				ab = row[1:-1].astype(np.float)/sampleCounts
				tax = row[-1].replace('(100)','')
				tax = np.asarray(tax.split(';'))[0:tax_index]
				ind = get_ind(indTable[location], otu)
				featureTable[r+1][0]=otu
				featureTable[r+1][1]=np.mean(ab)
				featureTable[r+1][2]=ind
				featureTable[r+1][-tax_index:]=tax
			for i,f in enumerate(features):
				print "For input table from zone {0} treatment {1} calculating feature {2}".format(location,t,f)
				fdist = get_info_per_samples(samplesFile, sampleNames, f)				
				for r,row in enumerate(abundances[1:-1,]):
					ab = row[1:-1].astype(np.float)/sampleCounts
					fcount = []
					stdCount = []
					featureValues = [float(fdist[s]) for s in sampleNames]
					avg = np.average(featureValues, weights = ab)
					std = math.sqrt(np.average((featureValues-avg)**2, weights=ab))
					featureTable[r+1][i+headerStart]=avg
					featureTable[r+1][i+headerStart+1]=std

			column_bias = headerStart+len(features)*2
			for i,m in enumerate(MEASURES):
				print "For input table from zone {0} treatment {1} calculating measure {2}".format(location,t,m.__name__)
				col = i+column_bias
				G = graphs[location+'_'+t]
				if "modul" in m.__name__:
					values = m(G,factor=factor)
				else:
					values = m(G)
				for r,row in enumerate(abundances[1:-1,]):
					otu = row[0]
					measureValue = 0
					if otu in G.nodes():
						measureValue = values[otu]
					elif 'OTU-'+otu in G.nodes():
						measureValue = values['OTU-'+otu]
					else:
						measureValue = NOT_A_NODE_VALUE				
					featureTable[r+1][col]=measureValue

			fileName = featureFile+'_{0}_{1}_{2}_factor{3}.txt'.format(edgetype,location,t,factor)
			tableFile = os.path.join(path,fileName)
			print "Saving table: ",tableFile

			np.savetxt(tableFile, featureTable, delimiter="\t", fmt='%s')
	return None





##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
#############                                                        #########
#################                                                 ############
####################                                    ######################
############################                       ###########################
#################################             ################################
###################################         ##################################
#######################################   ####################################
##############################################################################

def plot_degree_distribution_per_treatment(net_path, networkNames, figurePath, plot_sequence, edgetype):
	networks,treatments = get_network_fullnames(networkNames)
	graphs = get_multiple_graphs(networks, net_path, edgetype, False, False)
	data = {}
	locations = networkNames.keys()
	locations.sort(reverse=True)
	print locations, treatments

	mins = []
	maxs = []

	fig, axes = plt.subplots(len(locations),len(treatments))
	iterable = []
	for i,r in enumerate(locations):
		for j,c in enumerate(treatments):
			if len(locations)>1:
				iterable.append((axes[i][j],r,c))
			else:
				iterable.append((axes[j],r,c))

	dumpfile = open('\Users\Sarah\Desktop\LTSPnetworks\dumpfit.txt','w')
	colors = {treatment: OM_COLORS[treatment] for i,treatment in enumerate(treatments)}
	for ax,location,t in iterable:
		G = graphs[location+'_'+t]
		if locations.index(location) == 0:
			ax.set_title(t)
		if treatments.index(t) == 0:
			ax.set_ylabel(location.split('_')[1])
		if locations.index(location) == 0 and treatments.index(t) == 0:
			ax.set_xlabel('Node degree')
		ax.set_yscale('log')
		ax.set_xscale('log')
		# N = G.number_of_nodes()
		# ds = [] #each degree
		# fds = [] #each degree's frequency
		# degrees = sorted(nx.degree(G).values(),reverse=True)
		# # for d in set(degrees):
		# # 	ds.append(d)
		# # 	fds.append(float(degrees.count(d))/N)
		# # min_y = min(min(fds),min_y)

		degrees = sorted(nx.degree(G).values(),reverse=True)
		
		data = degrees
		if max(degrees)>2:
			fit = powerlaw.Fit(data,discrete=True,xmin=1)
			fit.plot_pdf(ax=ax,color=colors[t]) #plot probability distribution function
			dumpfile.write(location+','+t)
			dumpfile.write('\n')
			#print ','.join(['s',str(fit.power_law.sigma), str(fit.power_law.alpha)])
			fit_exp = fit.stretched_exponential
			beta,Lambda = fit.stretched_exponential.beta, fit.stretched_exponential.Lambda
			fit.stretched_exponential.plot_pdf(ax=ax, color=colors[t],linestyle='--',linewidth=2)
			#print 'beta,lambda', beta, Lambda
			for dist in ['exponential','power_law','lognormal']:	#,'truncated_power_law']:
				R, p = fit.distribution_compare(dist,'stretched_exponential',nested=None)
				print dist, R, p
			dumpfile.write('beta = '+str(beta)+', lambda = '+str(Lambda))
			dumpfile.write('\n')
			dumpfile.write('KS = '+str(fit_exp.KS())+',  a2 = '+str(fit_exp.Asquare))
			dumpfile.write('\n')
			dumpfile.write('\n')

			#print fit_exp.sigma
			i,j = ax.get_ylim()
			mins.append(i)
			maxs.append(i)

		#ax.plot(ds,[math.pow(d,-k) for d in ds],color=colors[t])
	# for ax in axes:
	# 	ax.set_ylim(min(mins),max(maxs))

	
	#add legend manually
	for t in treatments:
		ax.plot([], [], linestyle='-',label=str(t),color=colors[t])


	#ax.set_ylim([min_y,1])
	title = 'Degree distribution of '+' '.join(location)+' network with '+edgetype+' type of edges'
	figureTitle = fig.suptitle(title, horizontalalignment='center', fontsize=TITLE_FONT)
	lgd = ppl.legend(bbox_to_anchor=(1.05, 1), loc=2)

	fig.set_size_inches(BIG_FIG_HEIGHT,BIG_FIG_WIDTH)
	fig.savefig(figurePath, dpi=DPI, bbox_extra_artists=(lgd,), bbox_inches='tight')
	print "Saving the figure file: ", figurePath
	return None

def make_ecological_table(net_path, networkNames, filePath, edgetype, inputFolder, inputFileEnd):
	networks,treatments = get_network_fullnames(networkNames)
	#graphs = get_multiple_graphs(networks,net_path,edgetype, False, False)

	otuTable = {}
	for n in networks:
		otuTable[n] = np.loadtxt(os.path.join(inputFolder,n.replace('BAC_','')+inputFileEnd), dtype='S1000')

	table = np.zeros(shape=(len(INPUT_METRICS)+2, len(networkNames)*len(treatments)+1), dtype='S1000')
	i,j = 0,1 # i is row, j is column
	column = ['Zones','Treatments']
	column.extend([sm.__name__.replace('_',' ').capitalize() for sm in INPUT_METRICS])
	table[:,0]=column
	for location,treatments in networkNames.iteritems():
		table[i,j]=location
		for t in treatments:
			i+=1
			table[i,j]=t
			for im in INPUT_METRICS:
				print "For input table from zone {0} treatment {1} measuring {2}".format(location,t,im.__name__)
				i+=1
				S = otuTable[location+'_'+t]
				table[i,j]=im(S)
			j+=1
			i=0

	print "Saving table: ", filePath

	np.savetxt(filePath, table, delimiter="\t", fmt='%s')
	return None

def get_taxonomic_levels(featurePath,featureFile,location,treatments,tax_level,netcol):
	taxonomies = []
	edgetype = 'pos'
	for t in treatments:
		featureTableFile = os.path.join(featurePath,featureFile+'_{0}_{1}_{2}.txt'.format(edgetype,location,t))
		featureTable = np.loadtxt(featureTableFile,delimiter='\t', dtype='S1000')
		#get tax levels of only levels present in network
		netfeatureTable = featureTable[np.where(featureTable[:,netcol]!=NOT_A_NODE_VALUE)]
		taxcol = np.where(netfeatureTable[0,:]==tax_level)[0][0]
		taxonomies.extend(list(np.unique(netfeatureTable[1:,taxcol])))
	taxonomies = list(set(taxonomies))	
	taxonomies.sort()
	return taxonomies


def calculate_taxonomic_representation(net_path, networkNames, figurePath, featurePath, featureFile, tax_level, percentNodes, bcMinValue):
	edgetype = 'pos'
	colName = nx.betweenness_centrality.__name__.replace('_',' ').capitalize()
	BCtaxa = {}
	notBCtaxa ={}

	for location,treatments in networkNames.iteritems():
		BCtaxa = {}
		alltaxa = {}
		for t in treatments:
			featureTableFile = os.path.join(featurePath,featureFile+'_{0}_{1}_{2}.txt'.format(edgetype,location,t))
			featureTable = np.loadtxt(featureTableFile,delimiter='\t', dtype='S1000')
			centcol = np.where(featureTable[0,:]==colName)[0][0]
			taxcol = np.where(featureTable[0,:]==tax_level)[0][0]
			print featureTable.shape
			nodes = featureTable[np.where(featureTable[:,centcol]!=NOT_A_NODE_VALUE)]
			print nodes.shape
			nodes = nodes[np.where(nodes[:,taxcol]!="unclassified")]
			nodes = nodes[np.where(nodes[:,taxcol]!="uncultured")]
			nodes = nodes[1:,0]
			bcvalues = featureTable[1:,centcol]
			bcvalues = bcvalues[np.where(bcvalues!=NOT_A_NODE_VALUE)]
			bcvalues= list([float(k) for k in bcvalues])
			bcvalues.sort(reverse=True)
			if percentNodes<1:
				#cutoff = max(float(bcvalues[int(percentNodes*float(len(bcvalues)))-1]), bcMinValue)
 				cutoff = max(float(percentNodes*max(bcvalues)), bcMinValue)
			else:
				cutoff = float(bcvalues[int(percentNodes*float(len(bcvalues)))-1])

			for n in nodes:
				row = nm.findRow(n,featureTable)
				taxonlevel = featureTable[row][taxcol]
				value = float(featureTable[row][centcol])
				if value != NOT_A_NODE_VALUE and value >= cutoff:
					if taxonlevel in BCtaxa.keys():
						BCtaxa[taxonlevel].append(n)
					else:
						BCtaxa[taxonlevel] = [n]

				if taxonlevel in alltaxa.keys():
					alltaxa[taxonlevel].append(n)
				else:
					alltaxa[taxonlevel] = [n]

		taxonomies = alltaxa.keys()
		taxonomies.sort()

		Prob_all_seen_greater_1 = 1
		representation = np.zeros(shape=(len(taxonomies)+1,4), dtype='S1000')
		#representation = np.zeros(shape=(len(taxonomies)+1,3), dtype='S1000')
		representation[1:,0]=np.array(taxonomies)
		representation[0,:]=np.array([tax_level.capitalize(),"Number of taxa","Number of central taxa","Over-representation (p-value)"])
		#representation[0,:]=np.array([tax_level.capitalize(),"Number of taxa","Number of central taxa"])
		m = sum([len(v) for v in alltaxa.values()]) #total number of OTUs
		n = sum([len(v) for v in BCtaxa.values()]) #total number of central OTUs
		
		for i,taxonomy in enumerate(taxonomies):
			print taxonomy
			mi = len(alltaxa[taxonomy]) #number of OTUs with that taxonomy
			if taxonomy in BCtaxa.keys():
				yi = len(BCtaxa[taxonomy]) #number of central OTUs with that taxonomy
			else:
				yi = 0
			print m, n, mi, yi
			representation[i+1,0]=taxonomy
			representation[i+1,1]=mi
			representation[i+1,2]=yi
			if yi>0:
				prob = prob_hypergeometric(m,n,mi,yi)
				representation[i+1,3]=round(prob,2)
				#Prob_all_seen_greater_1 *= prob_hypergeometric(m,n,mi,yi,atleastone=True)
			else:
				representation[i+1,3]="None"

		#print Prob_all_seen_greater_1
		#save representation in a table
		f = open(os.path.join(net_path,figurePath,"representation_{0}_{1}.txt".format(location.split('_')[1],tax_level)),'w')
		np.savetxt(f, representation, delimiter="\t", fmt='%s')
		f.close()

	return None

def prob_hypergeometric(m,n,mi,yi,atleastone=False):
	p = 0
	if atleastone:
		for j in range(1,min(n,mi)+1):
			p += nCk(mi,j)*nCk(m-mi,n-j)/Decimal(nCk(m,n))
	else:
		if mi==yi:mi+=1
		for j in range(yi,min(n,mi)):
			p += nCk(mi,j)*nCk(m-mi,n-j)/Decimal(nCk(m,n))
	return p

#n choose k
def nCk(n,k):
    f = math.factorial
    return f(n) / f(k) / f(n-k)

def plot_venn_diagram(net_path, networkNames, figurePath, featurePath, featureFile, tax_level, percentNodes, bcMinValue):
	edgetype = 'pos'
	networks,treatments = get_network_fullnames(networkNames)
	colName = nx.betweenness_centrality.__name__.replace('_',' ').capitalize()
	taxonomy = TAXONOMY
	taxonomy.pop(0)
	fig, axes = plt.subplots(len(taxonomy))
	totaltax = {tax:0 for tax in taxonomy}
	unclassified = {tax:0 for tax in taxonomy}
	#print totaltax, unclassified

	for ax,tax_level in zip(axes, taxonomy):
		taxaSeen = {}
		for location,treatments in networkNames.iteritems():
			taxaSeen[location] = []
			for t in treatments:
				featureTableFile = os.path.join(featurePath,featureFile+'_{0}_{1}_{2}.txt'.format(edgetype,location,t))
				featureTable = np.loadtxt(featureTableFile,delimiter='\t', dtype='S1000')
				centcol = np.where(featureTable[0,:]==colName)[0][0]
				taxcol = np.where(featureTable[0,:]==tax_level)[0][0]
				nodes = featureTable[np.where(featureTable[:,centcol]!=NOT_A_NODE_VALUE)][1:,0]
				bcvalues = featureTable[1:,centcol]
				bcvalues = bcvalues[np.where(bcvalues!=NOT_A_NODE_VALUE)]
				bcvalues= list([float(k) for k in bcvalues])
				bcvalues.sort(reverse=True)
				if percentNodes<1:
 					#cutoff = max(float(bcvalues[int(percentNodes*float(len(bcvalues)))-1]), bcMinValue)
 					cutoff = max(float(percentNodes*max(bcvalues)), bcMinValue)
				else:
					cutoff = float(bcvalues[int(percentNodes*float(len(bcvalues)))-1])
				for n in nodes:
					row = nm.findRow(n,featureTable)
					taxonlevel = featureTable[row][taxcol]
					value = float(featureTable[row][centcol])
					if value != NOT_A_NODE_VALUE and value >= cutoff:
						taxaSeen[location].append(taxonlevel)
			unclassified[tax_level] += taxaSeen[location].count("unclassified")
			unclassified[tax_level] += taxaSeen[location].count("uncultured")
			#print unclassified, len(taxaSeen[location]),taxaSeen[location], totaltax[tax_level]
			totaltax[tax_level] = totaltax[tax_level] + len(taxaSeen[location])

			taxaSeen[location] = set(taxaSeen[location])
			if "unclassified" in taxaSeen[location]:
				taxaSeen[location].remove("unclassified")	
			if "uncultured" in taxaSeen[location]:
				taxaSeen[location].remove("uncultured")

		#print tax_level
		#print taxaSeen
		#ax.set_ylabel(tax_level)

		total = sum([len(taxa) for taxa in taxaSeen.values()])
		circleLabels = [k.split('_')[1] for k in taxaSeen.keys()]
		if total == 0:
			circleLabels=[]
		if len(networkNames)==3:
			v = venn3(subsets=taxaSeen.values(), set_labels = circleLabels, set_colors=('b','r','g'), alpha=0.3, ax=ax)#axlabel = tax_level)
		else:
			v = venn2(subsets=taxaSeen.values(), ax = ax)
	
	#save number of unclassified in a a table
	f = open(os.path.join(net_path,figurePath,"prop_unclassifed_OTUs_venn_{0}_{1}.txt".format(percentNodes,bcMinValue )),'w')
	f.write('\t'.join(['Taxonomic level','Number of unclassified OTUs','Total number of OTUs','Proportion of unclassified']))
	f.write('\n')
	for t in taxonomy:
		f.write('\t'.join([str(x) for x in [t, unclassified[t], totaltax[t], unclassified[t]/float(totaltax[t])]]))
		f.write('\n')
	f.close()


	title = "Venn diagram of econoze's central OTUs classified by taxonomic level"
	figureTitle = fig.suptitle(title, horizontalalignment='center', fontsize=TITLE_FONT)

	fig.set_size_inches(15,5*len(taxonomy))
	figureFile = os.path.join(net_path,figurePath,'Venn_diagram_plot_'+'_'.join([k.split('_')[1] for k in networkNames.keys()])+'_'+edgetype+'_'+str(percentNodes)+'_'+str(bcMinValue)+'.png')
	fig.savefig(figureFile, dpi=DPI,bbox_inches='tight')
	print "Saving the figure file: ", figureFile

	return None


##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
#############                                                        #########
#################                                                 ############
####################                                    ######################
############################                       ###########################
#################################             ################################
###################################         ##################################
#######################################   ####################################
##############################################################################


def plot_venn_otus_diagram(net_path, networkNames, figurePath, featurePath, featureFile, percentNodes, bcMinValue):
	edgetype = 'pos'
	networks,treatments = get_network_fullnames(networkNames)
	colName = nx.betweenness_centrality.__name__.replace('_',' ').capitalize()

	fig, ax = plt.subplots(1)

	otuSeen = {}
	for location,treatments in networkNames.iteritems():
		otuSeen[location] = []
		#centralities = {}
		for t in treatments:
			featureTableFile = os.path.join(featurePath,featureFile+'_{0}_{1}_{2}.txt'.format(edgetype,location,t))
			featureTable = np.loadtxt(featureTableFile,delimiter='\t', dtype='S1000')
			centcol = np.where(featureTable[0,:]==colName)[0][0]
			nodes = featureTable[np.where(featureTable[:,centcol]!=NOT_A_NODE_VALUE)][1:,0]
			taxcol = 0
			bcvalues = featureTable[1:,centcol]
			bcvalues = bcvalues[np.where(bcvalues!=NOT_A_NODE_VALUE)]
			bcvalues= list([float(k) for k in bcvalues])
			bcvalues.sort(reverse=True)
			if percentNodes<1:
				#cutoff = max(float(bcvalues[int(percentNodes*float(len(bcvalues)))-1]), bcMinValue)
				cutoff = max(float(percentNodes*max(bcvalues)), bcMinValue)
			else:
				cutoff = float(bcvalues[int(percentNodes*float(len(bcvalues)))-1])
			for n in nodes:
				row = nm.findRow(n,featureTable)
				taxonlevel = featureTable[row][taxcol]
				value = float(featureTable[row][centcol])
				if value != NOT_A_NODE_VALUE and value >= cutoff:
					otuSeen[location].append(taxonlevel)

		otuSeen[location] = set(otuSeen[location])

	#print otuSeen

	total = sum([len(taxa) for taxa in otuSeen.values()])
	circleLabels = [k.split('_')[1] for k in otuSeen.keys()]
	if total == 0:
		circleLabels=[]
	if len(networkNames)==3:
		v = venn3(subsets=otuSeen.values(), set_labels = circleLabels, set_colors=('b','r','g'), alpha=0.3, ax=ax)
	else:
		v = venn2(subsets=otuSeen.values(), ax = ax)
		
	title = "Venn diagram of econoze's central OTUs"
	figureTitle = fig.suptitle(title, horizontalalignment='center', fontsize=TITLE_FONT)

	fig.set_size_inches(5,5)
	figureFile = os.path.join(net_path,figurePath,'Venn_diagram_plot_OTUs_'+'_'.join(networkNames.keys())+'_'+edgetype+'_'+str(percentNodes)+'_'+str(bcMinValue)+'.png')
	fig.savefig(figureFile, dpi=DPI,bbox_inches='tight')
	print "Saving the figure file: ", figureFile

	return None



def centrality_plot(net_path, networkNames, figurePath, featurePath, featureFile, tax_level, percentNodes, bcMinValue):
	networks,treatments = get_network_fullnames(networkNames)
	edgetype = 'pos'
	colName = nx.betweenness_centrality.__name__.replace('_',' ').capitalize()

	if tax_level not in TAXONOMY:
		tax_level = TAXONOMY[1] #phylum

	fig, axes = plt.subplots(len(treatments))
	netNames = treatments
	max_y = 0


	for location,treatments in networkNames.iteritems():
		taxaSeen = []
		centralities = {}
		for t in treatments:
			featureTableFile = os.path.join(featurePath,featureFile+'_{0}_{1}_{2}.txt'.format(edgetype,location,t))
			featureTable = np.loadtxt(featureTableFile,delimiter='\t', dtype='S1000')
			centcol = np.where(featureTable[0,:]==colName)[0][0]
			taxcol = np.where(featureTable[0,:]==tax_level)[0][0]
			nodes = featureTable[np.where(featureTable[:,centcol]!=NOT_A_NODE_VALUE)][1:,0]
			taxonomies = get_taxonomic_levels(featurePath,featureFile,location, treatments, tax_level, centcol)
			centralities[t] = [[] for tax in taxonomies]
			bcvalues = featureTable[1:,centcol]
			bcvalues = bcvalues[np.where(bcvalues!=NOT_A_NODE_VALUE)]
			bcvalues= list([float(k) for k in bcvalues])
			bcvalues.sort(reverse=True)
			cutoff = max(float(bcvalues[int(percentNodes*float(len(bcvalues)))-1]), bcMinValue)
			for n in nodes:
				row = nm.findRow(n,featureTable)
				taxonlevel = featureTable[row][taxcol]
				value = float(featureTable[row][centcol])
				if value != NOT_A_NODE_VALUE and value >= cutoff:
					max_y = max(max_y,float(value))
					taxaSeen.append(taxonlevel)
					centralities[t][taxonomies.index(taxonlevel)].append(value)
		taxaSeen = set(taxaSeen)

		for i,tax in enumerate(taxonomies):
			if tax not in taxaSeen:
				taxonomies.pop(i)
				for t in treatments:
					centralities[t].pop(i)

		for ax,t in zip(axes,treatments):
			labels = [tax+' ('+str(len(centralities[t][i]))+')' for i,tax in enumerate(taxonomies)]
			ppl.boxplot(ax, centralities[t])
			xticks = ax.set_xticklabels(labels,rotation=15,fontsize=14)
			ax.set_ylabel('Treatment '+t,fontsize=14)

		for ax in axes:
			ax.set_autoscaley_on(False)
			m,x = ax.get_ylim()
			ax.set_ylim([0,max_y])
			#ax.set_yscale('log')
			ax.grid()

		title = "Centrality of taxonomic level '{0}' present in network {1}".format(tax_level,location,t)
		figureTitle = fig.suptitle(title, horizontalalignment='center', fontsize=TITLE_FONT)

		fig.set_size_inches(2.5*len(taxonomies),6*len(treatments))
		figureFile = os.path.join(net_path,figurePath,'centrality_plot_'+location+'_'+edgetype+'_'+str(percentNodes)+'_'+str(bcMinValue)+'.png')
		fig.savefig(figureFile, dpi=DPI,bbox_inches='tight')
		print "Saving the figure file: ", figureFile
	return None




def keystone_quantitative_feature_plot(net_path, networkNames, figurePath, featurePath, featureFile, features, percentNodes, bcMinValue):
	networks,treatments = get_network_fullnames(networkNames)
	edgetype = 'pos'
	colName = nx.betweenness_centrality.__name__.replace('_',' ').capitalize()
	#modName = nm.node_modularity.__name__.replace('_',' ').capitalize()

	fig, axes = plt.subplots(len(features))
	netNames = treatments

	for location in networkNames.keys():
		for ax,f in zip(axes,features):
			featureValues = []
			for i,t in enumerate(treatments):
				featureValues.append([])
				featureTableFile = os.path.join(featurePath,featureFile+'_{0}_{1}_{2}.txt'.format(edgetype,location,t))
				featureTable = np.loadtxt(featureTableFile,delimiter='\t', dtype='S1000')
				centcol = np.where(featureTable[0,:]==colName)[0][0]
				nodes = featureTable[np.where(featureTable[:,centcol]!=NOT_A_NODE_VALUE)][1:,0]
				featcol = np.where(featureTable[0,:]==f)[0][0]
				#modcol =  np.where(featureTable[0,:]==modName)[0][0]
				bcvalues = featureTable[1:,centcol]
				bcvalues = bcvalues[np.where(bcvalues!=NOT_A_NODE_VALUE)]
				bcvalues= list([float(k) for k in bcvalues])
				bcvalues.sort(reverse=True)
				cutoff = max(float(bcvalues[int(percentNodes*float(len(bcvalues)))-1]), bcMinValue)
				for n in nodes:
					row = nm.findRow(n,featureTable)
					bc = float(featureTable[row][centcol])
					#mod = int(featureTable[row][modcol])
					if bc != NOT_A_NODE_VALUE and bc >= cutoff: # and mod == 0:
						value = featureTable[row][featcol]
						featureValues[i].append(float(value))


			labels = [t+' ('+str(len(featureValues[i]))+')' for i,t in enumerate(treatments)]
			ppl.boxplot(ax, featureValues)
			xticks = ax.set_xticklabels(labels,rotation=15,fontsize=14)
			ax.set_ylabel(f,fontsize=14)

		for ax in axes:
			#ax.set_yscale('log')
			ax.grid()

		title = "Properties of high betweenness centrality OTUs present in network {0}".format(location)
		figureTitle = fig.suptitle(title, horizontalalignment='center', fontsize=TITLE_FONT)

		fig.set_size_inches(2.5*len(treatments),6*len(features))
		figureFile = os.path.join(net_path,figurePath,'high_bc_feature_plot_'+location+'_'+edgetype+'_'+str(percentNodes)+'_'+str(bcMinValue)+'.png')
		fig.savefig(figureFile, dpi=DPI,bbox_inches='tight')
		print "Saving the figure file: ", figureFile
	return None

def plot_diff_centralities(net_path, networkNames, figurePath, featurePath, featureFile, percentNodes, measures):
	edgetype = 'pos'
	networks,treatments = get_network_fullnames(networkNames)
	graphs = get_multiple_graphs(networks,net_path,edgetype, False, False)

	locations = networkNames.keys()
	locations.sort(reverse=True)

	measureNames = [m.__name__ for m in measures]
	measureAcronym = {m:m.split('_')[0].capitalize()[0]+m.split('_')[1].capitalize()[0] for m in measureNames}
	measurePrettyName = {m:m.replace('numpy','').replace('_',' ').capitalize() for m in measureNames}
#.__name__.replace('_',' ').replace('numpy','').capitalize()

	#for each network, we collect all BC, DC, CC, EC values and plot them together
	for location in locations:
		fig, axes = plt.subplots(len(measures),len(measures))
		iterable = []
		#first collect all the values
		centralities = {m.__name__:[] for m in measures}
		for t in treatments:
			G = graphs[location+'_'+t]
			for m in measures:
				measure = m.__name__
				values = [(n,v) for n,v in m(G).iteritems()]
				#sort by node name so that all centrality value line up
				values = sorted(values, key = lambda item: item[0])
				centralities[measure].extend(zip(*values)[1])

		for i,mx in enumerate(measureNames):
			for j,my in enumerate(measureNames):
				iterable.append((axes[i][j],mx,my))
		for ax,mx,my in iterable:
			if mx==my:
				ppl.hist(ax,centralities[my],grid='y',color='blue')
			else:
				ppl.scatter(ax,centralities[my],centralities[mx])
			#ax.tick_params(axis='x', which='both', labelbottom='off')
			#ax.tick_params(axis='y', which='both', labelleft='off')
			ax.tick_params(axis='both', which='major', labelsize=8)
			ax.tick_params(axis='both', which='minor', labelsize=8)
			ax.set_ylim(0,ax.get_ylim()[1])
			ax.set_xlim(0,ax.get_xlim()[1])
			ax.locator_params(nbins=5)

			if measureNames.index(mx)==0:
				ax.set_title(measurePrettyName[my],fontsize=12)
			else:
				pass
				#ax.tick_params(axis='y', which='both', labelleft='on')
			if measureNames.index(my)==0:
				ax.set_ylabel(measurePrettyName[mx],fontsize=12)
			else:
				pass
				#ax.tick_params(axis='x', which='both', labelbottom='on')
		#plt.show()

		#title = "Distribution of OTUs properties in all ecozones for treatment "+t
		#figureTitle = fig.suptitle(title, horizontalalignment='center', fontsize=TITLE_FONT)

		fig.set_size_inches(BIG_FIG_WIDTH,BIG_FIG_WIDTH)
		figureFile = os.path.join(net_path,figurePath,'compare_centralities_'+location+'_'+edgetype+'.png')
		fig.savefig(figureFile, dpi=DPI, bbox_inches='tight')
		print "Saving the figure file: ", figureFile

	return None



##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
#############                                                        #########
#################                                                 ############
####################                                    ######################
############################                       ###########################
#################################             ################################
###################################         ##################################
#######################################   ####################################
##############################################################################

def plot_scatter_bc(net_path, networkNames, figurePath, featurePath, featureFile, percentNodes, bcMinValue):
	edgetype = 'pos'
	networks,treatments = get_network_fullnames(networkNames)
	bcColName = nx.betweenness_centrality.__name__.replace('_',' ').capitalize()
	hzColName = "SoilHorizon avg"
	abColName = "Abundance"

	attributes = [hzColName]#abColName]
	atributeValues = {a:[] for a in attributes}
	locations = networkNames.keys()
	locations.sort(reverse=True)

	for attribute in attributes:
		fig, axes = plt.subplots(len(locations),len(treatments))
	 	iterable = []
		for i,r in enumerate(locations):
			for j,c in enumerate(treatments):
				if len(locations)>1:
					iterable.append((axes[i][j],r,c))
				else:
					iterable.append((axes[j],r,c))
		for ax,location,t in iterable:
			featureTableFile = os.path.join(featurePath,featureFile+'_{0}_{1}_{2}.txt'.format(edgetype,location,t))
			featureTable = np.loadtxt(featureTableFile,delimiter='\t', dtype='S1000')
			notBC = []
			BC = []
			atCol = np.where(featureTable[0,:]==attribute)[0][0]
			centcol = np.where(featureTable[0,:]==bcColName)[0][0]
			nodes = featureTable[np.where(featureTable[:,centcol]!=NOT_A_NODE_VALUE)][1:,0]
			bcvalues = featureTable[1:,centcol]
			bcvalues = bcvalues[np.where(bcvalues!=NOT_A_NODE_VALUE)]
			bcvalues= list([float(k) for k in bcvalues])
			bcvalues.sort(reverse=True)
# 			#cutoff = max(float(bcvalues[int(percentNodes*float(len(bcvalues)))-1]), bcMinValue)
 			cutoff = max(float(percentNodes*max(bcvalues)), bcMinValue)
			for n in nodes:
				row = nm.findRow(n,featureTable)
				value = float(featureTable[row][centcol])
				if value == NOT_A_NODE_VALUE:
					print "UH oh"
				if value >= cutoff:
					BC.append(float(featureTable[row][atCol]))
				else:
					notBC.append(float(featureTable[row][atCol]))
			if BC:
				minvalue = min(min(notBC),min(BC))
				maxvalue = max(max(notBC),max(BC))
			else:
				minvalue = min(notBC)
				maxvalue = max(notBC)

			binwidth=(maxvalue-minvalue)/NUM_BINS

			notBC.extend(BC) # want to plot histogram of whole data and color the BC ones on top
			print location,t
			if binwidth == 0:
				binwidth = float(1.0)/float(NUM_BINS)
				bins = np.arange(minvalue, maxvalue + binwidth*2, binwidth)
				ppl.hist(ax, notBC, color='black', grid='y', bins=bins, label = 'All OTUs')
				if len(BC)>0:
					ppl.hist(ax, BC, color='purple', grid='y', bins=bins, label = 'High BC OTUs')
				elif location != 'BAC_JP' and t != 'OM3':
					ppl.hist(ax,[0,0], color='purple', grid='y', label = "High BC OTUs")
			else:
				if attribute == hzColName:
					binwidth = float(1.0)/float(NUM_BINS)
				bins = np.arange(minvalue, maxvalue + binwidth, binwidth)
				ppl.hist(ax, notBC, color='black', grid='y', bins = bins, label = 'All OTUs')
				if len(BC)>0:
					ppl.hist(ax, BC, color='purple', grid='y', bins = bins, label = "High BC OTUs")
				elif location != 'BAC_JP' and t != 'OM3':
					ppl.hist(ax,[-1], color='purple', grid='y', label = "High BC OTUs")
				else: #HUGE HACK I KNOW
					ppl.hist(ax, notBC, color='purple', grid='y', bins = bins, label = 'High BC OTUs')
					ppl.hist(ax, notBC, color='black', grid='y', bins = bins)


			ax.locator_params(axis = "x",nbins=4)

			if attribute == abColName or attribute == bcColName:
				ax.set_yscale('symlog')
			if treatments.index(t)==0:
				ax.set_ylabel(location.split('_')[1])
			if locations.index(location)==0:
				ax.set_title(t)
			if attribute==hzColName:
				ax.set_xlim(1-binwidth,2+binwidth)

		lgd = ppl.legend(bbox_to_anchor=(1.05, 1), loc=2)
		
		fig.set_size_inches(BIG_FIG_HEIGHT,BIG_FIG_HEIGHT)
		figureFile = os.path.join(net_path,figurePath,'hist_BC_otus_'+attribute+'_'+'_'.join([l.split('_')[1] for l in locations])+'_'+edgetype+'_'+str(percentNodes)+'_'+str(bcMinValue)+'.png')
		fig.savefig(figureFile, dpi=DPI, bbox_extra_artists=([lgd]),bbox_inches='tight')
		print "Saving the figure file: ", figureFile

	return None

# def plot_scatter_bc(net_path, networkNames, figurePath, featurePath, featureFile, percentNodes, bcMinValue):
# 	edgetype = 'pos'
# 	networks,treatments = get_network_fullnames(networkNames)
# 	bcColName = nx.betweenness_centrality.__name__.replace('_',' ').capitalize()
# 	hzColName = "SoilHorizon avg"
# 	abColName = "Abundance"

# 	attributes = [bcColName,hzColName,abColName]
# 	locations = networkNames.keys()
# 	locations.sort(reverse=True)

# 	numSuperconnectors = []

# 	for t in treatments:
# 		ytitles = list(locations)
# 		xtitles = list(attributes)
# 		fig, axes = plt.subplots(len(attributes),len(locations))
# 		iterable = []
# 		for i,r in enumerate(attributes):
# 			for j,c in enumerate(locations):
# 				if len(attributes)>1:
# 					iterable.append((axes[i][j],r,c))
# 				else:
# 					iterable.append((axes[j],r,c))
# 		#ylim = {attribute:0 for attribute in attributes}
# 		#xlims = {attribute:(1,0) for attribute in attributes}
# 		for ax,attribute,location in iterable:
# 			featureTableFile = os.path.join(featurePath,featureFile+'_{0}_{1}_{2}.txt'.format(edgetype,location,t))
# 			featureTable = np.loadtxt(featureTableFile,delimiter='\t', dtype='S1000')
# 			notBC = []
# 			BC = []

# 			atCol = np.where(featureTable[0,:]==attribute)[0][0]
# 			centcol = np.where(featureTable[0,:]==bcColName)[0][0]
# 			nodes = featureTable[np.where(featureTable[:,centcol]!=NOT_A_NODE_VALUE)][1:,0]
# 			bcvalues = featureTable[1:,centcol]
# 			bcvalues = bcvalues[np.where(bcvalues!=NOT_A_NODE_VALUE)]
# 			bcvalues= list([float(k) for k in bcvalues])
# 			bcvalues.sort(reverse=True)
# 			#cutoff = max(float(bcvalues[int(percentNodes*float(len(bcvalues)))-1]), bcMinValue)
# 			cutoff = max(float(percentNodes*max(bcvalues)), bcMinValue)
# 			for n in nodes:
# 				row = nm.findRow(n,featureTable)
# 				value = float(featureTable[row][centcol])
# 				if value == NOT_A_NODE_VALUE:
# 					print "UH oh"
# 				if value >= cutoff:
# 					BC.append(float(featureTable[row][atCol]))
# 				else:
# 					notBC.append(float(featureTable[row][atCol]))
# 			if BC:
# 				minvalue = min(min(notBC),min(BC))
# 				maxvalue = max(max(notBC),max(BC))
# 			else:
# 				minvalue = min(notBC)
# 				maxvalue = max(notBC)
# 			binwidth=(maxvalue-minvalue)/NUM_BINS
# 			if attribute==bcColName:
# 				if len(BC)>0:
# 					numSuperconnectors.append([t, location.split('_')[1], len(BC), round(min(BC),4),round(max(BC),4)])
# 				else:
# 					numSuperconnectors.append([t, location.split('_')[1], 0, None, None])

# 			notBC.extend(BC) # want to plot histogram of whole data and color the BC ones on top
# 			if binwidth == 0:
# 				bindwidth = 1
# 				ppl.hist(ax, notBC, color='black', grid='y')
# 				if len(BC)>0:
# 					ppl.hist(ax, BC, color='purple', grid='y')
# 			else:
# 				bins = np.arange(minvalue, maxvalue + binwidth, binwidth)
# 				ppl.hist(ax, notBC, color='black', grid='y', bins = bins, label = 'All OTUs')
# 				if len(BC)>0:
# 					ppl.hist(ax, BC, color='purple', grid='y', bins = bins, label = "Superconnectors")

# 			ax.locator_params(axis = "x",nbins=4)

# 			if attribute == abColName or attribute == bcColName:
# 				ax.set_yscale('symlog')
# 			if location in ytitles:
# 				ytitles.remove(location)
# 				print location
# 				ax.set_title(location.split('_')[1])
# 			if attribute in xtitles:
# 				xtitles.remove(attribute)
# 				ax.set_xlabel(attribute + ' histogram')

# 		# 	print ylim
# 		# 	ylim[attribute] = max(ax.get_ylim()[1],ylim[attribute])
# 		# 	xlims[attribute] = (min(minvalue,xlims[attribute][0]), max(maxvalue,xlims[attribute][1]))
			
# 		# for ax,attribute,location in iterable:
# 		# 	ax.set_ylim(0,ylim[attribute])
# 		# 	ax.set_xlim(xlims[attribute][0],xlims[attribute][1])


# 		#save number of superconnectors in a a table
# 		f = open(os.path.join(net_path,figurePath,"superconnectors_OTUs_{0}_{1}.txt".format(percentNodes,bcMinValue )),'w')
# 		f.write('\t'.join(['Treatment','Ecozone','Number of superconnectors','lowest BC value','highest BC value']))
# 		f.write('\n')
# 		for num in numSuperconnectors:
# 			f.write('\t'.join([str(x) for x in num]))
# 			f.write('\n')
# 		f.close()

# 		lgd = ppl.legend(bbox_to_anchor=(1.05, 1), loc=2)

# 		title = "Distribution of OTUs properties in all ecozones for treatment "+t
# 		figureTitle = fig.suptitle(title, horizontalalignment='center', fontsize=TITLE_FONT)

# 		fig.set_size_inches(BIG_FIG_WIDTH,BIG_FIG_HEIGHT)
# 		figureFile = os.path.join(net_path,figurePath,'scatter_BC_otus_'+t+'_'+'_'.join([l.split('_')[1] for l in locations])+'_'+edgetype+'_'+str(percentNodes)+'_'+str(bcMinValue)+'.png')
# 		fig.savefig(figureFile, dpi=DPI, bbox_extra_artists=(lgd,figureTitle),bbox_inches='tight')
# 		print "Saving the figure file: ", figureFile

# 	return None


def module_structure(net_path, networkNames, filePath, edgetype, inputFolder, inputFileEnd,featurePath, featureFile, factor):
	networks,treatments = get_network_fullnames(networkNames)
	print networks, treatments
	graphs = get_multiple_graphs(networks,net_path,edgetype, False, False)
	#otuTable = {}
	modules = {}
	number_modules = {}
	for n in networks:
		#otuTable[n] = np.loadtxt(os.path.join(inputFolder,n.replace('BAC_','')+inputFileEnd), dtype='S1000')
		mods = nm.get_module_graphs(graphs[n],factor=factor)
		modules[n] = mods
		number_modules[n] = len(mods)
	print number_modules
	print sum(number_modules.values())


	if treatments != []:
		table = np.zeros(shape=(len(MODULE_METRICS)+len(MODULE_OTU_METRICS)+2, sum(number_modules.values())+len(treatments)+1), dtype='S1000')
		i,j = 0,1 # i is row, j is column
		column = ['Zones','Treatments']
		column.extend([sm.__name__.replace('_',' ').capitalize() for sm in MODULE_METRICS])
		column.extend([om.__name__.replace('_',' ').capitalize() for om in MODULE_OTU_METRICS])
		table[:,0]=column
		for location,treatments in networkNames.iteritems():
			table[i,j]=location
			for t in treatments:
				i+=1
				number_mods = number_modules[location+'_'+t]
				mods = modules[location+'_'+t]
				print t
				table[i,j]= t
				print table
				for sm in MODULE_METRICS:
					print "For network for zone {0} treatment {1} calculating metric {2}".format(location,t,sm.__name__)
					i+=1
					values = []
					for mod in mods:
						values.append(sm(mod))
					table[i,j:j+number_mods]=values
				for om in MODULE_OTU_METRICS:
					print "For network for zone {0} treatment {1} calculating metric {2}".format(location,t,om.__name__)
					i+=1
					values = []
					featureTableFile = os.path.join(featurePath,featureFile+'_{0}_{1}_{2}.txt'.format(edgetype,location,t))
					featureTable = np.loadtxt(featureTableFile,delimiter='\t', dtype='S1000')
					for mod in mods:
						values.append(om(mod,featureTable))
					table[i,j:j+number_mods]=values
				j+= max(1,number_mods)
				i=0
	else:
		print 'Can only do for multiple treatments. FIX ME'

	np.savetxt(filePath, table, delimiter="\t", fmt='%s')
	return None




def network_structure(net_path, networkNames, filePath, edgetype, inputFolder, inputFileEnd,featurePath, featureFile):
	networks,treatments = get_network_fullnames(networkNames)
	print networks, treatments
	graphs = get_multiple_graphs(networks,net_path,edgetype, False, False)
	#sys.exit()
	otuTable = {}
	if INPUT_METRICS:
		for n in networks:
			otuTable[n] = np.loadtxt(os.path.join(inputFolder,n.replace('BAC_','')+inputFileEnd), dtype='S1000')

	if treatments != []:
		table = np.zeros(shape=(len(INPUT_METRICS)+len(STRUCTURE_METRICS)+len(OTU_METRICS)+2, len(networkNames)*len(treatments)+1), dtype='S1000')
		i,j = 0,1 # i is row, j is column
		column = ['Zones','Treatments']
		column.extend([sm.__name__.replace('_',' ').capitalize() for sm in INPUT_METRICS])
		column.extend([sm.__name__.replace('_',' ').capitalize() for sm in STRUCTURE_METRICS])
		column.extend([om.__name__.replace('_',' ').capitalize() for om in OTU_METRICS])
		table[:,0]=column
		for location,treatments in networkNames.iteritems():
			table[i,j]=location
			for t in treatments:
				i+=1
				table[i,j]=t
				for im in INPUT_METRICS:
					print "For input table from zone {0} treatment {1} measuring {2}".format(location,t,im.__name__)
					i+=1
					S = otuTable[location+'_'+t]
					table[i,j]=im(S)
				for sm in STRUCTURE_METRICS:
					print "For network for zone {0} treatment {1} calculating metric {2}".format(location,t,sm.__name__)
					i+=1
					G = graphs[location+'_'+t]
					table[i,j]=sm(G)
				for om in OTU_METRICS:
					print "For network for zone {0} treatment {1} calculating metric {2}".format(location,t,om.__name__)
					i+=1
					G = graphs[location+'_'+t]
					featureTableFile = os.path.join(featurePath,featureFile+'_{0}_{1}_{2}.txt'.format(edgetype,location,t))
					featureTable = np.loadtxt(featureTableFile,delimiter='\t', dtype='S1000')
					table[i,j]=om(G,featureTable)
				j+=1
				i=0
	else:
		print 'Can only do for multiple treatments. FIX ME'

	np.savetxt(filePath, table, delimiter="\t", fmt='%s')
	return None

def plot_multiple(net_path, networkNames, measures, plotby, fraction, figurePath, figureName, edgetype, add_random, add_scalefree, max_y):

	networks,treatments = get_network_fullnames(networkNames)
	graphs = get_multiple_graphs(networks,net_path,edgetype, add_random, add_scalefree, LCC=True)
	data = {}
	for netName,G in graphs.iteritems():
		print 'Running simulation on {0}.'.format(netName)
		rand_lc_sizes, rand_sc_sizes = random_attack(G, fraction)
		data[netName] = {'random':(rand_lc_sizes, rand_sc_sizes)}
		for m in measures:
			targ_lc_sizes, targ_sc_sizes = target_attack(G, m, fraction)
			data[netName][m.__name__] = (targ_lc_sizes, targ_sc_sizes)
	networkNamesPlot = networkNames.keys()
	title = 'Robustness simulation on LCC of networks {0} with {1} type of edges'.format(','.join([n.replace('BAC_','') for n in networkNamesPlot]), edgetype)
	if add_random:
		networkNamesPlot.extend([RAND_NAME+n for n in networkNames.keys()])
	if add_scalefree:
		networkNamesPlot.extend([SCALE_NAME+n for n in networkNames.keys()])
	if plotby == 'by_treatment':
		multi_plot_robustness_by_treatment(data, figurePath, figureName, networkNamesPlot, treatments, measures, fraction, net_path, title, max_y)
	elif plotby == 'by_measure':
		multi_plot_robustness_by_measure(data, figurePath, figureName, networkNamesPlot, treatments, measures, fraction, net_path, title, max_y)
	return None


def random_attack(G,fraction):
	'''Measure the size of the largest component of the graph
	as nodes are removed randomly'''
	lc_sizes = [] #relative size of big component
	sc_sizes = [] #avg size of smaller components
	startSize = len(nm.get_components(G)[0])
	lc_sizes.append(1)
	sc_sizes.append(1)
	H = G.copy()

	nodes= G.nodes()
	np.random.shuffle(nodes)

	removal=int(len(nodes)*fraction)-1 #can't remove last node, otherwise there is nothing to measure!

	for n in nodes[:removal]:
		H.remove_node(n)
		components = nm.get_components(H)
		lc_sizes.append(len(components[0])/float(startSize)) #measure the relative size change
		if len(components)>1:
			sc_sizes.append(np.mean([len(c) for c in components[1:]]))
		else:
			sc_sizes.append(1.0)
	return lc_sizes,sc_sizes


def target_attack(G, measure,fraction):
	'''Measure the size of the largest component of the graph
	as nodes are removed given the measure (degree or centrality 
	measure). The order of the nodes to be removed IS NOT updated
	after each removal.
	'''
	lc_sizes = [] #relative size of big component
	sc_sizes = [] #avg size of smaller components
	startSize = len(nm.get_components(G)[0])
	lc_sizes.append(1)
	sc_sizes.append(1)
	H = G.copy()

	values = [(n,v) for n,v in measure(G).iteritems()]

	values = sorted(values, key = lambda item: item[1], reverse = True)

	removal=int(len(values)*fraction)-1 #can't remove last node in simulation, otherwise there is nothing to measure!

	for n in zip(*values)[0][:removal]:
		H.remove_node(n)
		components = nm.get_components(H)
		lc_sizes.append(len(components[0])/float(startSize))  #measure the relative size change
		if len(components)>1:
			sc_sizes.append(np.mean([len(c) for c in components[1:]]))
		else:
			sc_sizes.append(1.0)
	return lc_sizes,sc_sizes

def plot_robustness(data,filename):
	'''plots the simulations'''

	# plotting locations in rows and centralities in columns
	fig, axes = plt.subplots(1)
	measures = ['random'] + [m.__name__ for m in measures]
	colors = {measure: ppl.colors.set2[i] for i,measure in enumerate(measures)}

	for measure in measures:
		values = data[measure]
		x = range(len(values))
		ppl.plot(axes, x, 
			values,
			label=str(measure))
			#color=[colors[measure] for measure in measures])

	ppl.legend(axes)  
	figureFile = os.path.join(net_path,filename)
	fig.savefig(figureFile)
	print "Saving the figure file: ", figureFile
	return None


def plot_individual(path,networkNames,fraction):
	networks,treatments = get_network_fullnames(networkNames)
	graphs = get_multiple_graphs(networks,path)

	for netName,G in graphs.iteritems():
		rand_lc_sizes, rand_sc_sizes = random_attack(G, fraction)
		data = {}
		data['random']= (rand_lc_sizes, rand_sc_sizes)
		for m in measures:
			targ_lc_sizes, targ_sc_sizes = target_attack(G, m, fraction)
			data[m.__name__] = (targ_lc_sizes, targ_sc_sizes)
		plot_robustness(data, netName)
	return None

def multi_plot_robustness_by_treatment(multidata,figurePath,figureFile,rowLabels,colLabels, measures, fraction, net_path, title, max_y):
	'''plots the simulations in a multiplot: each row is a location and each column is a treatment'''

	# plotting locations in rows and treatments in columns
	fig, axes = plt.subplots(len(rowLabels),len(colLabels))
	netNames = rowLabels

	measures = ['random'] + [m.__name__ for m in measures]

	colors = {measure: ppl.colors.set1[i] for i,measure in enumerate(measures)}
	#print netNames, measures, len(rowLabels),len(colLabels), len(axes), colLabels*len(rowLabels)

	iterable = []
	for i,r in enumerate(rowLabels):

		for j,c in enumerate(colLabels):
			if len(rowLabels)>1:
				iterable.append((axes[i][j],r,c))
			else:
				iterable.append((axes[j],r,c))

	net_label_done = []
	treatment_label_done = []
	x_axis_label_done = False

	#to add to legend
	ppl.plot([], [], color='black', linestyle='-', label='relative size of LCC')
	ppl.plot([], [], color='black', linestyle='--',linewidth=2, label='avg size of other CC')

	min_yvalue = 1
	max_yvalue = 1
	for ax,net,treatment in iterable:

		for measure in measures:
			lc_values = multidata[net+'_'+treatment][measure][0]
			sc_values = multidata[net+'_'+treatment][measure][1]
			min_yvalue = min(min_yvalue ,min(lc_values))
			max_yvalue = max(max_yvalue ,max(sc_values))
			x = [float(r)*fraction for r in range(len(lc_values))]
			ppl.plot(ax,
				x, 
				lc_values,
				marker='.',
				linestyle='-',
				label=str(measure.replace('_',' ')),
				color=colors[measure])
			ppl.plot(ax,
				x, 
				sc_values,
				color=colors[measure],
				linestyle='--',linewidth=2)
			
		if treatment not in treatment_label_done:
			ax.set_title(treatment)
			treatment_label_done.append(treatment)
		
		ax.set_title(treatment)

		if net not in net_label_done:
			ax.set_ylabel(net)
			net_label_done.append(net)

		if not x_axis_label_done:
			x_axis_label_done = True
			ax.set_xlabel('Number of removed nodes')

	for ax in axes:
		ax.set_autoscaley_on(False)
		if max_y and max_yvalue > max_y:
			max_yvalue = max_y
		ax.set_ylim([min_yvalue,max_yvalue])

	figureTitle = fig.suptitle(title,
         horizontalalignment='center',
         fontsize=TITLE_FONT) 

	lgd = ppl.legend(bbox_to_anchor=(1.05, 1), loc=2)

	figureFile = os.path.join(net_path,figurePath,figureFile)
	fig.set_size_inches(10*len(colLabels),7*len(rowLabels))
	fig.savefig(figureFile, dpi=DPI,  bbox_extra_artists=(lgd,figureTitle), bbox_inches='tight')
	print "Saving the figure file: ", figureFile
	return None





##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
##############################                     ###########################
#############                                                        #########
#################                                                 ############
####################                                    ######################
############################                       ###########################
#################################             ################################
###################################         ##################################
#######################################   ####################################
##############################################################################




def multi_plot_robustness_by_measure(multidata,figurePath,figureFile,rowLabels,treatments,measures,fraction, net_path, title, max_y):
	'''plots the simulations in a multiplot: each row is a location and each column is a centrality measure'''

	rowLabels.sort(reverse=True)
	# plotting locations in rows and centralities in columns
	measures = ['random'] + [m.__name__ for m in measures]
	fig, axes = plt.subplots(len(rowLabels),len(measures))
	netNames = rowLabels

	colors = {treatment: OM_COLORS[treatment] for i,treatment in enumerate(treatments)}
	#print netNames, measures, len(rowLabels),len(colLabels), len(axes), colLabels*len(rowLabels)

	robustnessTable = np.zeros(shape=(len(treatments)+1,len(measures)+1), dtype='S1000')
	robustnessTable[1:,0]=np.array(treatments)
	robustnessTable[0,1:]=np.array([m.replace('_',' ').capitalize() for m in measures])


	iterable = []
	for i,r in enumerate(rowLabels):
		for j,c in enumerate(measures):
			if len(rowLabels)>1:
				iterable.append((axes[i][j],r,c))
			else:
				iterable.append((axes[j],r,c))

	net_label_done = []
	measure_label_done = []
	x_axis_label_done = False

	#to add to legend
	if len(rowLabels)==1:
		ppl.plot([], [], color='black', marker= '.', linestyle='-', label='relative size of LCC')
		ppl.plot([], [], color='black', linestyle='--',linewidth=2, label='avg size of other CC')

	min_yvalue = 1
	max_yvalue = 1
	i,j=0,0
	for ax,net,measure in iterable:
		#indices for robustness factor table
		j+=1 
		i=0
		for t in treatments:
			lc_values = multidata[net+'_'+t][measure][0]
			sc_values = multidata[net+'_'+t][measure][1]
			min_yvalue = min(min_yvalue ,min(lc_values))
			max_yvalue = max(max_yvalue ,max(sc_values))
			x = [float(r)/len(lc_values)*fraction for r in range(len(lc_values))]
			if len(lc_values)>3:
				ppl.plot(ax,
					x, 
					lc_values,
					marker='.',
					linestyle='-',
					label=str(t),
					color=colors[t])
			else:
				ppl.plot(ax,
					[], 
					[],
					marker='.',
					linestyle='-',
					label=str(t),
					color=colors[t])
			if len(rowLabels)==1:
				ppl.plot(ax,
					x, 
					sc_values,
					color=colors[t],
					linestyle='--',linewidth=2)

			#Output robustness
			i+=1
			lc_values.pop(0) #remove '1' added for plotting purposes
			x = 0
			for k,l in enumerate(lc_values):
				if l<=0.5:
					x = k
					break

			#print lc_values
			#print x

			robustnessTable[i,j]=str(round(float(x)/len(lc_values),2))


		ax.locator_params(nbins=5)
		ax.tick_params(axis='both', which='major', labelsize=8)
		ax.tick_params(axis='both', which='minor', labelsize=8)
		if measure not in measure_label_done:
			ax.set_title(measure.replace('numpy','').replace('_',' ').capitalize())
			measure_label_done.append(measure)
		if net not in net_label_done:
			ax.set_ylabel("Relative size of LCC of {0} networks".format(net.split('_')[1]))
			net_label_done.append(net)

		if not x_axis_label_done:
			x_axis_label_done = True
			ax.set_xlabel('fraction of nodes removed')
		if len(rowLabels)>1:
			ax.set_ylim(0,1)


	#save robustnessTable in a table
	f = open(os.path.join(net_path,figurePath,"robustness_by_measure_{0}.txt".format(net)),'w')
	np.savetxt(f, robustnessTable, delimiter="\t", fmt='%s')
	f.close()

	if len(rowLabels)==1:
		for ax in axes:
			ax.set_autoscaley_on(False)
			if max_y and max_yvalue > max_y:
				max_yvalue = max_y
			ax.set_ylim([min_yvalue,max_yvalue])

	# figureTitle = fig.suptitle(title,
 #         horizontalalignment='center',
 #         fontsize=TITLE_FONT)

	lgd = ppl.legend(bbox_to_anchor=(1.05, 1), loc=2)

	figureFile = os.path.join(net_path,figurePath,figureFile)
	#fig.tight_layout()
	fig.set_size_inches(BIG_FIG_HEIGHT,BIG_FIG_WIDTH)
	#fig.savefig(figureFile, dpi=DPI, bbox_extra_artists=(lgd,figureTitle), bbox_inches='tight')
	fig.savefig(figureFile, dpi=DPI, bbox_extra_artists=([lgd]), bbox_inches='tight')
	print "Saving the figure file: ", figureFile

	return None

def make_js_files(netpath, newnetpath, ecozone, treatment, featurePath, featureFile, edgetype):
    '''make a network in js format'''

    featureTableFile = os.path.join(featurePath,featureFile+'_{0}_{1}_{2}.txt'.format(edgetype,ecozone,treatment))
    #featureTable = np.loadtxt(featureTableFile,delimiter='\t', dtype='S1000')
    edgeFile = os.path.join(netpath,ecozone+'_'+treatment+'_edges.txt')

    hive.get_nodes(featureTableFile,removeNA=NOT_A_NODE_VALUE)
    hive.get_edges(edgeFile)
    
    sources, targets, nodes, nodeProperties, edgeProperties = hive.sources, hive.targets, hive.nodes, hive.nodeProperties, hive.edgeProperties

    newnodeFile = os.path.join(newnetpath,"{0}_{1}_{2}_nodes.js".format(ecozone,treatment,edgetype))
    newedgeFile = os.path.join(newnetpath,"{0}_{1}_{2}_edges.js".format(ecozone,treatment,edgetype))

    f = open(newnodeFile, 'w')
    networkType = {'pos':'copresences ','neg':'mutual exclusions ','both':''}
    f.write('var SVGTitle = "Ecozone {0} treatment {1} {2}Hive Panel"\n'.format(ecozone.split('_')[1],treatment,networkType[edgetype]))
    f.write('var nodes = [\n')
    for i,node in enumerate(nodes):
		if 'Otu' in node:
			line = '    {name: \'' + str(node) +'\''
			for p,v in nodeProperties.iteritems():
				if p=="class": p = "Class"
				line  += ', '+p+': \'' + str(v[i]) + '\''
			line += '},\n'
			f.write(line)
		else: pass
    f.write('];')
    f.close()

    f = open(newedgeFile,'w')
    f.write('var links = [\n')
    for i,(s,t) in enumerate(zip(sources, targets)):
    	s = s.replace('OTU-','')
    	t = t.replace('OTU-','')
    	if 'Otu' in s and 'Otu' in t:
    		if edgetype != 'both' and (s not in nodes or t not in nodes):
    			continue
    		if edgetype == 'pos' and edgeProperties['interactionType'][i]!='[copresence, copresence]':
    			continue
    		elif edgetype == 'neg' and edgeProperties['interactionType'][i]!='[mutualExclusion, mutualExclusion]':
    			continue
    		else:
		    	line = '  {source: nodes['+str(nodes.index(s))+'], target: nodes['+str(nodes.index(t))+']'
		        for p,v in edgeProperties.iteritems():
		            line  += ', '+p+': \'' + str(v[i]) + '\''
	        	line += '},\n'
	        	f.write(line)
        else: pass
    f.write('];')
    f.close()
    return None















