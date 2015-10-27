'''
Created on 27/03/2015

author: sperez8
'''

import sys, os
import argparse
import numpy as np
from network_simulation import *

#What to plot
import platform
if platform.system() == 'Windows':
	PATH = '\Users\Sarah\Desktop\LTSPnetworks'
else:
	PATH = '/Users/sperez/Desktop/LTSPnetworks'

print PATH
FOLDER = 'by_treatment'
WHOLE_FOLDER = 'by_zone'

INPUT_FOLDER = 'input'
INPUT_FILE_END = '_BAC-filtered-lineages_final.txt'

INDVAL_FOLDER = 'indtables'
INDVAL_FILE_END = '_indvals_combo_om_horizon.txt'

FEATURE_FILE = 'feature_and_node_measures_table'
FEATURES = ['SoilHorizon']
BC_FEATURES = ['Betweenness centrality','SoilHorizon avg','SoilHorizon std','Abundance']
MEASURES = [nx.betweenness_centrality, 
			nx.degree_centrality,
			nx.closeness_centrality,
			#nx.eigenvector_centrality_numpy,
			 ]
PERCENT_BC_NODES = 0.1
BC_MIN_VALUE = 0.005

TAX_LEVEL = 'phylum'

#TREATMENTS = ['OM3','OM2'] #use when testing
TREATMENTS = ['OM0','OM1','OM2','OM3']
TREATMENTS.sort()

PROP_TO_REMOVE = 1 #only removing this percent of nodes
MAX_Y_AXIS = 5.5
DEGREE_SEQUENCE = False

FACTOR = 2

def main(*argv):
	'''handles user input and runs plsa'''
	parser = argparse.ArgumentParser(description='This scripts analyzes co-occurrence networks')
	parser.add_argument('-path', help='Path where the networks are', default = PATH)
	parser.add_argument('-folder', help='Folder in path where the networks are', default = FOLDER)
	parser.add_argument('-factors', nargs='*', help='Different types of each network', default = TREATMENTS)
	parser.add_argument('-networks', nargs='*', help='Which network to use: SBS, IDF, etc.')
	parser.add_argument('-simulate', help='Simulates node removal', action = 'store_true')
	parser.add_argument('-calculate', help='Calculates networks properties', action = 'store_true')
	parser.add_argument('-modules', help='Calculates module properties', action = 'store_true')
	parser.add_argument('-factor', help='Sets in-out degree factor for modularity', default = FACTOR)
	parser.add_argument('-distribution', help='Plots degree distribution', action = 'store_true')
	parser.add_argument('-assess', help='Assess ecological properties', action = 'store_true')
	parser.add_argument('-maketable', help='Make OTU table with eclogical measures', action = 'store_true')
	parser.add_argument('-edgetype', help='Specify which types edges to use', default = 'both')
	#arguments used when running simulations
	parser.add_argument('-fraction', help='Fraction of nodes to remove', default = PROP_TO_REMOVE)
	parser.add_argument('-addrandom', help='Runs simulation on random network of same size', action = 'store_true')
	parser.add_argument('-addscalefree', help='Runs simulation on scale network of same size', action = 'store_true')
	parser.add_argument('-treatment', help='Makes a plot for each treatment', action = 'store_true')
	parser.add_argument('-measure', help='Makes a plot for each centrality measure', action = 'store_true')
	parser.add_argument('-showcomponents', help='Average size of large component fragments to show', default = MAX_Y_AXIS)
	parser.add_argument('-wholenetwork', help='Makes a plot for whole network, not per treatments', action = 'store_true')
	#arguments and plots for central OTUs
	parser.add_argument('-boxplot', help='Makes a boxplot per taxonomic level of otu centrality', action = 'store_true')
	parser.add_argument('-level', help='Selects taxonomic level at which to make the boxplot', default = TAX_LEVEL)
	parser.add_argument('-bcplot', help='Makes a boxplot per treatment high BC otu features', action = 'store_true')
	parser.add_argument('-percentnodes', help='Select the proportion of high bc nodes to plot', default = PERCENT_BC_NODES)
	parser.add_argument('-bcmin', help='Cutoff for BC values', default = BC_MIN_VALUE)
	parser.add_argument('-vennplot', help='Makes a venn diagram of high BC otu per ecozone', action = 'store_true')
	parser.add_argument('-makejs', help='Makes node and edges files in .js', action = 'store_true')
	parser.add_argument('-scatterplot', help='Makes scatter plot of BC nodes', action = 'store_true')
	parser.add_argument('-plotcentralities', help='Makes scatter plot different centralities', action = 'store_true')
	parser.add_argument('-taxarep', help='Calculates representation of phylum in central OTUs', action = 'store_true')
	args = parser.parse_args()
	
	#check that one of the options is true
	choices = [args.simulate,args.distribution,args.calculate,
				args.assess,args.maketable,args.boxplot,
				args.modules,args.bcplot,args.vennplot,
				args.makejs, args.scatterplot, args.taxarep,
				args.plotcentralities]
	if sum([1 for c in choices if c])>1 or sum([1 for c in choices if c])==0:
		print "\n***You must specify one of the three options to calculate porperties of, run simulations on or plot networks.***\n"
		parser.print_help()
		sys.exit()	

	if args.edgetype not in ['both','pos','neg']:
		print "\n***You must specify what edges you want to use to build the network: both, pos or neg.***\n"
		parser.print_help()
		sys.exit()

	factors = args.factors

	edgetype = args.edgetype
	net_path = os.path.join(args.path,args.folder)
	print net_path
	if args.folder == 'by_zone':
		networks = {('BAC_'+n if 'BAC_' not in n else n):[] for n in args.networks}
	else:
		networks = {('BAC_'+n if 'BAC_' not in n else n):factors for n in args.networks}

	factor = float(args.factor)

	folderNewPath = os.path.join(args.path,'panels','data')
	figurePath = os.path.join(args.path,'plots')
	samplesFile = os.path.join(args.path, 'Bacterialtags_info_edited.txt')
	featurePath = os.path.join(args.path,'tables')


	###depending on option specified, choose different things
	if args.calculate:
		print "\nCalculating structural properties on "+edgetype+" type of edges of networks:"
		print ", ".join(networks), '\n'
		filePath = os.path.join(figurePath,'table_of_measures_'+'_'.join(args.networks)+'_'+edgetype+'.txt')
		print filePath
		network_structure(net_path,networks,filePath,edgetype, os.path.join(args.path,INPUT_FOLDER),INPUT_FILE_END, featurePath, FEATURE_FILE)

	elif args.modules:
		print "\nCalculating structural properties on "+edgetype+" type of edges of modules in networks:"
		print ", ".join(networks), '\n'
		filePath = os.path.join(figurePath,'table_of_module_measures_'+'_'.join(args.networks)+'_'+edgetype+'_'+str(factor)+'.txt')
		print filePath
		module_structure(net_path,networks,filePath,edgetype, os.path.join(args.path,INPUT_FOLDER),INPUT_FILE_END, featurePath, FEATURE_FILE, factor)

	elif args.assess:
		print "\nCalculating ecological metrics of sample collection for the following networks:"
		print ", ".join(networks), '\n'
		filePath = os.path.join(figurePath,'ecological_measures_'+'_'.join(args.networks)+'.txt')
		make_ecological_table(net_path,networks,filePath,edgetype,os.path.join(args.path,INPUT_FOLDER),INPUT_FILE_END)

	elif args.maketable:
		print "\nMaking OTU table with ecological metrics for the following networks:"
		print ", ".join(networks), '\n'
		make_OTU_feature_table(net_path, networks, os.path.join(args.path,INPUT_FOLDER),INPUT_FILE_END,os.path.join(args.path,INDVAL_FOLDER),INDVAL_FILE_END, samplesFile, FEATURES, featurePath, FEATURE_FILE,edgetype,factor)

	elif args.distribution:
		if DEGREE_SEQUENCE:
			figureName = 'plot_distribution_'+'_'.join(args.networks)+'_'+edgetype+'_sequence'+'.png'
		else:
			figureName = 'plot_distribution_'+'_'.join(args.networks)+'_'+edgetype+'.png'
		figurePath = os.path.join(figurePath,figureName)
		print "\nPlotting the degree distribution on "+edgetype+" type of edges for networks ", ','.join(args.networks)
		plot_degree_distribution_per_treatment(net_path, networks, figurePath, DEGREE_SEQUENCE, edgetype)

	elif args.boxplot:
		edgetype = 'pos'
		percentNodes = float(args.percentnodes)
		bcMinValue = float(args.bcmin)
		level = args.level
		if level not in TAXONOMY:
			print level, "is not a taxonomic level"
			sys.exit()
		print "\nPlotting "+level+" centrality for OTUs in the following networks with "+edgetype+" type of edges:"
		print ", ".join(networks), '\n'
		centrality_plot(net_path,networks,figurePath,featurePath,FEATURE_FILE,level,percentNodes,bcMinValue)

	elif args.vennplot:
		edgetype = 'pos'
		percentNodes = float(args.percentnodes)
		bcMinValue = float(args.bcmin)
		level = args.level
		if level not in TAXONOMY:
			print level, "is not a taxonomic level"
			sys.exit()
		print "\nPlotting "+level+" venn diagram of central OTUs per ecozone with "+edgetype+" type of edges:"
		print ", ".join(networks), '\n'
		if level == "species":
			plot_venn_otus_diagram(net_path,networks,figurePath,featurePath,FEATURE_FILE,percentNodes,bcMinValue)
		else:	
			plot_venn_diagram(net_path,networks,figurePath,featurePath,FEATURE_FILE,level,percentNodes,bcMinValue)

	elif args.taxarep:
		edgetype = 'pos'
		percentNodes = float(args.percentnodes)
		bcMinValue = float(args.bcmin)
		level = args.level
		if level not in TAXONOMY:
			print level, "is not a taxonomic level"
			sys.exit()
		print "\nCalculating representation of "+level+" in central OTUs per ecozone"
		print ", ".join(networks), '\n'
		calculate_taxonomic_representation(net_path,networks,figurePath,featurePath,FEATURE_FILE,level,percentNodes,bcMinValue)

	elif args.scatterplot:
		edgetype = 'pos'
		percentNodes = float(args.percentnodes)
		bcMinValue = float(args.bcmin)
		print "\nPlotting scatterplot to compare central OTUs per ecozone with "+edgetype+" type of edges:"
		print ", ".join(networks), '\n'
		plot_scatter_bc(net_path,networks,figurePath,featurePath,FEATURE_FILE,percentNodes,bcMinValue)

	elif args.plotcentralities:
		edgetype = 'pos'
		percentNodes = 0.5 #float(args.percentnodes)
		print "\nPlotting centrality measures per ecozone with "+edgetype+" type of edges:"
		print ", ".join(networks), '\n'
		measures = MEASURES
		plot_diff_centralities(net_path,networks,figurePath,featurePath,FEATURE_FILE,percentNodes,measures)

	elif args.makejs:
		print "\nMaking node and edge file in .js format for following ecozones with "+edgetype+" type of edges:"
		print ", ".join(networks), '\n'
		for n in networks:
			for t in TREATMENTS:
				make_js_files(net_path,folderNewPath,n,t,featurePath,FEATURE_FILE,edgetype)

	elif args.bcplot:
		edgetype = 'pos'
		percentNodes = float(args.percentnodes)
		bcMinValue = float(args.bcmin)
		print "\nPlotting different features of high betweenness centrality OTUs in the following networks with "+edgetype+" type of edges:"
		print ", ".join(networks), '\n'
		keystone_quantitative_feature_plot(net_path,networks,figurePath,featurePath, FEATURE_FILE, BC_FEATURES, percentNodes,bcMinValue)


	elif args.simulate:
		if not args.treatment and not args.measure:
			print "\n***You must specify to plot by treatment or by measure.***\n"
			parser.print_help()
			sys.exit()
		if args.treatment:
			plot_by = 'by_treatment'
		elif args.measure:
			plot_by = 'by_measure'
		if args.addscalefree:
			add_scalefree = True
		else: 
			add_scalefree = False
		if args.addrandom:
			add_random = True
		else:
			add_random = False
		max_y = args.showcomponents
		if max_y.lower() == 'none':
			max_y == None
		else:
			max_y = float(max_y)

		fraction = float(args.fraction)
		if len(networks)>1:
			figureName = 'plot_'+'_'.join(args.networks)+'_'+edgetype+'_'+ plot_by+'.png'
		else:
			figureName = 'plot_'+'_'.join(args.networks)+'_'+edgetype+'_'+ plot_by +'_prop='+str(fraction)+'_maxy='+str(max_y)+'.png'
		measures = MEASURES

		print "\nSimulating and plotting the robustness on "+edgetype+" type of edges of networks:"
		print ", ".join(networks)
		print "and plotting "+str(fraction)+" fraction of nodes "+plot_by+" and with following measures:"
		print ", ".join([m.__name__ for m in measures])
		print "\n"
		plot_multiple(net_path, networks, measures, plot_by, fraction, figurePath, figureName, edgetype, add_random, add_scalefree, max_y)
	
if __name__ == "__main__":
	main(*sys.argv[1:])




