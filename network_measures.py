'''
created  10/06/2014

by sperez

'''

#library imports
import sys
import os
import argparse
import numpy as np
from math import pi
import copy
import scipy.stats

_cur_dir = os.path.dirname(os.path.realpath(__file__))
_root_dir = os.path.dirname(_cur_dir)
sys.path.insert(0, _root_dir)

import networkx as nx
from make_network import *

DECIMALS = 3 #for rounding measures
FACTOR = 1.5

SOILHORIZON_FEAT_NAME = 'SoilHorizon avg'

def number_of_nodes(G):
	return G.number_of_nodes()

def number_of_edges(G):
	return G.number_of_edges()

def number_of_nodes_of_largest_connected_component(G):
    return len(get_components(G)[0])

def number_of_edges_of_largest_connected_component(G):
    return get_LCC(G).number_of_edges()

def number_of_components(G):
    return nx.number_connected_components(G)

def size_of_big_components(G):
    cc = get_components(G)
    sizes = [str(len(c)) for c in cc if len(c) > 3]
    return ','.join(sizes)

def in_largest_connected_component(G):
    LCC = get_components(G)[0]
    members = {n:(1 if n in LCC else 0) for n in G.nodes()}
    return members

def node_degrees(G):
    return G.degree()

def average_degree(G):
    return round(np.mean(G.degree().values()), DECIMALS)

def connectance(G):
    return round(nx.density(G), DECIMALS)

def global_clustering_coefficient(G):
    return round(nx.average_clustering(G), DECIMALS)

def fraction_of_possible_triangles(G):
    return round(nx.transitivity(G), DECIMALS)

def size_of_largest_clique(G):
    return nx.graph_clique_number(G)

def degree_assortativity(G):
    return round(nx.degree_assortativity_coefficient(G), DECIMALS)

def diameter_of_largest_connected_component(G):
    H = get_LCC(G)
    return nx.diameter(H)

def average_path_on_largest_connected_component(G):
    H = get_LCC(G)
    return round(nx.average_shortest_path_length(H), DECIMALS)

def correlation_of_degree_and_betweenness_centrality(G):
    bc = nx.betweenness_centrality(G)
    d = nx.degree(G)
    bcn = []
    dn = []
    for n in G.nodes():
        bcn.append(bc[n])
        dn.append(d[n])

    r = scipy.stats.spearmanr(dn, bcn)
    return format_correlation(r[0],r[1])

def format_correlation(avg,std):
    return "{0} ({1})".format(round(avg,DECIMALS), round(std,5))


def get_components(G):
    '''gets connected components, sorts by size and returns a list of lists'''
    return sorted(nx.connected_components(G), key = len, reverse=True)

def get_LCC(G):
    '''gets connected subgraphs and returns LCC as a networkx graph'''
    LCC = None
    Ntemp = 0
    for graph in nx.connected_component_subgraphs(G):
        N = graph.number_of_nodes()
        if N > Ntemp:
            LCC = graph
            Ntemp = N
    return LCC


### Ecological measures


def remove_headers(S):
    return S[1:-1,1:-1].astype(np.float)

def normalize(S):
    col_sums = S.sum(axis=0)
    nS = S / col_sums[np.newaxis,:]
    return nS


def richness(S):
    S = remove_headers(S)
    return S.shape[0]

def shannon_diversity(S):
    D = 0
    S = normalize(remove_headers(S))
    D = -sum(np.mean(row) * np.log(np.mean(row)) for row in S)    
    return str(round(D,DECIMALS))


########## Measures using an OTU table and features

def correlation_of_edge_depth(G,featureTable):
    feature = SOILHORIZON_FEAT_NAME
    return compute_feature_correlation(G,feature,featureTable)


def correlation_of_degree_and_depth(G,featureTable):
    feature = SOILHORIZON_FEAT_NAME
    return compute_feature_degree_correlation(G,feature,featureTable)
    
def compute_feature_degree_correlation(G,feature,featureTable):
    col = np.where(featureTable[0,:]==feature)[0][0]
    d = nx.degree(G)
    degrees = []
    featureValues = []
    for n in d.keys():
        row = findRow(n,featureTable)
        if row:
            degrees.append(d[n])
            featureValues.append(featureTable[row][col])
        else: continue
    r = scipy.stats.spearmanr(degrees, featureValues)
    return format_correlation(r[0],r[1])


def compute_feature_correlation(G,feature,featureTable):
    col = np.where(featureTable[0,:]==feature)[0][0]
    iF = []
    jF = []
    for (i,j) in G.edges():
        irow = findRow(i,featureTable)
        jrow = findRow(j,featureTable)
        if irow and jrow:
            iF.append(featureTable[irow][col])
            jF.append(featureTable[jrow][col])
        else:
            continue
    r = scipy.stats.spearmanr(iF,jF)
    return format_correlation(r[0],r[1])

###FIX MEEE
def average_depth(G,feature,featureTable):
    col = np.where(featureTable[0,:]==feature)[0][0]
    iF = []
    jF = []
    for (i,j) in G.edges():
        irow = findRow(i,featureTable)
        jrow = findRow(j,featureTable)
        if irow and jrow:
            iF.append(featureTable[irow][col])
            jF.append(featureTable[jrow][col])
        else:
            continue
    r = scipy.stats.spearmanr(iF,jF)
    return format_correlation(r[0],r[1])

def findRow(otu,table):
    if 'OTU' in otu:
        row = np.where(table==otu.replace('OTU-',''))[0][0]
    elif 'Otu' in otu:
        row = np.where(table==otu)[0][0]
    else:
        row = None
        print "WARNING: Didn't find the otu: ", otu
        #sys.exit()
    return row





NODES = os.path.join(_root_dir, 'tests', 'test_nodes_friends.txt')
EDGES = os.path.join(_root_dir, 'tests', 'test_edges_friends.txt')


def main(*argv):
    '''handles user input and runs plsa'''
    parser = argparse.ArgumentParser(description='This script creates a networkx graph.')
    parser.add_argument('-n', help='The node file', default = NODES)
    parser.add_argument('-e', help='The edge file', default = EDGES)
    args = parser.parse_args()
    
    if (args.n == '' and args.e != '') or (args.n != '' and args.e == ''):
        print "\n***You must specify both a node and an edge file if specifying either.***\n"
        parser.print_help()
        sys.exit()
        
    nodeFile = args.n
    edgeFile = args.e

    G = import_graph(nodeFile,edgeFile)

    print "\nMade the networkx graph."

    return G


def compute_modularity_horizon(G,featureTable,modules = None):
    feature = SOILHORIZON_FEAT_NAME
    return compute_modularity_feature(G,feature,featureTable,modules = modules)

def compute_modularity_feature(G,feature,featureTable,factor=FACTOR,modules = None):
    col = np.where(featureTable[0,:]==feature)[0][0]
    modularity = node_modularity(G, modules = modules)
    # H = nx.Graph()
    module_features = {m:[] for m in set(modularity.values())}

    for node,mod in modularity.iteritems():
            row = findRow(node,featureTable)
            if row:
                module_features[mod].append(float(featureTable[row][col]))
            else: continue

    feature_values = {}
    for m,values in module_features.iteritems():
       avg = np.average(values)
       std = np.std(values)
       feature_values[m] = format_correlation(avg, std)
    return ';'.join([str(k)+':'+str(v) for k,v in feature_values.iteritems()])


def module_sizes(G,factor=FACTOR,modules = []):
    if not modules:
        modules = get_modules(G)
    if modules:
        modules.sort(key=lambda m: len(m),reverse=True) #order by size
        return [str(len(m)) for m in modules]
    else:
        return 'None'

def get_module_graphs(G,factor=FACTOR,modules = []):
    if not modules:
        modules = get_modules(G)
    if modules:
        modules.sort(key=lambda m: len(m),reverse=True) #order by size
        graphs = []
        for m in modules:
            H = nx.Graph()
            H.add_nodes_from(m)
            for s,t in G.edges():
                if s in m and t in m:
                    H.add_edge(s,t)
            graphs.append(H)
        return graphs
    else:
        return []

def module_connectance(G,factor=FACTOR,modules = []):
    if not modules:
        modules = get_modules(G)
    if modules:
        modules.sort(key=lambda m: len(m),reverse=True) #order by size
        connectances = []
        for i,m in enumerate(modules):
            n = len(m)
        return [str(e) for e in connectances]
    else:
        return 'None'

def node_modularity(G,factor=FACTOR, modules = None):
    '''gets modules using FAG-EC algorithm and returns a dictionary
    where keys are nodes and value is e module it belongs to where
    0 = no module
    1 = largest module
    i = ith largets module
    '''
    modularity = {}
    if not modules:
        modules = get_modules(G)
    modules.sort(key=lambda m: len(m),reverse=True) #order by size
    print "Found {0} modules with sizes: {1}".format(len(modules),','.join([str(len(m)) for m in modules]))
    for n in G.nodes():
        m = findSubgraph(modules,n)
        if m != None:
            m+=1  #start module index at 1
        else:
            m = 0 #for not in a module
        modularity[n]=m
    return modularity

def findSubgraph(subgraphs,n):
    for i,m in enumerate(subgraphs):
        if n in m:
            return i

def testModule(G,factor,subgraphNodes):
    isModule = False

    kin = 0
    kout = 0
    for s,t in G.edges():
        if s in subgraphNodes and t in subgraphNodes:
            kin+=1
        elif s in subgraphNodes or t in subgraphNodes:
            kout+=1

    #print 'degrees', kin,kout, G.degree(subgraphNodes)
    if kin>kout*factor:
        isModule = True
        #print "MODULE"

    return isModule

def edge_clustering(G):
    clusteringcoeffs = {}
    for s,t in G.edges():
        c = 0
        ns = set(nx.all_neighbors(G,s))
        nt = set(nx.all_neighbors(G,t))
        commons = ns.intersection(nt)
        ds = G.degree(s)
        dt = G.degree(t)
        c = (len(commons)+1)/float(min([ds,dt]))
        clusteringcoeffs[(s,t)] = c
    return clusteringcoeffs

def get_modules(G,factor=FACTOR):
    '''modularity algorithm from FAG-EC'''
    #initialize nodes as singleton clusters
    subgraphs = [[n] for n in G.nodes()]
    nonmergeable = []
    #get edge betweenness values and sort them by that value
    weights = edge_clustering(G)
    Sq = [(e,cc) for e,cc in weights.iteritems()]
    Sq.sort(key=lambda tup: tup[1],reverse=True)

    while len(Sq)>0:
        edge,cc = Sq.pop(0) #get mergeable edge with highest clustering coefficient
        s,t = edge[0],edge[1]
        mods = findSubgraph(subgraphs,s) #get index of the subgraph where node belongs
        modt = findSubgraph(subgraphs,t)

        if mods==modt: #already in the same subgraph
            continue
        #check if mergeable, ie. both not in nonmergeable modules
        elif findSubgraph(nonmergeable,s)==None and findSubgraph(nonmergeable,t)==None:
            ms = testModule(G,factor,subgraphs[mods])
            mt = testModule(G,factor,subgraphs[modt])
            if ms and mt: #if both modules, then non mergeable
                newmod = subgraphs[mods]
                nonmergeable.append(newmod)
                newmod = subgraphs[modt]
                nonmergeable.append(newmod)
            else: #merge
                newsubg = subgraphs.pop(mods)
                modt = findSubgraph(subgraphs,t) #need to do it again after poping
                s2 = subgraphs.pop(modt)
                newsubg.extend(s2)
                subgraphs.append(newsubg) #merging

        #THE FOLLOWING lines are consistent with the original algorithm,
        # which outputs all subgraphs, some modules some not.
        # elif findSubgraph(nonmergeable,s)==None:
        #     newmod = subgraphs[mods]
        #     nonmergeable.append(newmod)
        #     #set as non mergeable
        # elif findSubgraph(nonmergeable,t)==None:
        #     #set as non mergeable
        #     newmod = subgraphs[modt]
        #     nonmergeable.append(newmod)
        else: #both are modules, or one is a module so we continue
            continue

    return nonmergeable #return modules only, ie. not all nodes are returned.

if __name__ == "__main__":
    '''testing purposes'''
    #G = main(*sys.argv[1:])

    #to test modularity
    # G=nx.karate_club_graph()
    # m = node_modularity(G)
    # for i,n in m.iteritems():
    #     print i,n

    '''Use  commands like:
    nx.attribute_mixing_matrix(G, 'Gender')
    and
    nx.attribute_assortativity_coefficient(G, 'Gender')
    to find the assortativity between node attributes'''


