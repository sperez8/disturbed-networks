'''
created  10/06/2014

by sperez

Contains functions used by hive class to measure things like network properties
'''

import sys
import os
import argparse

# #need a specific version of networkx for read_gexf to work
# #import pkg_resources
# #pkg_resources.require("networkx==1.7")

import networkx as nx
import string
import numpy as np
from math import pi

def import_gexf(gexfFile):
    #parse graphml file
    G = nx.read_gexf(gexfFile)
    return G

def import_graphml(graphmlFile):
    #parse graphml file
    G = nx.read_graphml(graphmlFile)
    return G

def import_graph(nodeFile, edgeFile, edgetype, filterNonOtus, filterEdges = True):
    '''make a networkx graph from a csv or tsv'''
    
    nodes, nodeProperties = get_nodes(nodeFile)
    sources, targets, edgeProperties = get_edges(edgeFile)
    
    G = make_graph(sources, targets, nodes, filterEdges)
    for i,n in enumerate(nodes):
        for p,v in nodeProperties.iteritems():
            G.node[n][p] = v[i]

    for i,e in enumerate(zip(sources, targets)):
        s,t = e[0],e[1]
        if filterEdges:
            if s not in nodes or t not in nodes:
                continue
        if edgetype == 'pos' and 'mutualExclusion' in edgeProperties['interactionType'][i]:
            G.remove_edge(s,t)
        elif edgetype == 'neg' and 'copresence' in edgeProperties['interactionType'][i]:
            G.remove_edge(s,t)
        elif (filterEdges and s in nodes and t in nodes) or (not filterEdges):
            for p,v in edgeProperties.iteritems():
                G[s][t][p] = v[i]

    if filterNonOtus:
        for node in G.nodes():
            if 'Otu' not in node:
                G.remove_node(node)

    for n in G.nodes():
        if G.degree(n) == 0:
            G.remove_node(n)

    return G

def convert_gexf(gexfFile):
    G = import_gexf(gexfFile)
    fileName = gexfFile.split('.gexf')[0]
    convert_graph(G,fileName)

def convert_graphml(graphmlFile):
    G = import_graphml(graphmlFile)
    fileName = graphmlFile.split('.graphml')[0]
    convert_graph(G,fileName)

def convert_graph(G,fileName):
    sources, targets, edgeProperties = zip(*G.edges(data=True))
    nodeFile = fileName+'_nodes.csv'
    edgeFile = fileName+'_edges.csv'
    nf = open(nodeFile,'w')
    keys = []
    for node, nodeProperties in G.nodes(data=True):
        new_keys = nodeProperties.keys()
        if new_keys != keys:
            keys.extend(new_keys)
    
    keys = set(keys)
    
    #write header
    if keys:
        nf.write('Node'+','+','.join(keys))
    else:
        nf.write('Node')        
    
    for node, nodeProperties in G.nodes(data=True):
        row = []
        row.append(node)
        for k in keys:
            if k in nodeProperties.keys():
                row.append(str(nodeProperties[k]).replace(',', ';'))
            else:
                row.append('None')
        nf.write('\n' + ','.join([str(r) for r in row]))
    
    nf.close()
    
    ef = open(edgeFile,'w')
    keys = []
    for source,target, edgeProperties in G.edges(data=True):
        new_keys = edgeProperties.keys()
        if new_keys != keys:
            keys.extend(new_keys)
    
    keys = set(keys)
    
    if keys:
        ef.write('source' + ',' + 'target' + ',' +','.join(keys))
    else: 
        ef.write('source' + ',' + 'target')

    for source, target, edgeProperties in G.edges(data=True):
        row = []
        row.append(source)
        row.append(target)
        for k in keys:
            if k in edgeProperties.keys():
                row.append(str(edgeProperties[k]))
            else:
                row.append('None')
        ef.write('\n' + ','.join([str(r) for r in row]))
    
    print "writing nodefile", nodeFile
    print "writing edgefile", edgeFile
    return None

def make_graph(sources, targets, nodes, filterEdges= True):
    '''Makes a graph using the networkx package Graph instance'''
    G = nx.Graph()
    G.add_edges_from(zipper(sources,targets))
    if filterEdges:
        for n in G.nodes():
            if n not in nodes:
                G.remove_node(n)
    return G


def get_nodes(inputFile,removeNA=None):
    '''gets nodes and their properties from csv file'''
    
    delimiter = get_delimiter(inputFile)

    data = np.genfromtxt(inputFile, delimiter=delimiter, dtype='str', filling_values = 'None')
    
    #get properties and format as strings
    properties = data[0,1:]
    properties = format_properties(properties)
    
    if removeNA:
        colName = nx.betweenness_centrality.__name__.replace('_',' ').capitalize()
        col = np.where(data[0,:]==colName)[0][0]
        #remove first row with column names
        data = data[1:,]
        data = data[np.where(data[:,col]!=removeNA)]
    else:
        #remove first row with column names
        data = data[1:,]

    #get all the node data
    nodes = list(data[:,0])
    
    #take note of number of nodes
    totalNodes = len(nodes)
    
    #transform node properties into the numerical types if possible
    nodeProperties = {}

    for i, column in enumerate(data[:,1:].T):
        values = convert_type(list(column))
        nodeProperties[properties[i]] = values
    nodeProperties = nodeProperties

    return nodes, nodeProperties


def get_edges(inputFile):
    '''gets edges and their properties from csv file'''
    
    delimiter = get_delimiter(inputFile)
    data = np.genfromtxt(inputFile, delimiter=delimiter, dtype='str', filling_values = 'None')
    
    #get properties and format as strings
    properties = data[0,2:]
    properties = format_properties(properties)
    
    #remove first row with column names
    data = data[1:,]
    
    #get all the edge data
    sources = list(data[:,0])        
    targets = list(data[:,1])
    
    #take note of number of edges:
    totalEdges = len(sources)

    #transform edge properties into the numerical types if possible
    edgeProperties = {}

    for i, column in enumerate(data[:,2:].T):
        values = convert_type(list(column))
        edgeProperties[properties[i]] = values
    edgeProperties = edgeProperties
    
    #store the name of the edge properties
    edgePropertyList = edgeProperties.keys()

    return sources, targets, edgeProperties

def convert_type(data):
    def num(s):
        '''convert list of strings to corresponding int or float type'''
        try:
            return int(d)
        except ValueError:
            return float(d)
    
    try:
        convertedData = [num(d) for d in data]
        return convertedData
    except ValueError:
        return data

def get_delimiter(inputFile):
    '''detect if input file is a tab or comma delimited file
        and return delimiter.'''
    
    ext = os.path.splitext(os.path.basename(inputFile))[1]
    
    if 'tab' in ext or 'tsv' in ext:
        return '\t'
    elif 'csv' in ext:
        return ','
    elif 'txt' in ext:
        #detects delimiter by counting the number of tabs and commas in the first line
        f = open(inputFile, 'r')
        first = f.read()
        if first.count(',') > first.count('\t'):
            return ','
        elif first.count(',') < first.count('\t'):
            return '\t'
        else:
            print "Couldn't detect a valid file extension: ", inputFile
            return ','
    else:
        print "Couldn't detect a valid file extension: ", inputFile
        return ','


def format_properties(properties, debug = False):
    '''takes a list of property names and removes all punctuation and numbers'''
    
    numbers = {1:'one', 2:'two', 3:'three', 4:'four', 5:'five', 6:'six', 7:'seven', 8:'eight', 9:'nine', 10:'ten'}
    
    def convert_word(word):
        '''remove punctuation and numbers from a word'''
        w = word
        word = ''.join(word.split()) #removes all whitespace (tabs, newlines, spaces...)
        for c in string.punctuation + string.digits:
            word = word.replace(c,'')
        if w != word:
            if debug:
                print "The property \'{0}\' contains spaces, punctuation or digits and has been renamed '{1}'".format(w,word)
        return word
         
    newProperties = []
    i = 1
    for prop in properties:
        newProp = convert_word(prop)
        if not newProp:
            #if property isn't named, we give it one
            newProperties.append('unNamedProperty' + numbers[i] + '')
            i += 1
        elif newProp in newProperties:
            newProperties.append(newProp + 'second')
        else:
            newProperties.append(newProp)
            
    return newProperties

def zipper(*args):
    '''a revamped version of zip() method that checks that lists
    to be zipped are the same length'''
    for i,item in enumerate(args):
        if len(item) != len(args[0]):
            raise ValueError('The lists to be zipped aren\'t the same length.')
    
    return zip(*args)

def main(*argv):
    '''handles user input and runs plsa'''
    parser = argparse.ArgumentParser(description='This scripts converts networks to txt node and edge files')
    parser.add_argument('-input', help='Location of network file')
    parser.add_argument('-format', help='Input format of network')
    args = parser.parse_args()

    if args.format=='graphml':
        print "Converting graphml input file: ", args.input
        convert_graphml(args.input)
    if args.format=='gexf':
        convert_gexf(args.input)

if __name__ == "__main__":
    main(*sys.argv[1:])

# file = "C:\\Users\\Sarah\\Dropbox\\1-Hive panels\\Diseasome\\diseasome.gexf"
# f = open(file,'r')
# print f.readlines()
# convert_gexf(file)





