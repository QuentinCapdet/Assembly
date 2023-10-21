#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Quentin"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Quentin"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Quentin"
__email__ = "quentin.capdet@gmail.com"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.
    
    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences
    """
    with open(fastq_file, "r") as file:
        while True:
            try:
                next(file)  
                sequence = next(file).strip()  
                next(file)  
                next(file)  
                
                yield sequence
            except StopIteration:
                break



def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    for i in range(0, len(read) - kmer_size + 1):
        yield read[i:i + kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_dict = {}

    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1

    return kmer_dict 




def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    graph = nx.DiGraph()

    for kmer, count in kmer_dict.items():
        prefix = kmer[:-1]  
        suffix = kmer[1:]  

        if prefix not in graph:
            graph.add_node(prefix)
        if suffix not in graph:
            graph.add_node(suffix)

        graph.add_edge(prefix, suffix, weight=count)

    return graph



def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    #print(path_list)
    #print(delete_entry_node, delete_sink_node)
    #Graph.remove_nodes_from(path) #supprime les noeuds d’un chemin
    for path in path_list:
        if delete_entry_node == True and delete_sink_node == True:
            graph.remove_nodes_from(path)
        elif delete_entry_node == False and delete_sink_node == False:
                graph.remove_nodes_from(path[1:-1])
        elif delete_entry_node == True:
                graph.remove_nodes_from(path[:-1])
        elif delete_sink_node == True:
                graph.remove_nodes_from(path[1:])
    return graph



def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """

    # Calculate standard deviation for weight average and path length
    weight_stddev = statistics.stdev(weight_avg_list)
    length_stddev = statistics.stdev(path_length)
    #print(path_list)
    # If weight_stddev is greater than 0, select the path with the highest weight_avg
    if weight_stddev > 0:
        best_path_index = weight_avg_list.index(max(weight_avg_list))
    # If weight_stddev is 0 and length_stddev is greater than 0, select the longest path
    elif length_stddev > 0:
        best_path_index = path_length.index(max(path_length))
    # If both stddevs are 0, choose randomly among paths
    else:
        best_path_index = random.randint(0, len(path_list) - 1)

    best_path = path_list[best_path_index]
    #print(path_list[best_path_index])

    # Remove unwanted nodes
    for path in path_list:
       if path != best_path:
           remove_path = path.copy()
           if not delete_entry_node:
               remove_path.pop(0)
           if not delete_sink_node:
               remove_path.pop(-1)
           graph.remove_nodes_from(remove_path)

    return graph
    
    

def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    
    """ solve_bubble sera appelée lorsqu’une bulle sera  détectée entre un nœud ancêtre et un nœud descendant. 
    Il devra déterminer les chemins “simples” possible entre deux nœuds (un ancêtre et un descendant) et 
    calculer la longueur et le poids de ces chemins. Puis fait appel à select_best_path pour  faire le choix final. """
   
    path_data = []

    for path in list(nx.all_simple_paths(graph, ancestor_node, descendant_node)):
        path_avg_wgt = path_average_weight(graph, path)
        path_length = sum([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])
        path_data.append((path, path_avg_wgt, path_length))
    
    graph = select_best_path(graph, [p[0] for p in path_data], [p[1] for p in path_data], [p[2] for p in path_data])
   
    
    return graph

def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble = False
    
    for node in graph.nodes():
        predecessors = list(graph.predecessors(node))
        
        if len(predecessors) > 1:
            for i in range(len(predecessors) - 1):
                for j in range(i + 1, len(predecessors)):
                    ancestor_node = nx.lowest_common_ancestor(graph, predecessors[i], predecessors[j])
                    if ancestor_node != None:
                        bubble = True
                        break
            if bubble == True:
                break
    
    if bubble:
        graph = solve_bubble(graph, ancestor_node, node)
    
    return graph

def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    nodes_with_multi_pred = [node for node, degree in graph.in_degree() if degree > 1]
    
    for node in nodes_with_multi_pred:
        pred = list(graph.predecessors(node))
        
        if any(predecessor in starting_nodes for predecessor in pred):
            # Extract paths from starting nodes to the current node
            paths_from_starting = [list(nx.all_simple_paths(graph, source=start, target=node))[0] for start in starting_nodes if nx.has_path(graph, start, node)]
            
            # Calculate the weights of these paths
            path_wgt = []
            for path in paths_from_starting:
                wgts = [graph[path[i]][path[i+1]].get('weight', 1) for i in range(len(path)-1)]
                path_wgt.append(sum(wgts))
            
            # Determine the main path
            main_path_index = path_wgt.index(max(path_wgt))
            
            # Remove other paths
            for i, path in enumerate(paths_from_starting):
                if i != main_path_index:
                    for j in range(len(path) - 1):
                        graph.remove_edge(path[j], path[j+1])
    
    return graph

def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
     # Pour chaque nœud du graphe
    for node in graph.nodes():
        # Récupération des successeurs
        successors = list(graph.successors(node))
        
        # Si le nœud a plus d'un successeur
        if len(successors) > 1:
            # Récupération des poids des arêtes vers les successeurs
            weights = [graph[node][succ]['weight'] for succ in successors]
            
            # On supprime l'arête avec le poids le plus faible
            min_weight = min(weights)
            for succ in successors:
                if graph[node][succ]['weight'] == min_weight:
                    graph.remove_edge(node, succ)
                    break # On sort de la boucle dès la suppression

    return graph

def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    starting_nodes = []
    for node in graph.nodes():
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes.append(node)

    return starting_nodes

def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    ending_nodes = []
    for node in graph.nodes():
        if len(list(graph.successors(node))) == 0:
            ending_nodes.append(node)

    return ending_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contigs_list = []

    for start in starting_nodes:
        for end in ending_nodes:
            if nx.has_path(graph, start, end):
                for path in nx.all_simple_paths(graph, start, end):
                    contig = path[0]
                    for i in range(1, len(path)):
                        contig = contig + path[i][-1]
                    
                    contigs_list.append((contig, len(contig)))
    return contigs_list


def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    with open(output_file, "w") as file:
        for i, val in enumerate(contigs_list):
            contig = val[0]
            wrap = textwrap.fill(contig, width=80)
            file.write(f'>contig_{i} len={val[1]}\n{wrap}\n')



def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """

    # Get arguments
    args = get_arguments()

    #seq = read_fastq(args.fastq_file)
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    graph = simplify_bubbles(graph)
    
    start = get_starting_nodes(graph)
    end = get_sink_nodes(graph)

    graph = solve_entry_tips(graph, start)
    graph = solve_out_tips(graph, end)
    
    start = get_starting_nodes(graph)
    end = get_sink_nodes(graph)
    contigs = get_contigs(graph, start, end)

    
    save_contigs(contigs, args.output_file)
 
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    #if args.graphimg_file:
    #    draw_graph(graph, args.graphimg_file)


if __name__ == '__main__': # pragma: no cover
    main()
