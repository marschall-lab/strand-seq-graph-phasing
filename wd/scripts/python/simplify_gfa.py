# Library -----------------------------------------------------------------
import argparse
import pdb
import csv

import gfapy
from module_simplify_gfa.simplify_gfa_functions import *


# Parsing -----------------------------------------------------------------
parser = argparse.ArgumentParser(description='Simplify GFA.')
parser.add_argument('--input', nargs=1, type=str, required=True)

parser.add_argument('--output', nargs=1, type=str, required=True)

parser.add_argument('--segment-length-threshold', nargs=1, type=int, required=False)
parser.add_argument('--clusters', nargs=1, type=str, required=False)
args = parser.parse_args()
print(args)


thresh = None
clust_file = None

if args.segment_length_threshold is not None:
    thresh = args.segment_length_threshold[0]
    
if args.clusters is not None:
    clust_file = clusts_file = args.clusters[0]
    
if thresh is None and clust_file is None:
    raise ValueError('At least one of "clusters" and "segment_length_threshold" must be input')


# Import ------------------------------------------------------------------

graph = gfapy.Gfa.from_file(args.input[0], vlevel = 3)

def simplify_graph(graph, thresh=None, segments_to_keep_names=None):
    # remove self links is probably called more than necessary, but it fixed bugs.
    
    simplified_graph = remove_self_links(graph)
    simplified_graph = remove_solo_segments(simplified_graph, thresh, segments_to_keep_names)
    
    print('macro simplify')
    simplified_graph = remove_self_links(simplified_graph)
    simplified_graph = macro_simplify(simplified_graph, thresh, segments_to_keep_names)
    simplified_graph = remove_self_links(simplified_graph)
    simplified_graph = remove_solo_segments(simplified_graph, thresh, segments_to_keep_names)
    
    print('fine simplify')
    simplified_graph = remove_self_links(simplified_graph)
    simplified_graph = remove_nodes(simplified_graph, thresh, segments_to_keep_names) 
    
    print('cleanup and crimp')
    simplified_graph = remove_self_links(simplified_graph)
    simplified_graph = remove_dead_ends(simplified_graph, thresh, segments_to_keep_names)
    
    return simplified_graph
    
def simplify_and_crimp_graph(graph, thresh=None, segments_to_keep_names=None):
    simplified_graph= simplify_graph(graph, thresh, segments_to_keep_names)
    simplified_graph = crimp(simplified_graph)
    
    return simplified_graph

# Main --------------------------------------------------------------------

# Empty graph check:
if not graph.segments and not graph.dovetails:
    graph.to_file(args.output[0])

elif args.clusters is not None:
    # TODO
    clusts_file = args.clusters[0]

    clust_rname = dict()
    with open(clusts_file, 'r') as file:
        clusts = csv.reader(file, delimiter="\t")
        header = next(clusts)
        rname_ix = header.index('#rname')
        clust_ix = header.index('chrom_clust')
        for row in clusts:
            rname = row[rname_ix]
            clust = row[clust_ix]
            if clust not in clust_rname:
                clust_rname[clust] = set()
            
            clust_rname[clust].add(rname)
            print(row)
         
    all_clustered_rnames = []
    for x in clust_rname.values():
        all_clustered_rnames.extend(x)
    # Speed things up by doing an all over simplify before, which reduces a 
    # lot of redundant work when simplifying each cluster
    simplified_graph = simplify_graph(graph, thresh=thresh, segments_to_keep_names=all_clustered_rnames) 
    
    simplified_graphs = [simplify_and_crimp_graph(simplified_graph, thresh=thresh, segments_to_keep_names=keep_names) for keep_names in clust_rname.values()]
    simplified_graph = graph_from_graph_list(simplified_graphs, overlap=False)
else:
    simplified_graph = simplify_and_crimp_graph(graph, thresh=thresh)
   
print('export')

simplified_graph.to_file(args.output[0])

    
  
