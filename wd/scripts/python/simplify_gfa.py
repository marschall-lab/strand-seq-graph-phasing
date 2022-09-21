import argparse
import pdb
import csv

import gfapy
from module_simplify_gfa.simplify_gfa_functions import *

parser = argparse.ArgumentParser(description='Simplify GFA.')
parser.add_argument('--input', nargs=1, type=str, required=True)
parser.add_argument('--segment-length-threshold', nargs=1, type=int, required=True)
parser.add_argument('--output', nargs=1, type=str, required=True)

parser.add_argument('--clusters', nargs=1, type=str, required=False)
args = parser.parse_args()
print(args)

# TODO 
graph = gfapy.Gfa.from_file(args.input[0], vlevel = 3)
# graph = gfapy.Gfa.from_file('HG00211/assembly.homopolymer-compressed_graph_components/simplified_assembly.gfa', vlevel = 3)
# graph = gfapy.Gfa.from_file('input/HG002/assembly/assembly.homopolymer-compressed.gfa', vlevel = 3)

#TODO
thresh=args.segment_length_threshold[0]
# thresh=1000000
# thresh=99999999999

def simplify_graph(graph, thresh, segments_to_keep_names=None, macro_keep_only=False, macro_dead_end_names=None):
    print('macro simplify')
    simplified_graph = macro_simplify(graph, thresh, segments_to_keep_names, macro_keep_only, macro_dead_end_names)
    simplified_graph = remove_self_links(simplified_graph)
    simplified_graph = remove_solo_segments(simplified_graph, thresh, segments_to_keep_names)
    # simplified_graph.to_file('test_all1.gfa')
    
    print('fine simplify')
    simplified_graph = remove_self_links(simplified_graph)
    simplified_graph = remove_nodes(simplified_graph, thresh, segments_to_keep_names) 
    # simplified_graph.to_file('test_all2.gfa')
    
    print('cleanup and crimp')
    simplified_graph = remove_self_links(simplified_graph)
    simplified_graph = remove_dead_ends(simplified_graph, thresh, segments_to_keep_names)
    # simplified_graph.to_file('test_all3.gfa')
    simplified_graph = crimp(simplified_graph)
    
    return simplified_graph
    
# Empty graph check:
if not graph.segments and not graph.dovetails:
    graph.to_file(args.output[0])
elif args.clusters is not None:
    # TODO
    clusts_file = args.clusters[0]
    # clusts_file =  'HG002/SaaRclust/Clusters/MLclust.data'

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
    # could make this faster by only passing the dead ons on the conneted
    # components that contain the keep_names. then would have a different set of
    # dead ends for each keep, depending on that connected components are
    # coveres
    dead_end_segment_names =  [seg.name for seg in graph.segments if is_dead_end(seg)]
    print('calulated dead ends')
    # keep_names = clust_rname['17']
    simplified_graphs = [simplify_graph(graph, thresh=thresh, segments_to_keep_names=keep_names, macro_keep_only=True, macro_dead_end_names=dead_end_segment_names) for keep_names in clust_rname.values()]
    simplified_graph = graph_from_graph_list(simplified_graphs, overlap=False)
else:
    simplified_graph = simplify_graph(graph, thresh=thresh)
   
print('export')
# TODO
simplified_graph.to_file(args.output[0])
# simplified_graph.to_file('test_graph3.gfa')
    
  
  
# for seg in graph.segments: 
#     if seg.length <= 1000000: 
#         short_names.append(seg.name)
#         
# with open('short_names.txt', 'w') as f:
#     for item in short_names:
#         f.write("%s\n" % item)
