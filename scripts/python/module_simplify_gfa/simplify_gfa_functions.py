# TODO, make the X and Z util and crimping less yucky feeling
# TODO make it clear which functions happen in place and which don't
import itertools
import pdb

import gfapy

## Pass-Through Remove
def pass_through_remove(segment):
    graph=segment.gfa
    # Need name and orientation. Merging "from" nodes in "from" position to "to"
    # node in "to" positions
    # from----to----from-------to <- Use segment to setup orientation
    # L-------segment----------R
    # from---------------------to <- Final desired outcome

    dt_l = list()
    for edge in segment.dovetails_L:
        if is_self_link(edge):
            continue
        if edge.from_name == segment.name:
            dt_l.append(edge.complement())
        else:
            dt_l.append(edge.clone())
    dt_r = list()
    for edge in segment.dovetails_R:
        if is_self_link(edge):
            continue
        if edge.from_name == segment.name:
            dt_r.append(edge.clone())
        else:
            dt_r.append(edge.complement())

    dead_end = (not dt_r) or (not dt_l)

    if dead_end:
        return False

    for l, r in itertools.product(dt_l, dt_r):
        if l.from_name != r.to_name: # don't make self links? Maybe make self links?

            line = f'L\t{l.from_name}\t{l.from_orient}\t{r.to_name}\t{r.to_orient}\t0M'
            added = add_line_unique_safe(graph, line)
            if added:
                print(line)

    print(segment.name)
    graph.rm(segment.name)
    return True


###### Node Insertion #######

def edge_orientation_in_position(segment, side, from_to):
    assert(side in ['L', 'R'])
    assert(from_to in ['from', 'to'])
    # Only need one edge to get orientation, if all edges from same side
    edge, *_ = segment.dovetails_of_end(side)
    if from_to == 'from':
        out = edge.from_orient if segment.name == edge.from_name else edge.complement().from_orient
    else:
        out = edge.to_orient if segment.name == edge.to_name else edge.complement().to_orient
    return out

#####################
###### X ############
#####################

def __x_util(segment, side):
    # pdb.set_trace()
    assert(side in ['L', 'R'])

    segment_edges = segment.dovetails_of_end(side)
    detected = False
    if len(segment_edges) == 2:
        segment_end = gfapy.SegmentEnd(segment, side)
        n1, n2 = [e.other_end(segment_end) for e in segment_edges]
        n1_edges = n1.dovetails_of_end(n1.end_type)
        n2_edges = n2.dovetails_of_end(n2.end_type)
        n1_neigh = [e.other_end(n1) for e in n1_edges]
        n2_neigh = [e.other_end(n2) for e in n2_edges]
        if len(n1_neigh) == 2 and len(n2_neigh) == 2 and n1 not in n2_neigh and n2 not in n1_neigh:
            n13, = [x for x in n1_neigh if x != segment_end]
            n23, = [x for x in n2_neigh if x != segment_end]
            if n13 == n23:
                n3 = n13 # for naming consistency
                n3_edges = n3.dovetails_of_end(n3.end_type)
                # n3_neigh = getattr(n3, f"neighbours_{n3_side}")
                if len(n3_edges) == 2:
                    detected = True

    if not detected:
        return (detected, None)
    # Orientations
    # segment, n3 ~ from nodes | n1, n2 ~ to nodes
    segment_orient = edge_orientation_in_position(segment, side, 'from')
    n3_orient = edge_orientation_in_position(n3.segment, n3.end_type, 'from')
    n1_orient = edge_orientation_in_position(n1.segment, n1.end_type, 'to')
    n2_orient = edge_orientation_in_position(n2.segment, n2.end_type, 'to')
    from_out = [(segment.name, segment_orient, side),(n3.name, n3_orient, n3.end_type)]
    to_out = [(n1.name, n1_orient, n1.end_type),(n2.name, n2_orient, n2.end_type)]
    return (detected, [from_out,to_out])


#####################
###### Z ############
#####################

def __z_util(segment, side):
    # pdb.set_trace()
    assert(side in ['L', 'R'])
    segment_edges = segment.dovetails_of_end(side)
    detected = False
    if len(segment_edges) == 2:
        segment_end = gfapy.SegmentEnd(segment, side)
        n1, n2 = [e.other_end(segment_end) for e in segment_edges]
        n1_edges = n1.dovetails_of_end(n1.end_type)
        n2_edges = n2.dovetails_of_end(n2.end_type)
        n1_neigh = [e.other_end(n1) for e in n1_edges]
        n2_neigh = [e.other_end(n2) for e in n2_edges]
        if ((len(n1_neigh) == 2 and len(n2_neigh) == 1) or (len(n1_neigh) == 1 and len(n2_neigh) == 2)) and n1 not in n2_neigh and n2 not in n1_neigh:
            if len(n1_neigh) == 2:
                n3_source = n1
                n3_source_edges = n1_edges
            else:
                n3_source = n2
                n3_source_edges = n2_edges
            n3, = [e.other_end(n3_source) for e in n3_source_edges if e.other_end(n3_source) != segment_end]
            n3_edges = n3.dovetails_of_end(n3.end_type)
            if len(n3_edges) == 1:
                detected = True
    if not detected:
        return (detected, None)
    # Orientations
    # segment, n3 ~ from nodes | n1, n2 ~ to nodes
    segment_orient = edge_orientation_in_position(segment, side, 'from')
    n3_orient = edge_orientation_in_position(n3.segment, n3.end_type, 'from')
    n1_orient = edge_orientation_in_position(n1.segment, n1.end_type, 'to')
    n2_orient = edge_orientation_in_position(n2.segment, n2.end_type, 'to')
    from_out = [(segment.name, segment_orient, side),(n3.name, n3_orient, n3.end_type)]
    to_out = [(n1.name, n1_orient, n1.end_type),(n2.name, n2_orient, n2.end_type)]
    return (detected, [from_out,to_out])

#####################
###### INSERT #######
#####################
# X and Z insert might be amenable to some sort of scaffold that takes the
# detector and edge function and performs the node crimping, as both are
# examples of connecting 2 "from" nodes and 2 "to" nodes
def __two_two_crimp(f_util, prefix):
    n=0
    def _two_two_crimp(segment, side):
        nonlocal n
        detected, out = f_util(segment, side)
        if not detected:
            pass
        else:
            # pdb.set_trace()
            from_, to_ = out
            graph = segment.gfa

            # pdb.set_trace()
            for name, _, lr in from_+to_:
                edges = graph.segment(name).dovetails_of_end(lr)
                # with a for loop, disconnecting an edge modifies the list of
                # edges while iterating over that list of edges ~ can lead to
                # bugs. This is why while ~ pop is used
                while(edges):
                    edge = edges.pop()
                    edge.disconnect()

            tip_name = f'{prefix}tip-{n}'
            n += 1
            # new tip. The Sequence field is optional and can be *, meaning that
            # the nucleotide sequence of the segment is not specified
            graph.add_line(f'S\t{tip_name}\t*')
            for name, orient, lr in from_:
                graph.add_line(f'L\t{name}\t{orient}\t{tip_name}\t+\t0M')
            for name, orient, lr in to_:
                graph.add_line(f'L\t{tip_name}\t+\t{name}\t{orient}\t0M')

    return _two_two_crimp

######## Simple Tip Merging ##########
def is_simple_dead_end(seg):
    connectivity =  [len(seg.neighbours_L), len(seg.neighbours_R)]
    simple_dead_end = 0 in connectivity and 1 in connectivity
    return simple_dead_end #  and not seg.gfa.is_cut_segment(seg)

def has_two_simple_tips(seg, side):
    assert(side in ['L', 'R'])
    neighbours = seg.neighbours_of_end(side)
    return len(neighbours) == 2 and all([is_simple_dead_end(x) for x in neighbours])

# TODO refactor using segment ends?
def __e_crimp():
    n=0
    def _e_crimp(segment, side):
        nonlocal n
        if not has_two_simple_tips(segment, side):
            pass
        else:
            # pdb.set_trace()
            graph = segment.gfa
            tip_name = f'etip-{n}'
            n += 1

            seg1, seg2 = segment.neighbours_of_end(side)
            # Get orientation by taking opposite of side that actually has edges
            side1 = 'R' if len(seg1.neighbours_L) == 0 else 'L'
            side2 = 'R' if len(seg2.neighbours_L) == 0 else 'L'

            orient1 = edge_orientation_in_position(seg1, side1, 'from')
            orient1 = '-' if orient1 == '+' else '+'

            orient2 = edge_orientation_in_position(seg2, side2, 'from')
            orient2 = '-' if orient2 == '+' else '+'
            # new tip. The Sequence field is optional and can be *, meaning that
            # the nucleotide sequence of the segment is not specified
            graph.add_line(f'S\t{tip_name}\t*')
            graph.add_line(f'L\t{seg1.name}\t{orient1}\t{tip_name}\t+\t0M')
            graph.add_line(f'L\t{seg2.name}\t{orient2}\t{tip_name}\t+\t0M')
    return _e_crimp


x_insert = __two_two_crimp(__x_util, 'x')
z_insert = __two_two_crimp(__z_util, 'z')
e_insert = __e_crimp()


### Should these functions live here or in the main script? ###

# For hashability
def end_as_tuple(x):
    return (x.name, x.end_type)

def tuple_as_end(x, graph):
    return gfapy.SegmentEnd(graph.segment(x[0]), x[1])


# Self Links
def is_self_link(edge):
    return edge.from_segment == edge.to_segment

def remove_self_links(graph):
    # I have to write this function because so many of the gfapy builtins are buggy
    edges = graph.edges
    edges_to_keep = list()
    while(edges):
        edge = edges.pop()
        if not is_self_link(edge):
            edges_to_keep.append(edge.clone())

    new_graph = gfapy.Gfa(vlevel=3)
    for segment in graph.segments:
        new_graph.add_line(segment.clone())
    for edge in edges_to_keep:
        new_graph.add_line(edge)

    return(new_graph)

def get_connected_ends(segment_end):
    edges = segment_end.dovetails_of_end(segment_end.end_type)
    connecting_ends = set()
    for edge in edges:
        other_end = edge.other_end(segment_end)
        connecting_ends.add(end_as_tuple(other_end))
    return connecting_ends

def remove_if_not_cut(graph, seg_name):
    segment = graph.segment(seg_name)
    backup_lines = [segment.clone()] + [x.clone() for x in segment.all_references]
    num_cc_before = len(graph.connected_components())
    graph.rm(segment)
    if len(graph.connected_components()) > num_cc_before:
        for line in backup_lines:
            add_line_unique_safe(graph, line)
        return False
    else:
        return True

def add_line_unique_safe(graph, line):
    added = True
    try:
        graph.add_line(line)
    except gfapy.NotUniqueError:
        added = False
    return added

# The whole thing with dead end segments is stupid
def macro_simplify(
  graph,
  threshold=None,
  keep=None):
    #  keep ~ segment names
    if keep is None and  threshold is None:
        raise ValueError('no segments to remove')

    if threshold is None:
        segments_to_keep_names = keep
    else:
        segments_to_keep_names = set()
        for seg in graph.segments:
            if seg.length >= threshold:
                segments_to_keep_names.add(seg.name)

        if keep is not None:
            segments_to_keep_names = segments_to_keep_names.union(keep)

    segments_to_keep_neighbour_names = set()
    for segment_name in segments_to_keep_names:
        segment = graph.segment(segment_name)
        for neigh in segment.neighbours:
            segments_to_keep_neighbour_names.add(neigh.name)
        
    dead_end_segments_to_keep_names =  [seg.name for seg in graph.segments if is_dead_end(seg)]

    terminal_segment_names = segments_to_keep_neighbour_names.union(segments_to_keep_names).union(dead_end_segments_to_keep_names)

    simplified_graph = gfapy.Gfa(vlevel=3)

    for seg_name in terminal_segment_names:
        simplified_graph.add_line(graph.segment(seg_name).clone())

    for segment_name in segments_to_keep_neighbour_names:
        print(segment_name)
        segment = graph.segment(segment_name)
        for side in ['L', 'R']:
            segment_edges = segment.dovetails_of_end(side)
            segment_end = gfapy.SegmentEnd(segment, side)
            next_ends = get_connected_ends(segment_end)
            visited_ends = set()
            while(any(next_ends)):
                end = next_ends.pop()
                visited_ends.add(end)
                if end[0] not in terminal_segment_names :
                    other_seg_ends = get_connected_ends(tuple_as_end(end, graph).inverted())
                    other_seg_ends = {x for x in other_seg_ends if x not in visited_ends}
                    next_ends = next_ends.union(other_seg_ends)

            visited_terminal_ends = [x for x in visited_ends if x[0] in terminal_segment_names]
            for to_end in visited_terminal_ends:
                segment_orient = edge_orientation_in_position(segment_end.segment, segment_end.end_type, 'from')
                to_seg_end = tuple_as_end(to_end, graph)
                to_name = to_seg_end.name
                to_orient = edge_orientation_in_position(to_seg_end.segment, to_seg_end.end_type, 'to')
                line = f'L\t{segment_name}\t{segment_orient}\t{to_name}\t{to_orient}\t0M'
                add_line_unique_safe(simplified_graph, line)

    return simplified_graph

def refresh_graph(graph, removed_segment_names):
    # removing or disconnecting a segment isn't clean: Sometimes (Generally) it
    # fails to properly remove the segment and its associated edges, leaving
    # improper nodes that should have been deleted and duplicated edges. This is
    # not fixed by a simple IO step, unfortunately.
    graph2 = gfapy.Gfa(vlevel=3)
    for segment in graph.segments:
        if segment.name not in removed_segment_names:
            graph2.add_line(segment.clone())
    segment_names = [x.name for x in graph2.segments]
    for edge in graph.dovetails:
        if edge.from_name in segment_names and edge.to_name in segment_names:
            add_line_unique_safe(graph2, edge.clone())
    return graph2


# Solo segments
def is_solo(seg):
    connectivity =  [len(seg.neighbours_L), len(seg.neighbours_R)]
    return 0 == connectivity[0] and  0 == connectivity[1]

def remove_solo_segments(graph, thresh = None, keep = None):
    # solo_segment_names =  {seg.name for seg in graph.segments if is_solo(seg) and seg.length and seg.length < thresh}
    # 
    # 
    # if keep is not None:
    #     solo_segment_names = solo_segment_names.difference(keep)

    solo_segments_to_remove =  [seg for seg in graph.segments if is_solo(seg) and seg.length]
    if thresh is not None:
        solo_segments_to_remove = [seg for seg in solo_segments_to_remove if seg.length < thresh]    
    if keep is not None:
        solo_segments_to_remove = [seg for seg in solo_segments_to_remove if seg.name not in keep]
        
    solo_segment_names =  {seg.name for seg in solo_segments_to_remove}
    removed_nodes=set()
    for seg_name in solo_segment_names:
        print(seg_name)
        graph.rm(seg_name)
        removed_nodes.add(seg_name)

    return refresh_graph(graph, removed_nodes)

def is_dead_end(seg):
    connectivity =  [len(seg.neighbours_L), len(seg.neighbours_R)]
    return 0 in connectivity #  and not seg.gfa.is_cut_segment(seg)

def remove_dead_ends(graph, thresh = None, keep = None):
    # TODO I don't like the way this function is written, I would like to avoid
    # the repeated refresh_graph calls, hopefully by doing this additively
    # instead of subtractively
    dead_end_segments_to_remove =  [seg for seg in graph.segments if is_dead_end(seg) and seg.length]
    
    if thresh is not None:
        dead_end_segments_to_remove = [seg for seg in dead_end_segments_to_remove if seg.length < thresh]
    if keep is not None:
        dead_end_segments_to_remove = [seg for seg in dead_end_segments_to_remove if seg.name not in keep]

    dead_end_segments_to_remove_names =  {seg.name for seg in dead_end_segments_to_remove}
    
    
    removed_nodes=set()
    any_removed = True
    while any_removed:
        any_removed = False
        for name in dead_end_segments_to_remove_names:
            removed = remove_if_not_cut(graph, name)
            if removed:
                any_removed = True
                print(name)
                removed_nodes.add(name)
        graph = refresh_graph(graph, removed_nodes)
        
        dead_end_segments_to_remove =  [seg for seg in graph.segments if is_dead_end(seg) and seg.length]
        if thresh is not None:
            dead_end_segments_to_remove = [seg for seg in dead_end_segments_to_remove if seg.length < thresh]
        if keep is not None:
            dead_end_segments_to_remove = [seg for seg in dead_end_segments_to_remove if seg.name not in keep]
        dead_end_segments_to_remove_names =  {seg.name for seg in dead_end_segments_to_remove}
    

    return refresh_graph(graph, removed_nodes)

def remove_nodes(graph, threshold=None, keep = None):

    # iterate over list of names, not segments, as disconnecting segments while
    # looping over list of them them can lead to problems with the list being
    # modified during iteration
    if keep is None and  threshold is None:
        raise ValueError('no segments to keep')
      
    segments_to_keep_names = set()
    if threshold is not None:
        for seg in graph.segments:
            if seg.length >= threshold:
                segments_to_keep_names.add(seg.name)

    if keep is not None:
        segments_to_keep_names = segments_to_keep_names.union(keep)

    segments_to_remove_names = set()
    for name in graph.segment_names:
        if name not in segments_to_keep_names:
            segments_to_remove_names.add(name)
    # TODO is this step still needed?
    # Order segments_to_remove_names to remove segments adjacent to (nodes being
    # kept) last. this may help more accurately preserve the presence of tip
    # structures during pass through removal
    segments_to_keep_neighbour_names = set()
    for name in segments_to_keep_names:
        segment = graph.segment(name)
        for neigh in segment.neighbours:
            segments_to_keep_neighbour_names.add(neigh.name)

    segments_to_remove_names = list(segments_to_remove_names.difference(segments_to_keep_neighbour_names)) +\
                          list(segments_to_remove_names.intersection(segments_to_keep_neighbour_names))

    removed_segments = set()
    for name in segments_to_remove_names:
        removed = pass_through_remove(graph.segment(name))
        if removed:
            removed_segments.add(name)

    graph2 = refresh_graph(graph, removed_segments)
    return graph2

def unstar_edges(graph):
    for edge in graph.dovetails:
        if edge.field_to_s('overlap') == '*':
            edge.disconnect()
            edge.overlap = '0M'
            edge.connect(graph)

def crimp(graph):
    for segment in graph.segments:
        print(segment.name)
        print('x-left')
        x_insert(segment, 'L')
        print('x-right')
        x_insert(segment, 'R')
        print('z-left')
        z_insert(segment, 'L')
        print('z-right')
        z_insert(segment, 'R')
    # so after adding segments and removing edges, the references in general
    # become messed up, but it appears that some sort of i/o step properly fixes
    # things?
    graph2=gfapy.Gfa(str(graph), vlevel = 3)

    for segment in graph2.segments:
        print(segment.name)
        print('e-left')
        e_insert(segment, 'L')
        print('e-right')
        e_insert(segment, 'R')
    # For some reason, some edges are ending up with '*' overlap
    unstar_edges(graph2)
    # for edge in graph2.dovetails:
    #     if edge.field_to_s('overlap') == '*':
    #         edge.disconnect()
    #         edge.overlap = '0M'
    #         edge.connect(graph2)

    return graph2


################## Subgraph Simplification #######################

def graph_from_segment_list(l):
    new_graph = gfapy.Gfa(vlevel=3)
    for segment in l:
        refs = segment.all_references
        new_graph.add_line(segment.clone())
        for ref in refs:
            add_line_unique_safe(new_graph, ref.clone())
    return new_graph

def graph_from_graph_list(l, overlap=True):
    # overlap ~ allow the graphs being added to merge to together if they share nodes
    new_graph = gfapy.Gfa(vlevel=3)
    for graph in l:
        for line in graph.segments:
            added = add_line_unique_safe(new_graph, line.clone())
            while not added and not overlap:
                # TODO a more intelligent string ßto add? I am currenly adding a
                # number to try to not break regexs, as all the unitig names
                # appear to end in numbers
                line.name += '9'
                added = add_line_unique_safe(new_graph, line.clone())
        for line in graph.edges:
            add_line_unique_safe(new_graph, line.clone())

    return new_graph


def extract_components_containing_segment_names(graph, names):
    components_containing_names = [component for component in graph.connected_components() if any([seg.name in names for seg in component])]
    segment_list=sum(components_containing_names, [])
    sub_graph= graph_from_segment_list(segment_list)
    return sub_graph
