#!/bin/env python3

import sys
import os
#import pygraphviz as pgv

from gfapy import Gfa
from typing import List

#def gfa_to_dot(ingfa: Gfa) -> pgv.AGraph:#List[str]:
#    dot = pgv.AGraph(directed=True, landscape=True)
#    for 

def log(s: str):
    print(s, file=sys.stderr)

def gfa_to_dot(gfa, include_len=False, skip_missing=True):
    # Initialize DOT format
    dot_lines = ["digraph GFA_graph {"]
    dot_lines.append("  rankdir=LR;")
    dot_lines.append("  node [style=filled, fillcolor=lightblue];")
    dot_lines.append("  edge [color=darkgreen];")
    
    # Add segments as nodes
    for segment in gfa.segments:
        dot_lines.append(f'  {segment.name} [label={segment.name}];')
    
    # Add edges as connections
    for edge in gfa.edges:
        # seems that gfapy autocompletes these, not sure if possible to disable
        ## ignore lines outside the file
        ## for plotting an incomplete subset of a gfa
        #if skip_missing and not (edge.from_segment in gfa.segments and edge.to_segment in gfa.segments):
        #    continue

        dot_lines.append(f'  {edge.from_segment.name} -> {edge.to_segment.name};')
    
    dot_lines.append("}")
    return "\n".join(dot_lines)

#def save_to_file(txt: List[str], fname: str):
#    with open(fname, "wt") as f:
#        f.writelines(txt)

if __name__ == "__main__":
    # read in GFA, either from 1st arg or stdin
    try:
        #gfa = Gfa.from_file(sys.argv[1] if (len(sys.argv) > 1 and not sys.argv[1]=='-') else "/dev/stdin", vlevel=0) #sys.stdin.read(), vlevel=0)
        gfa = Gfa.from_file(sys.argv[1], vlevel=0) if (len(sys.argv) > 1 and not sys.argv[1]=='-') else Gfa(sys.stdin.read()[:-1], vlevel=0) # throws very weird errors with the terminating char for some reason
        log(f"Loaded GFA with {len(gfa.segments)} segments and {len(gfa.edges)} edges")
    except Exception as e:
        log(f"Error loading GFA file!")
        log(e)
        sys.exit(-1)

    dot = gfa_to_dot(gfa)

    # emit dot
    with open(sys.argv[2], "wt") if len(sys.argv) > 2 and not sys.argv[2]=='-' else sys.stdout as f:
        f.write(dot)
        f.write("\n")

    #TODO have fn to convert into pdf directly?
