#!/usr/bin/env python3
"""
Plot haplotype graph and the distribution of kmers in the nodes of the graph
"""
import argparse
import os


def cleanax(ax):
    """
    Remove the top and right spine and set plot to have tight layout
    """
    from matplotlib import pyplot as plt
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # ax.legend(bbox_to_anchor=(1.01, 1))
    plt.tight_layout(pad=0.1)
    return ax
# END


def bezierlink(p1, p2, z):
    import matplotlib.patches as patches
    from matplotlib.path import Path

    rs = p1[0]
    ry = p1[1]
    qs = p2[0]
    qy = p2[1]

    smid = (qs - rs) / 2  # Start increment
    hmid = (qy - ry) / 2  # Height increment
    verts = [(rs, ry),
             (rs + smid, ry),
             (rs + smid, ry + 2 * hmid),
             (rs + 2 * smid, ry + 2 * hmid),
             ]
    codes = [
        Path.MOVETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor="none", zorder=z)
    return patch
# END
# </editor-fold>


def plot_haplotype_graph(args):
    """
    Plot a given a haplotype graph
    :return:
    """
    from collections import defaultdict, deque, Counter
    from itertools import product
    import igraph as ig
    import string
    import math
    import numpy as np
    from numpy.polynomial import Polynomial as P
    from matplotlib import pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.collections import PatchCollection
    import pandas as pd
    from hometools.hometools import mylogger
    from hometools.plot import cleanax
    import os
    from hometools.plot import plotdensity

    class hapobject:
        """
        A haplotype block object
        """

        def __init__(self, id, start, end, genomes):
            self.id = id
            self.start = start
            self.end = end
            self.genomes = genomes
        # END

        def hasgen(self, gen):
            return gen in self.genomes
        # END
    # END

    def readthreads(f):
        listtoint = lambda x: list(map(int, x.replace('[', '').replace(']', '').split(', ')))
        threads = deque()
        with open(f, 'r') as fin:
            fin.readline()
            for line in fin:
                line = line.strip().split('\t')
                ts = listtoint(line[2])
                # Only select threads that cover at least three nodes
                if len(ts) >= 3:
                    threads.append([int(line[1]), listtoint(line[2])])
        return threads

    hapfin = args.hap.name
    thfin = args.threads.name
    outfin = args.out.name
    spread = args.spread
    W = args.W
    H = args.H

    cwd = os.getcwd()
    # emoutfin = "/home/ra98jam/d16/test/WhiteRose_results_2024_06_14/EM_results/EM.v03.WhR.w0.results.tsv"
    # hapfin = "../../data/haplotype_graph_2024_06_14.txt"


    logger = mylogger(__name__)
    # chrsize = 46102915 # Currently, set for Chr2
    # logger.info(f"Setting chromosome size to {chrsize}")
    samples = ['dm'] + [f'{i}_hap{j}' for i, j in product(string.ascii_uppercase[:10], range(5, 9))]

    # s = args.start
    # e = args.end
    # s = 0
    # e = chrsize

    os.chdir(cwd)

    hapoblist = deque()
    with open(hapfin, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            hapoblist.append(hapobject(int(line[3]), int(line[0]), int(line[1]), sorted(line[2].split(','))))

    # Create a Graph object for the haplotype blocks and use it to determine the Y-coordinate for visualisation of the blocks
    G = ig.Graph()
    G = ig.Graph.as_directed(G)
    G.add_vertices(range(len(hapoblist)))
    addededge = set()
    startnodes = deque()
    for hap in samples:
        hc = -1
        for h in hapoblist:
            if h.hasgen(hap):
                if hc == -1:
                    hc = h.id
                    startnodes.append(h.id)
                    continue
                if (hc, h.id) not in addededge:
                    G.add_edge(hc, h.id)
                    addededge.add((hc, h.id))
                hc = h.id

    # Get threads
    # f = '/home/ra98jam/d16/projects/potato_hap_example/results/threading/EMv04.WhR.N50-2.csv'
    threads = readthreads(thfin)
    threads = sorted(threads, key=lambda x: x[0])
    tnodes = sorted(set([t1 for t in threads for t1 in t[1]]))
    tedges = set([(z[0], z[1]) for t in threads for z in zip(t[1][:-1], t[1][1:])])
    nploidy = Counter([t1 for t in threads for t1 in t[1]])


    # Add ploidy levels in the graph
    # ploidy = deque()
    # emout = pd.read_csv(emoutfin, header=None, sep='\t')

    def get_ploidy(row):
        if all([r == 0 for r in row[9:]]):
            return 0
        # TODO: for a node, I assign based on the maximum probablity from EM.
        #  So, if for a node, EM outputs probablity of 51% for 1 copy and 49%
        #  for 2 copy, I assign ploidy of 1 to the node. Probably, this can be
        #  fine-tuned.
        return row[9:].to_list().index(max(row[9:]))

    # END

    # ploidy = emout.apply(get_ploidy, axis=1)
    G.vs['ploidy'] = 0
    ks = list(nploidy.keys())
    G.vs[ks]['ploidy'] = [nploidy[k] for k in ks]

    G.es['bad_edge'] = [G.vs[e[0]]['ploidy'] != G.vs[e[1]]['ploidy'] for e in G.get_edgelist()]

    """ TO BE ASSESSED

    # Delete edges between nodes with unequal ploidy
    to_del = [e for e in G.get_edgelist() if G.vs[e[0]]['ploidy'] != G.vs[e[1]]['ploidy']]
    G.delete_edges(to_del)

    # Delete vertices with ploidy 0
    G.delete_vertices([i for i, p in enumerate(ploidy) if p == 0])

    # Check whether the predicted ploidy is correct
    ploidy_color = {
        0: "lightgrey",
        1: "orange",
        2: "blue",
        3: "red",
        4: "black"
    }
    cor_plo = deque()  # Green for correct ploidy, white for wrong
    for v in G.vs:
        p = v['ploidy']
        gs = vars(hapoblist[v['name']])['genomes']
        if p == len([g for g in gs if 'A_' in g]):
            cor_plo.append(ploidy_color[v['ploidy']])
        else:
            cor_plo.append(ploidy_color[0])
    G.vs['cor'] = cor_plo

    # Save figures
    ig.plot(G, target='tmp.pdf', layout='kk', vertex_size=np.array(G.degree()) + 4, vertex_color=G.vs['cor'],
            vertex_frame_color="black", vertex_frame_width=0, edge_width=0.1, edge_arrow_size=0.2, bbox=[800, 800])
    
    """

    # G.vs['dist'] = -1
    G.vs['height'] = -1
    G.vs['threadid'] = -1
    G.vs['start'] = [h.start for h in hapoblist]
    # Add thread ID for each node
    for i, t in enumerate(threads):
        for n in t[1]:
            if G.vs[n]['threadid'] == -1:
                G.vs[n]['threadid'] = i

    ploidy_color = {
        0: "lightgrey",
        1: "orange",
        2: "blue",
        3: "red",
        4: "black"
    }

    # <editor-fold desc='Plot variant 1'>
    def plot_graph(hapoblist, G, ax):

        segs = deque()
        cs = deque()
        lwd = deque()
        for i, h in enumerate(list(hapoblist)):
            segs.append((((h.start + 10000), G.vs[i]['height']), ((h.end - 10000), G.vs[i]['height'])))
            cs.append(ploidy_color[G.vs[i]['ploidy']])
            lwd.append(0.4 if G.vs[i]['ploidy'] < 1 else 3)
            # if 'dm' in h.genomes:
            #     cs.append('lightgrey')
            # else:
            #     cs.append('black')

        line_segments = LineCollection(segs, colors=cs, linestyle='solid', linewidths=lwd)
        ax.add_collection(line_segments)
        # ax.set_ylim(math.floor(min(G.vs['height'])), math.ceil(max(G.vs['height'])))
        # hpos = [G.vs['height'][h.id] for start in range(s, e, 100000) for h in hapoblist if h.start == start + 1]
        # ax.set_ylim(math.floor(min(hpos)), math.ceil(max(hpos)))
        ax.set_ylim(math.floor(min(G.vs['height'])), math.ceil(max(G.vs['height'])))
        # starts = list(range(34800000, 35400000, 100000))

        # Link haplotype blocks
        offset = 10000
        # for start in range(s, e, 100000):
        for start in sorted(set(G.vs['start'])):
            chaps = [(h.id, h.genomes) for h in hapoblist if h.start == start]
            pcoll = deque()
            pcollcol = deque()
            pcolllwd = deque()
            pcollzord = deque()
            # kcoll = deque()
            # print(chaps)
            for c in chaps:
                yh = G.vs[c[0]]['height']
                # Get links to downstream nodes
                ss = G.successors(c[0])
                for ssi in ss:
                    p1 = (hapoblist[c[0]].end - offset, yh)
                    p2 = (hapoblist[ssi].start + offset, G.vs[ssi]['height'])
                    z = 0 if (c[0], ssi) not in tedges else 1
                    pcoll.append(bezierlink(p1, p2, z))
                    if (c[0], ssi) not in tedges:
                    # if (G.vs[c[0]]['threadid'] != G.vs[ssi]['threadid']) or (G.vs[c[0]]['threadid'] == -1)  or (G.vs[ssi]['threadid'] == -1):
                        pcollcol.append('lightgrey')
                        pcolllwd.append(0.4)
                        pcollzord.append(0)
                    else:
                        # pcollcol.append(ploidy_color[G.vs[c[0]]['ploidy']])
                        pcollcol.append(ploidy_color[1])
                        pcolllwd.append(1.5)
                        pcollzord.append(1)

                # # Get kmer-proportion
                # nodestat = nkstats.loc[nkstats[0] == c[0]]
                # # print(start, c[0], nodestat)
                # if nodestat.shape[0] == 0:
                #     continue
                # maxkprop = ((nodestat[1] / maxx) * 80000).iloc[0]
                # # print(start, c[0], maxkprop)
                # segments = np.array([[[hapoblist[c[0]].start + offset, yh+0.075], [hapoblist[c[0]].start + offset + maxkprop, yh + 0.075]],
                #                      [[hapoblist[c[0]].start + offset + maxkprop, yh + 0.075], [hapoblist[c[0]].end - offset, yh + 0.075]]])
                # lc = LineCollection(segments, colors=['blue', 'lightblue'], linewidth=2)
                # ax.add_collection(lc)
                #
                # kseqprop = (nodestat[3] * 80000).iloc[0]
                # # print(start, c[0], kseqprop)
                # segments = np.array([[[hapoblist[c[0]].start + offset, yh+0.15], [hapoblist[c[0]].start + offset + kseqprop, yh + 0.15]],
                #                      [[hapoblist[c[0]].start + offset + kseqprop, yh + 0.15], [hapoblist[c[0]].end - offset, yh + 0.15]]])
                # lc = LineCollection(segments, colors=['red', 'pink'], linewidth=2)
                # logger.debug('Plotting line')
                # ax.add_collection(lc)
                # ax.hlines(yh+np.array([0.25, 0.35, 0.45, 0.55, 0.65, 0.75]), xmin=start+offset, xmax=start+100000-offset, linestyles='dashed', colors='grey', linewidth=0.2)
                #
                # if outmode == 1:
                #     # Style 1
                #     try:
                #         xs = np.linspace(start+offset, start+100000-offset, nodestat[2].iloc[0])
                #         ys = sorted((np.array([int(x) for x in nodestat[5].iloc[0].split(',')])/hcnt)*0.1)
                #         ax.plot(xs, yh+0.25+ys, linewidth=1, color='black')
                #     except AttributeError:
                #         continue
                # elif outmode == 2:
                #     # Style 2
                #     try:
                #         kcnts = Counter([int(x) for x in nodestat[5].iloc[0].split(',')])
                #         mk = max(kcnts.values())
                #         xs = np.linspace(start+offset, start+100000-offset, 161)
                #         ax.plot(xs, yh+0.25+np.array([(kcnts[i]/mk)*0.5 for i in range(len(xs))]), color='black', linewidth=0.2)
                #         ax.vlines(xs[[(int(hcnt)*i) for i in range(1, 6)]], yh+0.25, yh+0.75, colors='grey', zorder=0, linewidth=0.2)
                #     except AttributeError:
                #         continue

            # p = PatchCollection(list(pcoll), facecolors='none', edgecolors='grey', linewidths=0.2)
            p = PatchCollection(list(pcoll), facecolors='none', edgecolors=pcollcol, linewidths=pcolllwd)
            ax.add_collection(p)
        return ax
    # END
    #
    #
    # # <editor-fold desc="Backup code">
    #
    # segs = deque()
    # cs = deque()
    # lwd = deque()
    # for i, h in enumerate(list(hapoblist)):
    #     segs.append((((h.start + 10000), G.vs[i]['height']), ((h.end - 10000), G.vs[i]['height'])))
    #     cs.append(ploidy_color[G.vs[i]['ploidy']])
    #     lwd.append(0.4 if G.vs[i]['ploidy'] < 1 else 2)
    #     # if 'dm' in h.genomes:
    #     #     cs.append('lightgrey')
    #     # else:
    #     #     cs.append('black')
    #
    # line_segments = LineCollection(segs, colors=cs, linestyle='solid', linewidths=lwd)
    # ax.add_collection(line_segments)
    # # ax.set_ylim(math.floor(min(G.vs['height'])), math.ceil(max(G.vs['height'])))
    # hpos = [G.vs['height'][h.id] for start in range(s, e, 100000) for h in hapoblist if h.start == start + 1]
    # ax.set_ylim(math.floor(min(hpos)), math.ceil(max(hpos)))
    #
    # # starts = list(range(34800000, 35400000, 100000))
    #
    # # Link haplotype blocks
    # offset = 10000
    # for start in range(s, e, 100000):
    #     chaps = [(h.id, h.genomes) for h in hapoblist if h.start == start + 1]
    #     pcoll = deque()
    #     pcollcol = deque()
    #     pcolllwd = deque()
    #     # kcoll = deque()
    #     # print(chaps)
    #     for c in chaps:
    #         yh = G.vs[c[0]]['height']
    #         # Get links to downstream nodes
    #         ss = G.successors(c[0])
    #         for ssi in ss:
    #             p1 = (hapoblist[c[0]].end - offset, yh)
    #             p2 = (hapoblist[ssi].start + offset, G.vs[ssi]['height'])
    #             pcoll.append(bezierlink(p1, p2))
    #             if (G.vs[c[0]]['threadid'] != G.vs[ssi]['threadid']) or (G.vs[ssi]['ploidy'] == 0):
    #                 pcollcol.append('lightgrey')
    #                 pcolllwd.append(0.4)
    #             else:
    #                 pcollcol.append(ploidy_color[G.vs[c[0]]['ploidy']])
    #                 pcolllwd.append(2)
    #
    #         # # Get kmer-proportion
    #         # nodestat = nkstats.loc[nkstats[0] == c[0]]
    #         # # print(start, c[0], nodestat)
    #         # if nodestat.shape[0] == 0:
    #         #     continue
    #         # maxkprop = ((nodestat[1] / maxx) * 80000).iloc[0]
    #         # # print(start, c[0], maxkprop)
    #         # segments = np.array([[[hapoblist[c[0]].start + offset, yh+0.075], [hapoblist[c[0]].start + offset + maxkprop, yh + 0.075]],
    #         #                      [[hapoblist[c[0]].start + offset + maxkprop, yh + 0.075], [hapoblist[c[0]].end - offset, yh + 0.075]]])
    #         # lc = LineCollection(segments, colors=['blue', 'lightblue'], linewidth=2)
    #         # ax.add_collection(lc)
    #         #
    #         # kseqprop = (nodestat[3] * 80000).iloc[0]
    #         # # print(start, c[0], kseqprop)
    #         # segments = np.array([[[hapoblist[c[0]].start + offset, yh+0.15], [hapoblist[c[0]].start + offset + kseqprop, yh + 0.15]],
    #         #                      [[hapoblist[c[0]].start + offset + kseqprop, yh + 0.15], [hapoblist[c[0]].end - offset, yh + 0.15]]])
    #         # lc = LineCollection(segments, colors=['red', 'pink'], linewidth=2)
    #         # logger.debug('Plotting line')
    #         # ax.add_collection(lc)
    #         # ax.hlines(yh+np.array([0.25, 0.35, 0.45, 0.55, 0.65, 0.75]), xmin=start+offset, xmax=start+100000-offset, linestyles='dashed', colors='grey', linewidth=0.2)
    #         #
    #         # if outmode == 1:
    #         #     # Style 1
    #         #     try:
    #         #         xs = np.linspace(start+offset, start+100000-offset, nodestat[2].iloc[0])
    #         #         ys = sorted((np.array([int(x) for x in nodestat[5].iloc[0].split(',')])/hcnt)*0.1)
    #         #         ax.plot(xs, yh+0.25+ys, linewidth=1, color='black')
    #         #     except AttributeError:
    #         #         continue
    #         # elif outmode == 2:
    #         #     # Style 2
    #         #     try:
    #         #         kcnts = Counter([int(x) for x in nodestat[5].iloc[0].split(',')])
    #         #         mk = max(kcnts.values())
    #         #         xs = np.linspace(start+offset, start+100000-offset, 161)
    #         #         ax.plot(xs, yh+0.25+np.array([(kcnts[i]/mk)*0.5 for i in range(len(xs))]), color='black', linewidth=0.2)
    #         #         ax.vlines(xs[[(int(hcnt)*i) for i in range(1, 6)]], yh+0.25, yh+0.75, colors='grey', zorder=0, linewidth=0.2)
    #         #     except AttributeError:
    #         #         continue
    #
    #     # p = PatchCollection(list(pcoll), facecolors='none', edgecolors='grey', linewidths=0.2)
    #     p = PatchCollection(list(pcoll), facecolors='none', edgecolors=pcollcol, linewidths=pcolllwd)
    #     ax.add_collection(p)
    # # </editor-fold>

    # </editor-fold>



    # <editor-fold desc='Process Graph'>
    win2node = defaultdict(deque)
    node2win = {}
    for h in hapoblist:
        win2node[int(h.start/100000)].append(h.id)
        node2win[h.id] = int(h.start/100000)

    G.vs['height'] = -10000
    tpos = set(list(range(-(spread-1), spread)))
    for i, thread in enumerate(threads):
        tnodes = thread[1]
        twin = [node2win[t] for t in tnodes]
        theights = set([G.vs[t]['height'] for w in twin for t in win2node[w]])
        th = np.random.choice(list(tpos - theights), 1)[0]
        # print([th if G.vs[t]['height'] == -10000 else G.vs[t]['height'] for t in tnodes])
        G.vs[tnodes]['height'] = [th if G.vs[t]['height'] == -10000 else G.vs[t]['height'] for t in tnodes] # Assign height to nodes for which the height has not been set by a longer thread

    # roots = [v.index for v in G.vs if v['dist'] == 0]
    roots = [v.index for v in G.vs if v['start'] == 1]
    # possible heights
    pheights = sorted(set(np.ceil(np.arange(len(roots)) - (len(roots) / 2))) - set(G.vs[roots]['height']))
    for root in roots:
        if G.vs[root]['height'] == -10000:
            G.vs[root]['height'] = pheights.pop(0)

    # G.vs[roots]['height'] = np.ceil(np.arange(len(roots)) - (len(roots) / 2))
    # maxx = max(G.vs['dist'])
    # for i in range(1, maxx):
    for i in sorted(set(G.vs['start'])):
        ns = [v.index for v in G.vs if v['start'] == i]
        heights = sorted(set(np.ceil(np.arange(len(ns)) - (len(ns) / 2))) - set([G.vs[v]['height'] for v in ns]))
        hord = sorted(range(len(heights)), key=lambda x: np.mean(G.vs[G.predecessors(ns[x])]['height']))
        nh = deque()
        cnt = 0
        for j in ns:
            if G.vs[j]['height'] == -10000:
                nh.append(heights[hord.index(cnt)])
                cnt += 1
            else:
                nh.append(G.vs[j]['height'])
        G.vs[ns]['height'] = list(nh)
    # </editor-fold>

    # <editor-fold desc="Get Plot">

    fig, ax = plt.subplots(figsize=[W, H])
    ax = plot_graph(hapoblist, G, ax)
    ax = cleanax(ax)
    maxl = max(G.vs['start'])+100000
    ax.set_xlim(0, maxl)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_visible(False)

    ax.ticklabel_format(axis='x', useOffset=False, style='plain')
    xticks = ax.get_xticks()
    if maxl >= 1000000000:
        xticksl = xticks / 1000000000
        ax.set_xlabel('Chromosome position (in Gbp)')
    elif maxl >= 1000000:
        xticksl = xticks / 1000000
        ax.set_xlabel('Chromosome position (in Mbp)')
    elif maxl >= 1000:
        xticksl = xticks / 1000
        ax.set_xlabel('Chromosome position (in Kbp)')
    else:
        xticksl = xticks
        ax.set_xlabel('Chromosome position')
    ax.set_xticks(xticks[:-1])
    ax.set_xticklabels(xticksl[:-1])
    plt.tight_layout()
    # </editor-fold>

    plt.savefig(outfin, dpi=300)
    logger.info(f'Finished plotting {outfin}')
    logger.handlers.clear()
    plt.close()
    return
# END


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot haplotype graph and K-mer stats for each node in the graph.')
    parser.add_argument('hap', help='Haplotype graph', type=argparse.FileType('r'))
    parser.add_argument('threads', help='Threads from EM', type=argparse.FileType('r'))
    # parser.add_argument('nkstats', help='Kmer stats for each node or window', type=argparse.FileType('r'))
    parser.add_argument('out', help='Output plot file name (add .pdf or .png)', type=argparse.FileType('w'))
    parser.add_argument('--spread', dest='spread', help='Regulates vertical spread of the selected nodes (minimum allowed value 3)', type=int, default=3)
    parser.add_argument('-W', help='Plot width', type=int, default=8)
    parser.add_argument('-H', help='Plot height', type=int, default=6)
    # parser.add_argument('-s', dest='start', help='Start position for the visualisation', type=int, default=0)
    # parser.add_argument('-e', dest='end', help='End position for the visualisation', type=int, default=1000001)
    # parser.add_argument('--mode', dest='mode', help='Kmer-count mode', type=int, default=1, choices=[1, 2])
    args = parser.parse_args()
    plot_haplotype_graph(args)
# END
