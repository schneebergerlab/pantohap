#!/usr/bin/env python3

import argparse

import pandas as pd


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


def bezierlink(p1, p2):
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
    patch = patches.PathPatch(path, facecolor="none")
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

    logger = mylogger(__name__)
    chrsize = 46102915 # Currently, set for Chr2
    logger.info(f"Setting chromosome size to {chrsize}")
    samples = ['dm'] + [f'{i}_hap{j}' for i, j in product(string.ascii_uppercase[:10], range(5, 9))]

    hapfin = args.hap.name
    nkfin = args.nkstats.name
    outfin = args.out.name

    s = args.start
    e = args.end
    # s = 34800000
    # e = 35400000

    # hapfin = '/home/ra98jam/d16/projects/potato_hap_example/data/haplotype_graph.txt'
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

    # to_delete_ids = [v.index for v in G.vs if v.degree() == 0]
    # G.delete_vertices(to_delete_ids)
    G.vs['dist'] = -1
    G.vs['height'] = -1
    G.vs['start'] = [h.start for h in hapoblist]

    # roots = [v.index for v in G.vs if v['dist'] == 0]
    roots = [v.index for v in G.vs if v['start'] == 1]
    G.vs[roots]['height'] = np.ceil(np.arange(len(roots)) - (len(roots) / 2))
    # maxx = max(G.vs['dist'])
    # for i in range(1, maxx):
    for i in sorted(set(G.vs['start']))[1:]:
        ns = [v.index for v in G.vs if v['start'] == i]
        heights = np.ceil(np.arange(len(ns)) - (len(ns) / 2))
        hord = sorted(range(len(ns)), key=lambda x: np.mean(G.vs[G.predecessors(ns[x])]['height']))
        G.vs[ns]['height'] = [heights[hord.index(j)] for j in range(len(ns))]

    G.vs['height'] += np.random.uniform(-0.1, 0.1, len(G.vs))

    # For each node draw a track for number of kmers
    ## read node_k_stats
    nkstats = pd.read_table(nkfin, header=None, sep='\t')
    # nkstats = pd.read_table(
    #     '/home/ra98jam/d16/projects/potato_hap_example/results/kmer_analysis/kmer_size_21/node_k_stats.txt',
    #     header=None, sep='\t')
    maxx = max(nkstats[1])
    nkstats.fillna(0, inplace=True)
    kmercnts = deque()
    for i in nkstats[5]:
        for j in str(i).split(','):
            kmercnts.appendleft(int(j))
    kmercnts = Counter(kmercnts)
    xs = sorted(list(kmercnts.keys()))
    ys = [kmercnts[x] for x in xs]
    pfit = P.fit(xs, ys, 16)

    # Get most frequent count
    fx, fy = pfit.linspace(100)  # generate 100 sample points on this graph
    hcnt = fx[np.where(fy == max(fy))[0][0]]
    logger.info(f"Selected haplotig kmer count: {hcnt}")

    fig, ax = plt.subplots(figsize=[8, 6])
    plt.tight_layout()
    ax = cleanax(ax)
    # ax.set_xlim(0, chrsize)
    ax.set_xlim(s, e)

    segs = deque()
    cs = deque()
    for i, h in enumerate(list(hapoblist)):
        segs.append((((h.start + 10000), G.vs[i]['height']), ((h.end - 10000), G.vs[i]['height'])))
        if 'dm' in h.genomes:
            cs.append('lightgrey')
        else:
            cs.append('black')

    line_segments = LineCollection(segs, colors=cs, linestyle='solid', linewidths=1)
    ax.add_collection(line_segments)
    ax.set_ylim(math.floor(min(G.vs['height'])), math.ceil(max(G.vs['height'])))

    # starts = list(range(34800000, 35400000, 100000))

    # Link haplotype blocks
    offset = 10000
    for start in range(s, e, 100000):
        chaps = [(h.id, h.genomes) for h in hapoblist if h.start == start + 1]
        pcoll = deque()
        kcoll = deque()
        # print(chaps)
        for c in chaps:
            yh = G.vs[c[0]]['height']
            # Get links to downstream nodes
            ss = G.successors(c[0])
            for ssi in ss:
                p1 = (hapoblist[c[0]].end - offset, yh)
                p2 = (hapoblist[ssi].start + offset, G.vs[ssi]['height'])
                pcoll.append(bezierlink(p1, p2))
            # Get kmer-proportion

            nodestat = nkstats.loc[nkstats[0] == c[0]]
            # print(start, c[0], nodestat)
            if nodestat.shape[0] == 0:
                continue
            maxkprop = ((nodestat[2] / maxx) * 80000).iloc[0]
            # print(start, c[0], maxkprop)
            segments = np.array([[[hapoblist[c[0]].start + offset, yh+0.075], [hapoblist[c[0]].start + offset + maxkprop, yh + 0.075]],
                                 [[hapoblist[c[0]].start + offset + maxkprop, yh + 0.075], [hapoblist[c[0]].end - offset, yh + 0.075]]])
            lc = LineCollection(segments, colors=['blue', 'lightblue'], linewidth=2)
            ax.add_collection(lc)

            kseqprop = (nodestat[3] * 80000).iloc[0]
            # print(start, c[0], kseqprop)
            segments = np.array([[[hapoblist[c[0]].start + offset, yh+0.15], [hapoblist[c[0]].start + offset + kseqprop, yh + 0.15]],
                                 [[hapoblist[c[0]].start + offset + kseqprop, yh + 0.15], [hapoblist[c[0]].end - offset, yh + 0.15]]])
            lc = LineCollection(segments, colors=['red', 'pink'], linewidth=2)
            logger.debug('Plotting line')
            ax.add_collection(lc)
            ax.hlines(yh+np.array([0.25, 0.35, 0.45, 0.55, 0.65, 0.75]), xmin=start+offset, xmax=start+100000-offset, linestyles='dashed', colors='grey', linewidth=0.2)

            # # Style 1
            # try:
            #     xs = np.linspace(start+offset, start+100000-offset, nodestat[2].iloc[0])
            #     ys = sorted((np.array([int(x) for x in nodestat[5].iloc[0].split(',')])/hcnt)*0.1)
            #     ax.plot(xs, yh+0.25+ys, linewidth=1, color='black')
            # except AttributeError:
            #     continue

            # Style 2
            try:
                kcnts = Counter([int(x) for x in nodestat[5].iloc[0].split(',')])
                mk = max(kcnts.values())
                xs = np.linspace(start+offset, start+100000-offset, 161)
                ax.plot(xs, yh+0.25+np.array([(kcnts[i]/mk)*0.5 for i in range(len(xs))]), color='black', linewidth=0.2)
                ax.vlines(xs[[(int(hcnt)*i) for i in range(1, 6)]], yh+0.25, yh+0.75, colors='grey', zorder=0, linewidth=0.2)
            except AttributeError:
                continue

        p = PatchCollection(list(pcoll), facecolors='none', edgecolors='grey', linewidths=0.2)
        ax.add_collection(p)
    plt.tight_layout()
    plt.savefig(outfin, dpi=150)
    logger.info(f'Finished plotting {outfin}')
    logger.handlers.clear()
    plt.close()
    return
# END


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot haplotype graph and K-mer stats for each node in the graph.')
    parser.add_argument('hap', help='Haplotype graph', type=argparse.FileType('r'))
    parser.add_argument('nkstats', help='Kmer stats for each node or window', type=argparse.FileType('r'))
    parser.add_argument('out', help='Output plot file name (add .pdf or .png)', type=argparse.FileType('w'))
    parser.add_argument('-s', dest='start', help='Start position for the visualisation', type=int, default=0)
    parser.add_argument('-e', dest='end', help='End position for the visualisation', type=int, default=1000001)
    args = parser.parse_args()
    plot_haplotype_graph(args)