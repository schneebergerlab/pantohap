import pandas as pd
import os
from collections import defaultdict, deque, Counter, OrderedDict
from matplotlib import pyplot as plt

import socket

if socket.gethostname() == 'LMBIDP1-LXSCH01':
    # Indir for local compute
    indir = '/home/ra98jam/d16/projects/potato_hap_example/data/'
    pwd = '/home/ra98jam/d16/projects/potato_hap_example/results/threading/thread_fastas/'
else:
    # Indir for cluster compute
    indir = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/data/'
    pwd = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/potato_hap_example/results/threading/thread_fastas/'

os.chdir(pwd)
# <editor-fold desc="Reader">

def samtocoords(f):
    import sys
    import logging
    from pandas import DataFrame
    from collections import deque
    logger = logging.getLogger('SAM reader')
    rc = {}  # Referece chromosomes
    rcs = {}  # Selected chromosomes
    al = deque()  # Individual alignment
    try:
        with open(f, 'r') as fin:
            for l in fin:
                if l[:3] == '@SQ':
                    c, s = 0, 0
                    for h in l.strip().split()[1:]:
                        h = h.split(':')
                        if h[0] == 'SN': c = h[1]
                        if h[0] == 'LN': s = int(h[1])
                    rcs[c] = s
                    continue
                elif l[0] == '@':
                    continue

                l = l.split('\t')[:6]
                # if l[1] == '2064': break
                if l[2] == '*':
                    logger.warning(l[
                                       0] + ' do not align with any reference sequence and cannot be analysed. Remove all unplaced scaffolds and contigs from the assemblies.')  # Skip rows corresponding to non-mapping sequences (contigs/scaffolds)
                    continue

                if 'M' in l[5]:
                    logger.error(
                        'Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: ' + l[5])
                    sys.exit()
                cgt = [[int(j[0]), j[1]] for j in [i.split(';') for i in
                                                   l[5].replace('S', ';S,').replace('H', ';H,').replace('=',
                                                                                                        ';=,').replace(
                                                       'X', ';X,').replace('I', ';I,').replace('D', ';D,').split(',')[
                                                   :-1]]]
                if len(cgt) > 2:
                    if True in [True if i[1] in ['S', 'H'] else False for i in cgt[1:-1]]:
                        logger.error(
                            f"Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: {cgt}")
                        sys.exit()

                bf = '{:012b}'.format(int(l[1]))

                rs = int(l[3])
                re = rs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'D']])

                if bf[7] == '0':  # forward alignment
                    if cgt[0][1] == '=':
                        qs = 1
                    elif cgt[0][1] in ['S', 'H']:
                        qs = cgt[0][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                elif bf[7] == '1':  # inverted alignment
                    if cgt[-1][1] == '=':
                        qs = 1
                    elif cgt[-1][1] in ['S', 'H']:
                        qs = cgt[-1][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                    qs, qe = qe, qs

                al.append([
                    rs,
                    re,
                    qs,
                    qe,
                    abs(re - rs) + 1,
                    abs(qs - qe) + 1,
                    format((sum([i[0] for i in cgt if i[1] == '=']) / sum(
                        [i[0] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])) * 100, '.2f'),
                    1,
                    1 if bf[7] == '0' else -1,
                    l[2],
                    l[0],
                    "".join([str(i[0]) + i[1] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])
                ])
                rcs[l[2]] = 1
            rcs = list(rcs.keys())
            for k in list(rc.keys()):
                if k not in rcs: logger.warning(l[
                                                    0] + ' do not align with any query sequence and cannot be analysed. Remove all unplaced scaffolds and contigs from the assemblies.')
    except Exception as e:
        logger.error('Error in reading SAM file: ' + str(e))
        sys.exit()
    al = DataFrame(list(al))
    al[6] = al[6].astype('float')
    al.sort_values([9, 0, 1, 2, 3, 10], inplace=True, ascending=True)
    al.index = range(len(al.index))
    return al
# END


def readsambam(fin, ftype='B'):
    import pysam
    import logging
    import sys
    import numpy as np
    import pandas as pd
    logger = logging.getLogger('Reading BAM/SAM file')
    try:
        if ftype == 'B':
            findata = pysam.AlignmentFile(fin, 'rb')
        elif ftype == 'S':
            return samtocoords(fin)
        else:
            raise ValueError("Wrong parameter")
    except ValueError as e:
        logger.error("Error in opening BAM/SAM file. " + str(e))
        sys.exit()
    except OSError as e:
        logger.error("Error in reading input file." + str(e))
        sys.exit()
    except Exception as e:
        logger.error("Unexpected error in opening BAM/SAM file. " + str(e))
        sys.exit()

    try:
        qry_prim = {}
        ref_prim = {}
        cgdict = {1: 'I', 2: 'D', 7: '=', 8: 'X'}
        coords = {}
        index = 0
        for aln in findata:
            index += 1
            ## Check whether every sequence has at least one primary alignment
            if aln.reference_name is not None:
                if aln.reference_name not in ref_prim.keys():
                    ref_prim[aln.reference_name] = False
            if aln.query_name not in qry_prim.keys():
                qry_prim[aln.query_name] = False
            if aln.reference_name is not None:
                if not ref_prim[aln.reference_name]:
                    if aln.flag < 256:
                        ref_prim[aln.reference_name] = True
            if not qry_prim[aln.query_name]:
                if aln.flag < 256:
                    qry_prim[aln.query_name] = True

            ## Pass non-alinging chromosomes
            if aln.cigarstring is None:
                logger.warning(aln.query_name + ' do not align with any reference chromosome and cannot be analysed')
                continue

            ## Check CIGAR:
            if False in [False if i[0] not in [1, 2, 4, 5, 7, 8] else True for i in aln.cigartuples]:
                logger.error(
                    "Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: " + str(
                        aln.cigarstring))
                sys.exit()
            if len(aln.cigartuples) > 2:
                if True in [True if i[0] in [4, 5] else False for i in aln.cigartuples[1:-1]]:
                    logger.error(
                        "Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + aln.cigarstring)
                    sys.exit()

            ## Parse information from the aln object
            astart = aln.reference_start + 1
            aend = aln.reference_end
            is_inv = True if np.binary_repr(aln.flag, 12)[7] == '1' else False
            if not is_inv:
                if aln.cigartuples[0][0] in [4, 5]:
                    bstart = aln.cigartuples[0][1] + 1
                else:
                    bstart = 1
                bend = bstart + aln.query_alignment_length - 1
            else:
                if aln.cigartuples[-1][0] in [4, 5]:
                    bend = aln.cigartuples[-1][1] + 1
                else:
                    bend = 1
                bstart = bend + aln.query_alignment_length - 1
            alen = abs(aend - astart) + 1
            blen = abs(bend - bstart) + 1
            iden = format((sum([i[1] for i in aln.cigartuples if i[0] == 7]) / sum(
                [i[1] for i in aln.cigartuples if i[0] in [1, 2, 7, 8]])) * 100, '.2f')
            adir = 1
            bdir = -1 if is_inv else 1
            achr = aln.reference_name
            bchr = aln.query_name
            cg = "".join([str(i[1]) + cgdict[i[0]] for i in aln.cigartuples if i[0] not in [4, 5]])
            coords[index] = [astart, aend, bstart, bend, alen, blen, iden, adir, bdir, achr, bchr, cg]

        ## Give warning for chromosomes which do not have any primary alignment
        for k, v in ref_prim.items():
            if not v:
                logger.warning(
                    'No primary alignment found for reference sequence ' + k + '. This could mean that the entire chromosome ' + k + ' is reapeated.')
        for k, v in qry_prim.items():
            if not v:
                logger.warning(
                    'No primary alignment found for query sequence ' + k + '. This could mean that the entire chromosome ' + k + ' is reapeated.')

        ## Return alignments
        coords = pd.DataFrame.from_dict(coords, orient='index')
        coords.sort_values([9, 0, 1, 2, 3, 10], inplace=True, ascending=True)
        coords.index = range(len(coords.index))
        coords[6] = coords[6].astype('float')
        return coords
    except Exception as e:
        logger.error("Error in reading BAM/SAM file. " + str(e))
        sys.exit()
# END


def readpaf(paf):
    import logging
    import sys
    from collections import deque
    import pandas as pd
    coords = deque()
    logger = logging.getLogger('Reading BAM/SAM file')
    try:
        with open(paf, 'r') as fin:
            for line in fin:
                line = line.strip().split()
                astart = int(line[7]) + 1
                aend = int(line[8])
                adir = 1
                bdir = 1 if line[4] == '+' else -1
                bstart = int(line[2]) + 1 if bdir == 1 else int(line[3])
                bend = int(line[3]) if bdir == 1 else int(line[2]) + 1
                alen = abs(aend - astart) + 1
                blen = abs(bend - bstart) + 1 if bdir == 1 else bstart - bend + 1
                cg = [i.split(":")[-1] for i in line[12:] if i[:2] == 'cg']
                if len(cg) != 1:
                    logger.error("CIGAR string is not present in PAF at line {}. Exiting.".format("\t".join(line)))
                    sys.exit()
                cg = cg[0]
                ## Check CIGAR:
                if not all([True if i[1] in {'I', 'D', 'H', 'S', 'X', '='} else False for i in cgtpl(cg)]):
                    logger.error(
                        "Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: " + str(
                            cg))
                    sys.exit()
                if len(cgtpl(cg)) > 2:
                    if any([True if i[1] in {'H', 'S'} else False for i in cgtpl(cg)]):
                        logger.error(
                            "Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + str(
                                cg))
                        sys.exit()

                iden = round((sum([int(i[0]) for i in cgtpl(cg) if i[1] == '=']) / sum(
                    [int(i[0]) for i in cgtpl(cg) if i[1] in {'=', 'X', 'D', 'I'}])) * 100, 2)
                achr = line[5]
                bchr = line[0]
                coords.append([astart, aend, bstart, bend, alen, blen, iden, adir, bdir, achr, bchr, cg])
        coords = pd.DataFrame(coords)
        coords.sort_values([9, 0, 1, 2, 3, 10], inplace=True, ascending=True)
        coords.index = range(len(coords.index))
        coords[6] = coords[6].astype('float')
        return coords
    except FileNotFoundError:
        logger.error("Cannot open {} file. Exiting".format(paf))
        sys.exit()
    except ValueError as e:
        logger.error("Error in reading PAF: {}. Exiting".format(e))
        sys.exit()
# END


def readcoords(coordsfin, ftype, f, cigar=False):
    import logging
    import pandas as pd
    import numpy as np
    import sys
    logger = logging.getLogger('Reading Coords')
    logger.debug(ftype)
    if ftype == 'T':
        logger.info("Reading input from .tsv file")
        try:
            coords = pd.read_table(coordsfin, header=None)
        except pd.errors.ParserError:
            coords = pd.read_table(coordsfin, header=None, engine="python")
        except Exception as e:
            logger.error("Error in reading the alignment file. " + str(e))
            sys.exit()
    elif ftype == 'S':
        logger.info("Reading input from SAM file")
        try:
            coords = readsambam(coordsfin, ftype='S')
        except Exception as e:
            logger.error("Error in reading the alignment file. " + str(e))
            sys.exit()
    elif ftype == 'B':
        logger.info("Reading input from BAM file")
        try:
            coords = readsambam(coordsfin, ftype='B')
        except Exception as e:
            logger.error("Error in reading the alignment file" + str(e))
            sys.exit()
    elif ftype == 'P':
        logger.info("Reading input from PAF file")
        try:
            coords = readpaf(coordsfin)
        except Exception as e:
            logger.error("Error in reading the alignment file" + str(e))
            sys.exit()
    else:
        logger.error("Incorrect alignment file type specified.")
        sys.exit()

    if not cigar:
        if coords.shape[1] >= 12:
            coords = coords.iloc[:, 0:11]
        coords.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr"]
    else:
        if coords.shape[1] > 12:
            coords = coords.iloc[:, 0:12]
        coords.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr",
                          'cigar']

    # Sanity check input file
    try:
        coords.aStart = coords.aStart.astype('int')
    except ValueError:
        logger.error('astart is not int')
        sys.exit()

    try:
        coords.aEnd = coords.aEnd.astype('int')
    except ValueError:
        logger.error('aend is not int')
        sys.exit()

    try:
        coords.bStart = coords.bStart.astype('int')
    except ValueError:
        logger.error('bstart is not int')
        sys.exit()

    try:
        coords.bEnd = coords.bEnd.astype('int')
    except ValueError:
        logger.error('abend is not int')
        sys.exit()

    try:
        coords.aLen = coords.aLen.astype('int')
    except ValueError:
        logger.error('alen is not int')
        sys.exit()

    try:
        coords.bLen = coords.bLen.astype('int')
    except ValueError:
        logger.error('blen is not int')
        sys.exit()

    try:
        coords.iden = coords.iden.astype('float')
    except ValueError:
        logger.error('iden is not float')
        sys.exit()

    try:
        coords.aDir = coords.aDir.astype('int')
    except ValueError:
        logger.error('aDir is not int')
        sys.exit()

    if any(coords.aDir != 1):
        logger.error('aDir can only have values 1')
        sys.exit()

    try:
        coords.bDir = coords.bDir.astype('int')
    except ValueError:
        logger.error('bDir is not int')
        sys.exit()

    for i in coords.bDir:
        if i not in [1, -1]:
            logger.error('bDir can only have values 1/-1')
            sys.exit()

    try:
        coords.aChr = coords.aChr.astype(str)
    except:
        logger.error('aChr is not string')
        sys.exit()

    try:
        coords.bChr = coords.bChr.astype(str)
    except:
        logger.error('bChr is not string')
        sys.exit()

    # Filter small alignments
    if f:
        coords = coords.loc[coords.iden > 90]
        coords = coords.loc[(coords.aLen > 100) & (coords.bLen > 100)]

    ## check for bstart > bend when bdir is -1
    check = np.unique(coords.loc[coords.bDir == -1, 'bStart'] > coords.loc[coords.bDir == -1, 'bEnd'])
    if len(check) > 1:
        logger.error(
            'Inconsistent start and end position for inverted alignment in query genome. For inverted alignments, either all bstart < bend or all bend > bstart')
        sys.exit()
    elif len(check) == 0:
        logger.info('No Inverted alignments present.')
    elif check[0]:
        pass
    else:
        logger.info('For inverted alignments, bstart was less than bend. Swapping them.')
        coords.loc[coords.bDir == -1, 'bStart'] = coords.loc[coords.bDir == -1, 'bStart'] + coords.loc[
            coords.bDir == -1, 'bEnd']
        coords.loc[coords.bDir == -1, 'bEnd'] = coords.loc[coords.bDir == -1, 'bStart'] - coords.loc[
            coords.bDir == -1, 'bEnd']
        coords.loc[coords.bDir == -1, 'bStart'] = coords.loc[coords.bDir == -1, 'bStart'] - coords.loc[
            coords.bDir == -1, 'bEnd']
    coords.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)
    return coords
# END

inttr = lambda x: [int(x[0]), x[1]]
def cgtpl(cg, to_int=False):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    for i in "MIDNSHPX=":
        cg = cg.replace(i, ';'+i+',')
    if to_int:
        return [inttr(i.split(';')) for i in cg.split(',')[:-1]]
    else:
        return [i.split(';') for i in cg.split(',')[:-1]]
# END

# </editor-fold>

def get_chrs_len_from_fai(fin):
    clen = {}
    with open(fin, 'r') as f:
        for line in f:
            line = line.strip().split()
            clen[line[0]] = int(line[1])
    return clen


# Using the sorted contigs list generated using the D-genesis webtool
sort_pos = pd.read_table("RussetBurbank.pseudo_contigs_RussetBurbank.genome_assoc.tsv")
# Remove unplaced scaffolds
sort_pos = sort_pos.loc[['chr' in i for i in sort_pos.Target]]
selected_pseudo_contigs = sort_pos.Query.unique()
sort_pos_grp = sort_pos.groupby('Target')
coords = readcoords('rb.contigs_to_genome.paf', 'P', f=False, cigar=True)
coords = coords.loc[['chr' in i for i in coords.aChr]]
coords = coords.loc[coords.bChr.isin(selected_pseudo_contigs)]
rchrs_len = get_chrs_len_from_fai('RussetBurbank.genome.fa.fai')
qchrs_len = get_chrs_len_from_fai('RussetBurbank.pseudo_contigs.fa.fai')


al = coords.copy()
al['aStart aEnd bStart bEnd'.split()] = al['aStart aEnd bStart bEnd'.split()].astype('float')
# ragpdata = {}
# qagpdata = {}
# if ragp is not None: ragpdata = readagp(ragp.name)
# if qagp is not None: qagpdata = readagp(qagp.name)
# TODO: Ensure that there are no errors when there are no alignments for a sequence
# alachr = set(np.unique(al.aChr))
# albchr = set(np.unique(al.bChr))
# rchrs = sorted([k for k in rchrs_len.keys() if k in alachr], key=lambda x: rchrs_len[x], reverse=True)
# qchrs = sorted([k for k in qchrs_len.keys() if k in albchr], key=lambda x: qchrs_len[x], reverse=True)
# qchrs = sorted(qchrs_len.keys(), key=lambda x: qchrs_len[x], reverse=True)
# rchrs = []
# [rchrs.append(x) for x in sort_pos.Target.values if x not in rchrs]
rchrs = [r for r in sorted(rchrs_len) if 'chr' in r]
qchrs = []
_ = [qchrs.append(x) for r in rchrs for x in sort_pos_grp.get_group(r).Query.values if x not in qchrs]
rcumsum = deque([0])
qcumsum = deque([0])
roffdict = OrderedDict()
qoffdict = OrderedDict()

for i, r in enumerate(rchrs):
    roffdict[r] = sum([rchrs_len[rchrs[j]] for j in range(i)])

for i, q in enumerate(qchrs):
    qoffdict[q] = sum([qchrs_len[qchrs[j]] for j in range(i)])

al['aStart'] += [roffdict[c] for c in al['aChr']]
al['aEnd'] += [roffdict[c] for c in al['aChr']]
al['bStart'] += [qoffdict[c] for c in al['bChr']]
al['bEnd'] += [qoffdict[c] for c in al['bChr']]

#
#     cumsum = sum([rchrs_len[rchrs[j]] for j in range(i)])
#     al.loc[al['aChr'] == rchrs[i], 'aStart'] += cumsum
#     al.loc[al['aChr'] == rchrs[i], 'aEnd'] += cumsum
#     # if ragp is not None:
#     #     for k in range(len(ragpdata[rchrs[i]])):
#     #         ragpdata[rchrs[i]][k][0] += cumsum
#     #         ragpdata[rchrs[i]][k][1] += cumsum
#     rcumsum.append(cumsum)
#     print(rchrs[i], al.aStart.min())
#
# for i in range(1, len(qchrs)):
#     cumsum = sum([qchrs_len[qchrs[j]] for j in range(i)])
#     al.loc[al['bChr'] == qchrs[i], 'bStart'] += cumsum
#     al.loc[al['bChr'] == qchrs[i], 'bEnd'] += cumsum
#     # if qagp is not None:
#     #     for k in range(len(qagpdata[qchrs[i]])):
#     #         qagpdata[qchrs[i]][k][0] += cumsum
#     #         qagpdata[qchrs[i]][k][1] += cumsum
#     qcumsum.append(cumsum)
#
al_data = deque()
for row in al.itertuples(index=False):
    # if row[4] < 1000 or row[5] < minsize: next
    al_data.append([row[0], row[1]])
    al_data.append([row[2], row[3]])
    if row[8] == 1: al_data.append('r')
    if row[8] == -1: al_data.append('b')
al_data = list(al_data)
xticks = [roffdict[r] + (rchrs_len[r]) / 2 for r in rchrs]
yticks = [qoffdict[q] + (qchrs_len[q]) / 2 for q in qchrs]
# yticks = [qcumsum[j] + (qchrs_len[qchrs[j]]) / 2 for j in range(len(qchrs))]

# Pericentromere coords
centro = pd.read_csv(f'{indir}/assemblies/R_48chrs_pericentromere.bed', header=None, sep='\t')
centro[0] = [c.rsplit('_', maxsplit=1)[0] for c in centro[0].astype('str')]
centro[1] += [roffdict[r] for r in centro[0]]
centro[2] += [roffdict[r] for r in centro[0]]

# logger.info('starting drawing')
width = 18
height = 12
figure = plt.figure(figsize=[width, height])
ax = plt.subplot(1, 1, 1)
ax.margins(x=0, y=0)
ax.set_xlim([0, sum([rchrs_len[k] for k in rchrs])])
ax.set_ylim([0, sum([qchrs_len[k] for k in qchrs])])
ax.plot(*al_data, linewidth=0.5, zorder=1)
ax.set_xticks(xticks)
# for r in rcumsum:
#     ax.axvline(r, linestyle='--', color='lightgrey', alpha=1, linewidth=0.1)

for r in np.array(list(roffdict.values())):
    ax.axvline(r, linestyle='-', color='grey', alpha=1, linewidth=0.5)
for r in np.array(list(roffdict.values()))[4::4]:
    ax.axvline(r, linestyle='-', color='black', alpha=1, linewidth=2)

for row in centro.itertuples(index=False):
    ax.axvspan(row[1], row[2], color='lightgrey', alpha=0.5, zorder=0, linewidth=0)


ax.set_xticklabels(rchrs, rotation=90)
ax.set_xticklabels(rchrs, rotation=90)
ax.tick_params(left=False, labelleft=False)

ax.set_xlabel("RB chromosomes")
ax.set_ylabel("Predicted pseudo-contigs")
plt.tight_layout()
plt.savefig('dotplot.png')
plt.savefig('dotplot.pdf')
plt.close()
