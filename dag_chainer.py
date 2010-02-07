#!/usr/bin/env python
# -*- coding: UTF-8 -*-

CONSTANT_MATCH_SCORE=None
MAX_MATCH_SCORE=50.0

import math
from subprocess import Popen, PIPE
import os
op = os.path
PATH = op.dirname(op.abspath(__file__))
import sys
sys.path.insert(0, PATH)
from dagtools import DagLine
import collections
import operator

try:
    from processing import Process, Pipe as mPipe
except ImportError:
    from multiprocessing import Process, Pipe as mPipe

def scoringF(evalue, constant_match=CONSTANT_MATCH_SCORE, max_match=MAX_MATCH_SCORE):
    if not constant_match is None:
        return constant_match

    matchScore = 10 * -math.log10(evalue);
    matchScore = int(matchScore +.5) / 10
    return max_match if matchScore > max_match else matchScore


def get_dag_line(fh):
    line = fh.readline()
    if not line: return None, None
    if line[0] == "#": 
        get_dag_line.header = parse_pyheader(line, asstring=True)
        line = fh.readline()
    return DagLine(line), get_dag_line.header
# `get_dag_line.header` is sort of a hack to save the last seen header.
# so we know every dagline that follows belong to this group.
get_dag_line.header = None

JS="^" # magic separator. (J)oin (S)tring

def parse_pyheader(header, asstring=False):
    cols = ('id', 'dagscore', 'a_seqid', 'b_seqid', 'dir', 'ngenes')
    #1  17397.0 athaliana_1 athaliana_1 f   432
    if asstring:
        header = header.replace(JS, "!!")
    li = header[1:-1].split('\t')
    if asstring:
        return JS.join(li)
    li[0], li[-1] = int(li[0]), int(li[-1])
    li[1] = float(li[1])
    return dict(zip(cols, li))


def get_merge_gene(fh, header=[None]):
    if header[0] is None:
        header[0] = fh.readline()
    line = fh.readline()
    genes = []
    while line and line[0] != "#": 
        d = DagLine(line)
        genes.append(d)
        line = fh.readline()
    if len(genes) == 0: return None, None
    l = header[0]
    header_string = parse_pyheader(header[0], asstring=True)
    # save the next header.
    header[0] = line

    # header string is joined with JS
    reverse = JS + "r" + JS in header_string

    a_start = min(g.a_start for g in genes)
    a_end   = max(g.a_end for g in genes)

    b_start = min(g.b_start for g in genes)
    b_end   = max(g.b_end for g in genes)
    if reverse: b_start, b_end = b_end, b_start

    d = {'a_seqid': genes[0].a_seqid,
         'b_seqid': genes[0].b_seqid,
         'a_accn': 'a' + header_string,
         'b_accn': 'b' + header_string,
         'a_start': a_start, 
         'b_start': b_start, 
         'a_end': a_end, 
         'b_end': b_end, 
         'evalue': 1e-250}
    return DagLine.from_dict(d), header_string
    

def parse_file(dag_file, evalue_cutoff, ignore_dist, merge_genes=False):
    """ if dag_file is "-", then the stuff is read from stdin. """

    accn_info = {}
    matches = {}
    fh = open(dag_file) if dag_file != "-" else sys.stdin
    dag = True
    while dag:
        if merge_genes:
            dag, dag_header = get_merge_gene(fh)
            if dag is None: break

        else:
            dag, dag_header = get_dag_line(fh)
            if dag is None: break

            if dag.evalue >= evalue_cutoff: continue
            if dag.a_seqid == dag.b_seqid:
                if abs(dag.a_start - dag.b_start) < ignore_dist: continue
                if dag.a_accn == dag.b_accn: continue
        
        if not dag.a_accn in accn_info:
            mid = int((dag.a_start + dag.a_end + 0.5) / 2)
            a_feat = {'accn': dag.a_accn, 'start': dag.a_start, 'end': dag.a_end, 'mid': mid, 'seqid': dag.a_seqid}
            accn_info[dag.a_accn] = a_feat
        else:
            a_feat = accn_info[dag.a_accn]

        if not dag.b_accn in accn_info:
            mid = int((dag.b_start + dag.b_end + 0.5) / 2)
            b_feat = {'accn': dag.b_accn, 'start': dag.b_start, 'end': dag.b_end, 'mid': mid, 'seqid': dag.b_seqid}
            accn_info[dag.b_accn] = b_feat
        else:
            b_feat = accn_info[dag.b_accn]
    
        # always sort by seqid and order. 
        if dag.a_seqid > dag.b_seqid:
            a_feat, b_feat = b_feat, a_feat

        elif dag.a_seqid == dag.b_seqid and a_feat['mid'] > b_feat['mid']:
            a_feat, b_feat = b_feat, a_feat


        seqid_key = (a_feat['seqid'], b_feat['seqid'])
        if not seqid_key in matches: matches[seqid_key] = {}
        these_matches = matches[seqid_key]

        if a_feat['accn'] < b_feat['accn']:
            accn_key = a_feat['accn'], b_feat['accn']
        else:
            accn_key = b_feat['accn'], a_feat['accn']

        if accn_key in these_matches:
            if dag.evalue < these_matches[accn_key]['evalue']: these_matches[accn_key]['evalue'] = dag.evalue
        else:
            these_matches[accn_key] = {'A': a_feat, 'B': b_feat, 'evalue': dag.evalue}
            these_matches[accn_key]['diag_str'] = dag_header

    get_dag_line.header = None
    return matches


def parse_cheader(header):
    """ dagchainer.cpp sends a header line: ">Alignment #%d  score = %.1f\n" 
    we just want the 2 numbers.
    """
    stuff = header.split()
    return int(stuff[1][1:]), float(stuff[-1][:-1])


def run_dag_chainer(a_seqid, b_seqid, filename, matches, reverse, options,
                    child_conn,
                   dagchainer=os.path.join(os.path.abspath(os.path.dirname(__file__)), "dagchainer")):
    """
    calls dagchainer and yields groups of matches
    """
    o = options
    cmd = "%(dagchainer)s -G %(gap_length)s -O %(gap_init)s -E %(gap_extend)s -S " +\
          "%(min_score)s -D %(max_dist)s  -F %(filename)s %(reverse)s" # > %(tmp_file)s";
    cmd = cmd % dict(gap_length=o.gap_dist, gap_init=o.gap_init, 
                     gap_extend=o.gap_extend, min_score=o.min_score, 
                     max_dist=o.gap_max, filename="-", reverse=reverse,
                    dagchainer=dagchainer)
    num2pair = matches.values()
    process = Popen(cmd, stdin=PIPE, stdout=PIPE, bufsize=8*4096, shell=True)
    write = process.stdin.write
    for i, pair in enumerate(num2pair):
        write("%i\t%i\t%i\t%.4f\n" % (i, pair['A']['mid'], pair['B']['mid'], scoringF(pair['evalue'])))
    process.stdin.close()

    header = None
    all_data = [] # added instead of yield to allow parallelization.
    data = []
    for line in process.stdout:
        if line[0] == ">":
            if header is None:
                header = parse_cheader(line[1:].strip())
            else:
                if len(data) >= o.min_aligned_pairs:
                    #yield header, data
                    dag_num, dag_score = header
                    all_data.append((dag_num, dag_score, data))
                header = parse_cheader(line[1:].strip())
                data = []
            continue

        #index, pair_id, pos1, pos2, match_score, dag_chain_score = line.strip().split()
        pair_id, dag_chain_score = line.rstrip("\n").split(" ")
        pair = num2pair[int(pair_id)]
        data.append({'pair': pair, 'dag_score': float(dag_chain_score)})

    if len(data) >= o.min_aligned_pairs:
        dag_num, dag_score = header
        all_data.append((dag_num, dag_score, data))
    child_conn.send(all_data)
    child_conn.close()
    
def print_alignment(dir, diag_num, dag_score, group, opts, out):
    # dir is 'f' or 'r'
    # diag_id dagchainer score, a_seqid, b_seqid, dir, npairs
    header_fmt = "#%i\t%.1f\t%s\t%s\t%s\t%i" 

    d = group[0]['pair']
    print >>out, header_fmt % \
            (diag_num, dag_score, d['A']['seqid'], 
             d['B']['seqid'], dir, len(group))

    for pair_dict in group:
        A = pair_dict['pair']['A']
        B = pair_dict['pair']['B']
        print >>out, "%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%e\t%d" % (\
                     A['seqid'], A['accn'], A['start'], A['end'],
                     B['seqid'], B['accn'], B['start'], B['end'],
                     pair_dict['pair']['evalue'], pair_dict['dag_score'])

def run_and_print(all_matches, opts, out=sys.stdout): 
    filename = "-" # tells dagchainer to read from stdin.
    # if out is False, it means we dont want to print, and so we 
    # dont print.
    print_genes = bool(out)
    merge_ids = []

    for (a_seqid, b_seqid), matches in sorted(all_matches.iteritems()):

        parent_connf, child_connf = mPipe()
        pf = Process(target=run_dag_chainer, args=(a_seqid, b_seqid, filename, matches, "", opts, child_connf))
        pf.start()

        parent_connr, child_connr = mPipe()
        pr = Process(target=run_dag_chainer, args=(a_seqid, b_seqid, filename, matches, "-r", opts, child_connr))
        pr.start()

        for dag_num, dag_score, group in parent_connf.recv():
            if print_genes:
                print_alignment('f', dag_num, dag_score, group, opts, out)
            else:
                # for merged merge we just keep the direction and the 'accn' where
                # the 'accn' is actually just the diag_id for the case of a merge
                # run. this is used later to merge the merge with the genes.
                # since the 'A' and 'B' accn are the same (except for the starting
                # letter, just keep 'A'.
                merge_ids.append(('f', a_seqid, b_seqid,
                                 [g['pair']['A']['accn'][1:] for g in group], 
                                 len(group)))

        for dag_num, dag_score, group in parent_connr.recv():
            if print_genes:
                print_alignment('r', dag_num, dag_score, group, opts, out)
            else:
                merge_ids.append(('r', a_seqid, b_seqid,
                                 [g['pair']['A']['accn'][1:] for g in group],
                                 len(group)))

        pr.join()
        pf.join()

    return merge_ids

######################
## merge diags stuff ##
######################
"""
# all_matches
{('athaliana_2', 'athaliana_3'): {('AT2G47870', 'AT3G62950'): {'A': {'seqid': 'athaliana_2', 'start': 19610409, 'accn': 'AT2G47870', 'end': 19610720, 'mid': 19610564}, 'evalue': 2.3403500000000001e-09, 'B': {'seqid': 'athaliana_3', 'start': 23277224, 'accn': 'AT3G62950', 'end': 23277909, 'mid': 23277566}, 'diag_str': '1^1945.0^athaliana_2^athaliana_3^f^158'}, ('AT2G40820', 'AT3G56480'): {'A': {'se 

# merge
[('f', [('a8^252.0^athaliana_1^athaliana_1^f^18', 'b8^252.0^athaliana_1^athaliana_1^f^18'), ('a4^100.0^athaliana_1^athaliana_1^r^5', 'b4^100.0^athaliana_1^athaliana_1^r^5')]), ('f', [('a7^486.0^athaliana_1^athaliana_1^f^43', 'b7^486.0^athaliana_1^athaliana_1^f^43'), ('a1^6750.0^athaliana_1^athaliana_1^f^414', 'b1^6750.0^athaliana_1^athaliana_1^f^414')]), ('f', [('a11^154.0^athaliana_1^athaliana_1^
"""
def merge_merge(merge, all_matches, opts, out):
    """ merge the merge genes with the 
    original dag data sent in"""
    #cols = ('id', 'dagscore', 'a_seqid', 'b_seqid', 'dir', 'ngenes')
    by_diag = matches_by_diag_id(all_matches)
    # wnat the longest first. then we remove shorter ones that are completely contained
    # in the larger ones.
    seen = {}
    merge.sort(key=operator.itemgetter(4), reverse=True)
    for i, (direction, a_seqid, b_seqid, diag_str_list, llen) in enumerate(merge):
        # so here we have a list of merge-diags merged into a single diag... 
        dags = []

        # and we go through and merge them into dags[].
        # TODO: need to sort this out better. can merge-diags with opposite directions
        # be merged??? when?
        for diag_str in diag_str_list:
            if not diag_str in by_diag: continue
            these_dags = by_diag[diag_str]
            #for d in by_diag[diag_str]:
            #    these_dags.append(d)  
            unseen_dags = [t for t in these_dags if not (t.a_accn, t.b_accn) in seen]
            for ud in unseen_dags: seen[(ud.a_accn, ud.b_accn)] = None

            if len(unseen_dags) >= opts.min_aligned_pairs / 2. \
                       or len(unseen_dags) / float(len(these_dags)) > 0.60:
                dags.extend((str(ud) for ud in unseen_dags))

        if len(dags) >= opts.min_aligned_pairs:
            header = "#" + "\t".join([str(i), "100.0", a_seqid, b_seqid, direction,
                                 str(len(dags))])
            print >>out, header
            print >>out, "\n".join(dags)

def matches_by_diag_id(matches):
    """ take the structure returned by parse_file and return a dictionary
    where the keys are the diag_ids and the values are the dag-pair."""
    by_diag = collections.defaultdict(list)
    for seqid_pair, accn_pair_dict in matches.iteritems():
        for accn_pair_key, accn_pair in accn_pair_dict.iteritems():
            assert accn_pair['diag_str'] is not None
            by_diag[accn_pair['diag_str']].append(DagLine.from_pair_dict(accn_pair))
      
    return dict(by_diag) 

def adjust_opts_for_merge(opts):
    opts.min_aligned_pairs = 1
    opts.min_score = int(opts.min_aligned_pairs * 0.5 * opts.max_match_score)
    opts.gap_dist = opts.gap_dist_merge if opts.gap_dist_merge != 0 else 4 * opts.gap_dist 
    opts.gap_max = opts.gap_max_merge if opts.gap_max_merge != 0 else 5 * opts.gap_max
    assert 0 < opts.gap_dist < opts.gap_max

if __name__ == "__main__":
    import optparse, sys
    p = optparse.OptionParser()

    p.add_option('-i', dest='dag', help="""dag file with format
a_seqid<tab>a_accn<tab>a_start<tab>a_end<tab>b_seqid<tab>b_accn<tab>b_start<tab>b_end<tab>e-value""")
    p.add_option('-o', dest='gap_init', type='float', default=0, 
                help="gap open penalty")
    p.add_option('-e', dest='gap_extend', type='float', default=-3, 
                help="gap extension penalty")

    p.add_option('-x', dest="min_score", type='float', default=None,
                 help="minimum alignment score. alignment stops when score " + \
                 " below this value")
    p.add_option('-g', dest='gap_dist', type='float', default=10000, 
                help="averge distance expected between 2 syntenic genes")
    p.add_option('--gm', dest='gap_dist_merge', type='float', default=0, 
                help="averge distance expected between 2 syntenic merged genes"
                " only applicable if --merge is specified")

    p.add_option('-D', dest='gap_max', type='float', default=200000,
                help="maximum distance between 2 matches")
    p.add_option('--Dm', dest='gap_max_merge', type='float', default=0, 
                help="maximum distance expected between 2 syntenic merged genes"
                " only applicable if --merge is specified")

    p.add_option('-E', dest='evalue', type='float', default=1e-5,
                help="maximum evalue.")
    p.add_option('-A', dest='min_aligned_pairs', type='int', default=6,
                help="minimum number of pairs to be considered a diagonal")

    p.add_option('-I', dest='ignore_dist', default=25, type='int',
                help="ignore hits on teh same chromosome within this distance"
                " removes the self-self diagonal")

    p.add_option('-M', dest='max_match_score', type='float', default=50,
                help="maximum score to be assigned to a match")
    p.add_option('--merge', dest='merge', default=None,
                 help="""\
                 path to a file to send the output.
                 when this is specified, the the output is sent to the
                 specified file. then dagchainer is re-run with --gm and
                 --Dm (corresponding to -g and -D in this help menu. but
                 run with each diagonal in the original output as a single
                 gene the resulting 'merged'-genes are run as normal
                 through dagchainer.
                 the final ouput file with the name of this value + ".merge",
                 will contain genes merged into merge groups.
                 """)

    # dag file can also be sent in as the first arg.
    opts, maybe_dag = p.parse_args() 

    if not (opts.dag or maybe_dag):
        sys.exit(p.print_help())
    if not opts.dag: opts.dag = maybe_dag[0]

    if opts.min_score is None:
        opts.min_score = int(opts.min_aligned_pairs * 0.5 * opts.max_match_score)


    # so here, we run the original dag_chainer. and save to opts.merge.
    all_matches = parse_file(opts.dag, opts.evalue, opts.ignore_dist, merge_genes=False)
    out_file = open(opts.merge, 'wb') if opts.merge else sys.stdout
    run_and_print(all_matches, opts, out=out_file)




    if opts.merge:
        out_file.close()
        # then we read that output file, merging the genes.
        merged_dags = parse_file(opts.merge, opts.evalue, opts.ignore_dist, merge_genes=True)

        # doh! still need to re-parse original. but it's pretty short compared to the original
        # opts.dag so it's pretty fast.
        unmerged_dags = parse_file(opts.merge, opts.evalue, opts.ignore_dist, merge_genes=False)

        # then we run and print without printing.
        adjust_opts_for_merge(opts)
        # run it to get the merge, but dont print.
        merge = run_and_print(merged_dags, opts, out=False)
        out_file = open(opts.merge + ".merge", 'wb') 
        merge_merge(merge, unmerged_dags, opts, out_file)

