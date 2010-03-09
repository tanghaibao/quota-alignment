#!/usr/bin/env python
# -*- coding: UTF-8 -*-


"""
%prog gff_file

convert the common gff format into the bed format that we use
"""

import os.path as op
import sys

try:
    import gt
except:
    print >>sys.stderr, "gt module not found, try download and install genometools`" \
            "<http://genometools.org/>"
    sys.exit(1)


def gff_to_bed(gff_file, bed_fh=sys.stdout):

    fi = gt.FeatureIndexMemory()
    fi.add_gff3file(gff_file)

    for seqid in fi.get_seqids():
        for feat in fi.get_features_for_seqid(seqid):
            has_cds = False
            if "ID" in feat.attribs and feat.type!="Chr": 
                id = feat.attribs["ID"]

            for subf in feat:
                if subf.type == 'CDS': 
                    has_cds = True
                    break

            # only print the protein-coding gene models
            if has_cds:
                print >>bed_fh, "\t".join(str(x) for x in \
                    (feat.seqid, feat.start, feat.end, id))


if __name__ == "__main__":

    import optparse

    parser = optparse.OptionParser(__doc__)
    (options, args) = parser.parse_args()

    if len(args) != 1:
        sys.exit(parser.print_help())

    gff_file = args[0]

    gff_to_bed(gff_file)

