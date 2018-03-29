#!/usr/bin/env python
import argparse
import sys
from simplesam import Reader, Writer
from gffutils import FeatureDB, create_db
from transcoorder import transcript_sam_to_genomic_sam


def main(ext_args=None):
    from transcoorder import __version__
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        'gtf', type=str, help='GTF file containing transcripts')
    parser.add_argument(
        'bam',
        type=Reader,
        nargs='*',
        help="SAM or BAM files aligned to transcriptome")
    parser.add_argument(
        '-o',
        '--out',
        type=Writer,
        default='-',
        help="output file for genomic SAM (default: stdout)")
    parser.add_argument(
        '-t',
        '--tag-name',
        type=str,
        default='ZT',
        help=
        "SAM tag name for storing transcript identifier. default: %(default)s")
    parser.add_argument(
        '--version',
        action="version",
        version=__version__,
        help="display version number")
    # print help usage if no arguments are supplied
    if len(sys.argv) == 1 and not ext_args:
        parser.print_help()
        sys.exit(1)
    elif ext_args:
        args = parser.parse_args(ext_args)
    else:
        args = parser.parse_args()

    try:
        db = FeatureDB(args.gtf + '.db')
    except ValueError:
        sys.stderr.write("building sqlite database for %s..." % args.gtf)
        db = create_db(
            args.gtf,
            args.gtf + '.db',
            disable_infer_transcripts=True,
            disable_infer_genes=True)

    args.out.write()
    for read in args.bam:
        try:
            transcript = db[read.rname]
            transcript_sam_to_genomic_sam(bam, db, transcript)
        except KeyError:
            pass


if __name__ == "__main__":
    main()
