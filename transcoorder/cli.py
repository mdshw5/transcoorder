#!/usr/bin/env python
import argparse
import sys
from simplesam import Reader, Writer, bam_read_count
from gffutils import FeatureDB, create_db
from pyfaidx import Fasta
from tqdm import tqdm
from transcoorder import transcript_sam_to_genomic_sam, build_sam_header_from_fasta, cache_gtf_features


def main(ext_args=None):
    from transcoorder import __version__
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        'gtf', type=str, help='GTF file containing transcripts')
    parser.add_argument(
        'bam',
        type=argparse.FileType('r'),
        help="SAM or BAM files aligned to transcriptome")
    parser.add_argument(
        'fasta', type=Fasta, help="FASTA format assembly coresponding to GTF")
    parser.add_argument(
        '-o',
        '--out',
        type=argparse.FileType('w'),
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
        '--debug', action="store_true", help="enable debugging")
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

    header = build_sam_header_from_fasta(args.fasta)
    with Reader(args.bam) as bamfile, Writer(args.out, header) as outfile:
        try:
            read_count = len(bamfile)
        except NotImplementedError:
            read_count = None
        with tqdm(total=read_count, unit='read') as pbar:
            for read in bamfile:
                features = cache_gtf_features(db, read.rname)
                if features == None and args.debug:
                    sys.stderr.write("%s not found in %s\n" % (read.rname,
                                                               args.gtf))
                else:
                    transcript, genome_offset, transcript_coords = features
                    read = transcript_sam_to_genomic_sam(
                        read, transcript, genome_offset, transcript_coords)
                    if read:
                        outfile.write(read)
                pbar.update(1)


if __name__ == "__main__":
    main()
