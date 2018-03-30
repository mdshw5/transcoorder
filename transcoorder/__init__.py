from simplesam import DefaultOrderedDict
from collections import OrderedDict

__version__ = '0.0.1'


def build_sam_header_from_fasta(fasta):
    """ Make a new header for the genomic SAM file from the chromosomes
    in the assembly file """
    header = DefaultOrderedDict(OrderedDict)
    header['@HD']['VN:1.0'] = ['SO:unknown']
    for rec in fasta:
        header['@SQ']['SN:{name}'.format(name=rec.name)] = [
            'LN:{len}'.format(len=len(rec))
        ]
    return header


def transcript_sam_to_genomic_sam(sam, db, transcript):
    """ Transforms a simplesam.Sam from transcriptome coordinates to
    genomic coordinates.

    Arguments:
        sam: simplesam.Sam object with transcriptome coordinates
        transcript:  gffutils.Feature object representing transcript of interest
    Returns:
        simplesam.Sam object with genomic coordinates

    """
    # set ZT tag for grouping reads by their original transcript
    sam['ZT'] = sam.rname

    exons = db.children(transcript, featuretype='exon', order_by='start')
    genome_coords = [(exon.start, exon.end) for exon in exons]
    genome_offset = genome_coords[0][0] - 1
    transcript_coords = [(start - genome_offset, end - genome_offset)
                         for start, end in genome_coords]

    # determine the exon read starts in
    sam.rname = transcript.seqid
    for coord in transcript_coords:
        if sam.pos <= coord[1] and sam.pos >= coord[0]:
            # if read falls entirely within exon
            if coord[1] - sam.pos > len(sam):
                sam.pos += genome_offset
                return sam
