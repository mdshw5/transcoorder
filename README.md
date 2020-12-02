# Transcoorder

Convert SAM/BAM files with transcript coordinates to SAM files with genomic coordinates

```
$ transcoord -h
usage: transcoord [-h] [-o OUT] [-t TAG_NAME] [--debug] [--version] gtf bam fasta

positional arguments:
  gtf                   GTF file containing transcripts
  bam                   SAM or BAM files aligned to transcriptome
  fasta                 FASTA format assembly coresponding to GTF

optional arguments:
  -h, --help            show this help message and exit
  -o OUT, --out OUT     output file for genomic SAM (default: stdout)
  -t TAG_NAME, --tag-name TAG_NAME
                        SAM tag name for storing transcript identifier. default: ZT
  --debug               enable debugging
  --version             display version number
```

Note: This is an incomplete work in progress script and you shouldn't rely on its output. I 
made this as a proof of concept, and never had time to finish. It currently runs very slowly 
and only handles reads that are mapped entirely within an exon.