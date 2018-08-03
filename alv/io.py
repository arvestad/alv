from Bio import AlignIO
from .alignment import aaAlignment, dnaAlignment, codonAlignment
from .colorize import aaPainter, dnaPainter, codonPainter

def read_alignment(filename, args):
    '''
    Factory function. Read the alignment with BioPython's support, and 
    return an appropriate alv alignment, by looking at input options.
    '''
    seqtype = args.type
    input_format = args.format

    alignment = AlignIO.read(args.infile, input_format)
    if seqtype == 'aa':
        return aaAlignment(alignment), aaPainter(args)
    if seqtype == 'dna':
        return dnaAlignment(alignment), dnaPainter(args)
    elif seqtype == 'codon':
        return codonAlignment(alignment), codonPainter(args)
    elif seqtype == 'guess':
        raise Exception('Not implemented yet')
    else:
        raise Exception('Unknown option')

    return alignment, painter




            
