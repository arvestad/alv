from .io import get_alv_objects
from .alignmentterminal import AlignmentTerminal

def view(msa, seqtype='guess', color_scheme='guess', genetic_code=1, width=80, dotted=False):
    '''
    Present a colorised version of a multiple sequence alignment (MSA).

    Options:
    * msa           An alignment as BioPython object (from AlignIO).
    * seqtype       String describing the type of sequence. One of 'aa', 'dna',
                    'rna', 'codon', or 'guess' (which almost always works). 
    * color_scheme  String indicating what coloring scheme to use. Ignore unless
                    you want to use 'taylor' or 'hydrophobicity'.
    * genetic_code  An integer indicating which genetic code to use. In the 
                    range 1 to 33 as of this writing, with 1 being the standard
                    code. See NCBI's list at 
                       https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    * width         Integer indicating the desired width of alignment view. 
                    Default: 60 characters.
    * dotted        If True, all characters in the first sequence of the MSA is
                    written and remaining sequences are dotted everwhere except
                    for where they differ from the first sequence.
    '''
    alignment, painter = get_alv_objects(msa, seqtype, color_scheme, genetic_code)
    terminal = AlignmentTerminal(width)
    terminal.output_alignment(alignment, painter, width - 12, dotted)


def glimpse(msa, n_seq=20, seqtype='guess', color_scheme='guess', genetic_code=1, width=60, dotted=False):
    '''
    Works like "view", but extracts a conserved region and outputs a random sample of the sequences or all of them if it is a small alignment.

    Options:
    * msa           An alignment as BioPython object (from AlignIO).
    * n_seq         The number of sequences to output.
    * seqtype       String describing the type of sequence. One of 'aa', 'dna',
                    'rna', 'codon', or 'guess' (which almost always works). 
    * color_scheme  String indicating what coloring scheme to use. Ignore unless
                    you want to use 'taylor' or 'hydrophobicity'.
    * genetic_code  An integer indicating which genetic code to use. In the 
                    range 1 to 33 as of this writing, with 1 being the standard
                    code. See NCBI's list at 
                       https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    * width         Integer indicating the desired width of alignment view. 
                    Default: 60 characters.
    * dotted        If True, all characters in the first sequence of the MSA is
                    written and remaining sequences are dotted everwhere except
                    for where they differ from the first sequence.
    '''
    alignment, painter = get_alv_objects(msa, seqtype, color_scheme, genetic_code)
    terminal = AlignmentTerminal(width)
    terminal.height = n_seq + 2
    terminal.output_glimpse(alignment, painter, width, dotted)

