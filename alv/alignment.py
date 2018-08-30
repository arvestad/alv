import Bio.Seq
import itertools
from collections import Counter
from .exceptions import AlvEmptyAlignment

class BaseAlignment:
    '''
    Base class for alignments. It contains the interface to use
    to access all alignments, and implementations for the subclasses
    with single-letter column widths.
    '''
    def __init__(self, alignment):
        self.al = alignment   # Holder of the BioPython alignment object
        self.type = None
        self.column_width = 1
        self._update_seq_index()
        self.columns = self._summarize_columns()
        self.basic_info = {'Number of sequences': len(self.al),
                           'Alignment width': self.al.get_alignment_length()}        


    def _update_seq_index(self):
        self.seq_indices = { r.id : i for i, r in enumerate(self.al)} # Get a dictionary mapping accession to row index in alignment
        
    def al_width(self):
        '''
        The number of columns in the alignment.
        '''
        return self.al.get_alignment_length()

    def trim_accessions(self, start, stop):
        for record in self.al:
            acc = record.id
            new_acc = acc[start:stop]
            record.id = new_acc
        self._update_seq_index()

    def abbreviate_accessions(self, n_chars):
        for record in self.al:
            acc = record.id
            new_acc = acc[0:n_chars] + '*' + acc[-n_chars:]
            record.id = new_acc
        self._update_seq_index()
        
    def accessions(self):
        '''
        The accessions in arbitrary order?
        '''
        return map(lambda r: r.id, self.al)

    def sorted_accessions(self):
        '''
        Accessions in alphabetical order.
        '''
        return sorted(self.accessions())

    def sort_by_identity(self, acc):
        '''
        Return accessions in order by similarity (as percent identity) to sequence 'acc'.
        '''
        if acc not in self.seq_indices:
            raise ValueError('Accession "' + acc + '" not found in the alignment')
        pivot_sequence_record = self.al[self.seq_indices[acc]]
        sorted_accessions = map(lambda rec: rec.id,
                                sorted(self.al,
                                       key=lambda rec: percent_identity(pivot_sequence_record.seq, rec),
                                       reverse=True))
        return sorted_accessions

    def accession_widths(self, accessions=None):
        '''
        Compute the space needed for all accessions, plus one for
        a delimiter.
        '''
        max_accession_length = 5        # initial guess, or minimum
        for record in self.al:
            if not accessions or record.id in accessions:
                if len(record.id) > max_accession_length:
                    max_accession_length = len(record.id)
        return max_accession_length


    def block_width(self, terminal_width, args_width):
        '''
        For wide alignments, we need to break it up in blocks.
        This method calculates how many characters to output in a block.
        
        Take the margin size (for accessions) into account and avoid ending up
        with blocks of size 1.
        '''
        if args_width == 0:
            al_width = self.al.get_alignment_length()
            left_margin = 1 + self.accession_widths() # Add 1 for a space to the right of the accessions
            return self._compute_block_width(terminal_width, al_width, left_margin)
        else:
            return args_width

    def _compute_block_width(self, terminal_width, al_width, left_margin):
        '''
        Helper for block_width().
        '''
        last_block_width = al_width % (terminal_width - left_margin)
        n_blocks = al_width // (terminal_width - left_margin) # Not counting last block, if uneven

        if last_block_width == al_width:
            # Short alignment, fits in one block
            return al_width
        elif last_block_width > 10:
            return terminal_width - left_margin
        else:
            sacrifice = max([1, (10 - last_block_width) // n_blocks])
            return terminal_width - left_margin - sacrifice

    def blocks(self, block_width):
        al_width = self.al.get_alignment_length()
        if al_width == 0:
            raise AlvEmptyAlignment()
        else:
            for start in range(0, al_width, block_width):
                end = min(al_width, start + block_width)
                yield AlignmentBlock(start, end)

    def apply_painter(self, acc, block, painter):
        '''
        Colorize a subsequence according to some style (painter).
        This implementation works for aa and plain DNA/RNA.
        '''
        idx = self.seq_indices[acc]
        seq_record = self.al[idx, block.start:block.end]
        colored_seq = ''
        for col_no, c in enumerate(seq_record.seq):
            colored_seq += painter.colorizer(c, self.columns[block.start + col_no])
        return colored_seq + painter.eol()
        
    def _summarize_columns(self):
        '''
        Count the different elements in each column.
        '''
        columns = []
        for col_no in range(self.al.get_alignment_length()):
            columns.append(Counter(self.al[:, col_no]))
        return columns

    def get_basic_info(self):
        '''
        Generator for pairs of description and value, e.g., "alignment width" and an integer.
        '''
        for descr, val in self.basic_info.items():
            yield descr, val

class aaAlignment(BaseAlignment):
    def __init__(self, alignment):
       	super().__init__(alignment)
        self.type = 'aa'

class dnaAlignment(BaseAlignment):
    def __init__(self, alignment):
       	super().__init__(alignment)
        self.type = 'dna'

class codonAlignment(BaseAlignment):
    '''
    Alignment of coding DNA. A column has a width of three nucleotides.
    '''
    def __init__(self, alignment):
        super().__init__(alignment)
        self.type = 'codon'
        self.column_width = 3
        self.genetic_code = 1   # The standard code
        self.basic_info['Genetic code'] = self.genetic_code

    def block_width(self, terminal_width, args):
        '''
        Refinement of the superclass' implemention to ensure that all blocks have 
        a width that is a multiple of three.
        '''
        nominal_width = super().block_width(terminal_width, args)
        remainder = nominal_width % 3
        return nominal_width - remainder # Safe, because the superclass method guarantees width >= 10.

    def apply_painter(self, acc, block, painter):
        '''
        Colorize a subsequence according to some style (painter).
        This implementation works for codons.
        '''
        idx = self.seq_indices[acc]
        seq_record = self.al[idx, block.start:block.end]
        seq = str(seq_record.seq)
        colored_seq = ''

        for codon_col, pos in enumerate(range(0, len(seq), 3)):
            c = seq[pos:pos+3]
            colored_seq += painter.colorizer(c, self.columns[block.start // 3 + codon_col])
        return colored_seq + painter.eol()
        
    def _summarize_columns(self):
        '''
        Specialization of base method for codon columns. Do not focus on the amino acids, but look at
        amino acid columns.
        '''
        columns = []
        for pos in range(0, self.al.get_alignment_length(), 3):
            codon_column = map(lambda r: str(r.seq), self.al[:, pos:pos+3])
            aa_column = map(lambda codon: self._translate(codon), codon_column)
            columns.append(Counter(aa_column))
            #columns.append(Counter(codon_column))
        return columns

    def _translate(self, codon):
        try:
            if codon == '---':
                return '-'
            else:
                aa = Bio.Seq.translate(codon, table = self.genetic_code)
                return aa
        except:
            # For when we have weird codons/alignments
            return 'X'

    def set_genetic_code(self, code):
        self.genetic_code = code
        self.basic_info['Genetic code'] = code


class AlignmentBlock:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        

def percent_identity(seq1, seq2):
    identical = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 == c2:
            identical +=1
    return identical / len(seq1)
