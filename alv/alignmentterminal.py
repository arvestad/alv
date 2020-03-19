from .alignment import AlignmentBlock
import random
import sys
import shutil

class AlignmentTerminal:
    '''
    This class encapsulates knowledge about the terminal and how to draw on it.
    '''
    def __init__(self, args):
        self.width, self.height = shutil.get_terminal_size() # First item in this tuple is the width, the other is the height.
        self.random_sample_size = args.random_accessions
        self.sorting = args.sorting
        if args.sort_by_id:
            self.sorting = 'by identity'
            self.order = args.sort_by_id
        elif args.sorting_order:
            self.sorting = 'fixed'
            self.order = args.sorting_order.split(',')
            if len(self.order) == 0:
                raise Exception('Bad order specification: no accessions in input')
        if args.select_matching:
            self.selection = args.select_matching
        else:
            self.selection = False

    def get_accession_list(self, al):
        '''
        Choose accessions and their order based on the user-choices.
        '''
        if self.sorting == 'alpha':
            accessions = al.sorted_accessions()
        elif self.sorting == 'fixed':
            accessions = self.order
        elif self.sorting == 'by identity':
            accessions = al.sort_by_identity(self.order)
        else:
            accessions = al.accessions()

        if self.selection:
            chosen_accessions = []
            for acc in accessions:
                if self.selection in acc:
                    chosen_accessions.append(acc)
            return chosen_accessions
        else:
            return list(accessions)


    def _output_block(self, al, painter, width, chosen_accessions, block, dotted):
        '''
        Write one block of an aligmnent to the terminal.
        The accessions to the left, sequences to the right, and at the bottom a tick string.
        '''
        if dotted:
            template_acc=chosen_accessions[0]
            chosen_accessions = chosen_accessions[1:]
            colored_subseq = al.apply_painter(template_acc, block, painter)
            print("{0:{width}}{1}".format(template_acc, colored_subseq, width=self.left_margin))
            
        for acc in chosen_accessions:
            if dotted:
                colored_subseq = al.apply_dotter(acc, block, painter, template_acc)
            else:
                colored_subseq = al.apply_painter(acc, block, painter)
            print("{0:{width}}{1}".format(acc, colored_subseq, width=self.left_margin))
        print(make_tick_string(self.left_margin, block.start, block.end, 20, 7))


    def _setup_left_margin(self, al, chosen_accessions):
        '''
        Compute the width of the left margin from the accession lengths
        '''
        self.left_margin = 1 + al.accession_widths(chosen_accessions)
        assert self.left_margin < self.width - 10


    def output_alignment(self, al, painter, width, dotted=False):
        '''
        Output alignment al to stdout in blocks of width at most w with colors from painter.

        Args:
          al -- the alignment object
          painter -- Painter object that decides colors for symbols
          width -- How many columns to use for one alignment block
          dotted -- Flag. Output periods (.) if position identical to first sequence?
        '''
        chosen_accessions =  self.get_accession_list(al)

        if self.random_sample_size and len(chosen_accessions) > self.random_sample_size:
            chosen_accessions = random.sample(chosen_accessions, self.random_sample_size)

        self._setup_left_margin(al, chosen_accessions)

        columns_per_block = al.block_width(self.width, width)
        for block in al.blocks(columns_per_block):
            self._output_block(al, painter, width, chosen_accessions, block, dotted)


    def output_glimpse(self, al, painter, width, dotted=False):
        '''
        Output a single-screen glimpse of the alignment. An attempt at finding
        the most interesting (guessed to be the most conserved part of a random
        sample of sequences) of the alignment is done.
        '''
        chosen_accessions =  self.get_accession_list(al)

        n_seqs_to_view = min(len(chosen_accessions), self.height - 2) # -2 to make room for tick line and next prompt
        chosen_accessions = random.sample(chosen_accessions, n_seqs_to_view)

        self._setup_left_margin(al, chosen_accessions)

        n_columns = al.block_width(self.width, width) # This many alignment columns
        conserved_block = al.get_conserved_block(n_columns)

        self._output_block(al, painter, width, chosen_accessions, conserved_block, dotted)


def calc_tick_indices(start, end, distance, min_distance):
    '''
    Return a list of indices for which we want a tick mark at the bottom of the alignment.
    The goal is to have an index for the starting position of a block (leftmost column number),
    and then a tick mark on even multiples of 20 (or what is given by 'distance'), for example:
        53   60                  80                 100
    Care is needed so that space is left between first and second indices, and min_distance indicates
    how much.
    '''
    first_even_pos = (start // distance + 1) * distance
    if first_even_pos - start < min_distance:
        first_even_pos += distance # Compensate a bit
    positions = range(first_even_pos, end, distance)
    return positions

def make_one_tick(position, space, tickmark):
    '''
    Return a string which is 'space' wide and contains a number (the position)
    followed by an up-arrow.
    '''
    return '{0:>{width}}{tickmark}'.format(position, width=space-1,tickmark=tickmark)

def make_tick_string(left_margin, start, end, distance, min_distance):
    '''
    Construct the index bar which is printed at the bottom of an alignment block.

    left_margin is how much space is allowed for accessions.
    start is the column number of the beginning of an alignment block.
    end is the last column of an alignment block.
    distance is the desired distance between up-arrows
    min_distance is the space we allow for position numbers plus an up-arrow
    '''
    even_indices = calc_tick_indices(start, end, distance, min_distance)

    tickmark = '^'
    if sys.stdout.encoding == 'UTF-8':
        tickmark = '↑'

    # Initial space
    index_bar = ' ' * (left_margin - min_distance + 1) # Account for space needed by indices

    # Add first column index
    index_bar += make_one_tick(start, min(left_margin+1, min_distance), tickmark)

    last_pos = start
    for pos in even_indices:
        spacer = pos - last_pos - 1
        index_block = '{0:>{width}}{tickmark}'.format(pos, width=spacer, tickmark=tickmark)
        index_bar += index_block
        last_pos = pos
    return index_bar
