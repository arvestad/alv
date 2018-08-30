from alv.get_terminal_size import get_terminal_size

class AlignmentTerminal:
    '''
    This class encapsulates knowledge about the terminal and how to draw on it.
    '''
    def __init__(self, args):
        self.width = get_terminal_size()[0] # First item in this tuple is the width, the other is the height.
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

    def output_alignment(self, al, painter, width):
        '''
        Output alignment al to stdout in blocks of width at most w with colors from painter.
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
        else:
            chosen_accessions = list(accessions)

        self.left_margin = 1 + al.accession_widths(chosen_accessions)
        assert self.left_margin < self.width - 10
            
        columns_per_block = al.block_width(self.width, width)
        for block in al.blocks(columns_per_block):
            for acc in chosen_accessions:
                colored_subseq = al.apply_painter(acc, block, painter)
                print("{0:{width}}{1}".format(acc, colored_subseq, width=self.left_margin))
            print(make_tick_string(self.left_margin, block.start, block.end, 20, 7))

#            print(' ' * self.left_margin, '↑',  block.start, sep='') # print index of first column


    # def print_one_sequence_block(self, record, left_margin, start, block_width):
    #     colored_string = colorize_sequence_string(rec.seq[start : start + block_width])
    #     print("{0:{width}}{1}".format(rec.id, colored_string, width=left_margin))


        
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

def make_one_tick(position, space):
    '''
    Return a string which is 'space' wide and contains a number (the position)
    followed by an up-arrow.
    '''

    return '{0:>{width}}↑'.format(position, width=space-1) 

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

    # Initial space
    index_bar = ' ' * (left_margin - min_distance + 1) # Account for space needed by indices

    # Add first column index
    index_bar += make_one_tick(start, min(left_margin+1, min_distance))

    last_pos = start
    for pos in even_indices:
        spacer = pos - last_pos - 1
        index_block = '{0:>{width}}↑'.format(pos, width=spacer)
        index_bar += index_block
        last_pos = pos
    return index_bar
