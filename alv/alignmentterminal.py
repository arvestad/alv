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

    def output_alignment(self, al, painter, args):
        '''
        Output alignment al to stdout in blocks of width at most w with colors from painter.
        '''
        self.left_margin = 1 + al.accession_widths()
        assert self.left_margin < self.width - 10

        if self.sorting == 'alpha':
            accessions = al.sorted_accessions()
        elif self.sorting == 'fixed':
            accessions = self.order
        elif self.sorting == 'by identity':
            accessions = al.sort_by_identity(self.order)
        else:
            accessions = al.accessions()
        accessions = list(accessions)
            
        columns_per_block = al.block_width(self.width, args)
        for block in al.blocks(columns_per_block, args):
            for acc in accessions:
                colored_subseq = al.apply_painter(acc, block, painter)
                print("{0:{width}}{1}".format(acc, colored_subseq, width=self.left_margin))
            print(' ' * self.left_margin, block.start, sep='')


    def print_one_sequence_block(self, record, left_margin, start, block_width):
        colored_string = colorize_sequence_string(rec.seq[start : start + block_width])
        print("{0:{width}}{1}".format(rec.id, colored_string, width=left_margin))


    #
    # Needs reimplementation: Loop over columns until a block is full
    # 

    # for start in range(0, al_width, block_width):
    #     for record in alignment:
    #         painter.print_one_sequence_block(record, left_margin, start, block_width)
    #     print(' ' * left_margin, start, sep='')

        
