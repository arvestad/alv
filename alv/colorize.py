import Bio.Seq
from colorama import init as colorama_init, Fore, Back, Style
import sys

class Painter:
    '''
    Base class for all painters
    '''
    def __init__(self, args):
        if args.keep_colors_when_redirecting:
            colorama_init(strip=False)
        else:
            colorama_init()
        self.restrictions = []
        if args.majority:
            self.restrictions.append(restrict_to_majority)
        if args.no_indels:
            self.restrictions.append(restrict_to_no_indels)

    def color_for_bad_data(self):
        return Back.WHITE + Fore.RED, Fore.BLACK + Back.WHITE

    def color_for_stop(self):
        return Back.BLACK + Fore.RED, Fore.BLACK

    def indel_color(self):
        return Back.WHITE, None

    def eol(self):
        return Style.RESET_ALL

# Restriction functions
#
# These look at columns, given as a Counter object, and decides whether the column
# should be colored or not. Returns True for 'color it', or False for 'dont color it'.
#
def restrict_to_majority(column):
    maj = column.most_common(1) # Single most common
    if maj >= 0.5 * len(column.elements()):
        return True
    else:
        return False

def restrict_to_no_indels(column):
    return column['-'] == 0

        

class aaPainter(Painter):
    '''
    Put paint of amino acids.
    '''
    def __init__(self, args):
        super().__init__(args)
        
    def _color_lookup(self, c):
        if c in 'AILMFWVC':
            return Back.BLUE, Back.WHITE
        elif c in 'KR':
            return Back.RED, Back.WHITE
        elif c in 'ED':
            return Back.MAGENTA, Back.WHITE
        elif c in 'NQST':
            return Back.GREEN, Back.WHITE
        elif c in 'G':
            return Back.YELLOW, Back.WHITE
        elif c in 'P':
            return Back.YELLOW, Back.WHITE
        elif c in 'HY':
            return Back.CYAN, Back.WHITE
        elif c in 'X':
            return Back.WHITE, Back.WHITE
        elif c in "!?":
            return self.color_for_bad_data()
        elif c in "*":
            return self.color_for_stop()
        elif c in '-_.:':
            return self.indel_color()
        else:
            return Back.WHITE, None

    def colorizer(self, c, column):
        if all(map(lambda r: r(column), self.restrictions)): # True also if the list self.restrictions is empty
            before_color, after_color = self._color_lookup(c)
            colored_item = before_color + c
            if after_color:
                return colored_item + after_color
            else:
                return colored_item
        else:
            before_color, after_color = self.indel_color()
            if after_color:
                return before_color + c + after_color
            else:
                return before_color + c
        

class dnaPainter(Painter):
    '''
    Put paint of nucleotides.
    '''
    def __init__(self, args):
        super().__init__(args)
        
    def colorizer(self, c, column):
        if c in 'TUtu':           # Handles RNA too
            return Back.BLUE + c
        elif c in 'Aa':
            return Back.GREEN + c
        elif c in 'Cc':
            return Back.YELLOW + c
        elif c in 'Gg':
            return Back.RED + c
        elif c in '!*':
            before, after = self.color_for_bad_data()
            return before + c + after
        elif c in '-.:':
            before, after = self.indel_color()
            if after:
                return before + c + after
            else:
                return before + c
        else:
            return Back.WHITE + c

class codonPainter(Painter):
    '''
    Put paint of codons.
    '''
    def __init__(self, args):
        super().__init__(args)
        self.aa_painter = aaPainter(args)

    def colorizer(self, c, column):
        '''
        c is a expected to be a codon, or a single letter.
        Single letters may occur in MACSE alignments, due to frameshifts.
        '''
        try:
            if all(map(lambda r: r(column), self.restrictions)): # True also if the list self.restrictions is empty
                if len(c) != 3:
                    return Back.WHITE + Fore.RED + Style.BRIGHT + c + Fore.BLACK + Style.NORMAL
                elif '!' in c:
                    before_color, after_color = self.color_for_bad_data()
                    assert after_color # Weird if we don't have a color reset.
                    return before_color + c + after_color
                else:
                    if c == '---':
                        aa = '-'
                        before_color, after_color = self.indel_color()
                    else:
                        aa = Bio.Seq.translate(c)
                        before_color, after_color = self.aa_painter._color_lookup(aa)
                        colored_item = before_color + c
                        if after_color:
                            return colored_item + after_color
                        else:
                            return colored_item
            else:
                before_color, after_color = self.indel_color()
                if after_color:
                    return before_color + c + after_color
                else:
                    return before_color + c
                
                
        except Bio.Data.CodonTable.TranslationError:
            return c
        except Exception as e:
            print('Warning:', str(e), file=sys.stderr)
            return c                # Temporarily. Adding colors later


