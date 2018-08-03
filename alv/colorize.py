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

    def eol(self):
        return Style.RESET_ALL




class aaPainter(Painter):
    '''
    Put paint of amino acids.
    '''
    def __init__(self, args):
        super().__init__(args)
        
    def _color_lookup(self, c):
        if c in 'AILMFWVC':
            return Back.BLUE, None
        elif c in 'KR':
            return Back.RED, None
        elif c in 'ED':
            return Back.MAGENTA, None
        elif c in 'NQST':
            return Back.GREEN, None
        elif c in 'G':
            return Back.YELLOW, None
        elif c in 'P':
            return Back.YELLOW, None
        elif c in 'HY':
            return Back.CYAN, None
        elif c in "*!":
            return Back.BLACK + Fore.RED, Back.WHITE + Fore.BLACK
        else:
            return Back.WHITE, None

    def colorizer(self, c):
        before_color, after_color = self._color_lookup(c)
        colored_item = before_color + c
        if after_color:
            return c + after_color
        else:
            return colored_item
        

class dnaPainter(Painter):
    '''
    Put paint of nucleotides.
    '''
    def __init__(self, args):
        super().__init__(args)
        
    def colorizer(self, c):
        if c in 'TUtu':           # Handles RNA too
            return Back.BLUE + c
        elif c in 'Aa':
            return Back.GREEN + c
        elif c in 'Cc':
            return Back.YELLOW + c
        elif c in 'Gg':
            return Back.RED + c
        elif c in '!*':
            return Back.BLACK + Fore.RED + Style.BRIGHT + c + Fore.BLACK + Style.NORMAL
        else:
            return Back.WHITE + c

class codonPainter(Painter):
    '''
    Put paint of codons.
    '''
    def __init__(self, args):
        super().__init__(args)
        self.aa_painter = aaPainter(args)

    def colorizer(self, c):
        '''
        c is a expected to be a codon, or a single letter.
        Single letters may occur in MACSE alignments, due to frameshifts.
        '''
        try:
            if len(c) != 3:
                return Back.WHITE + Fore.RED + Style.BRIGHT + c + Fore.BLACK + Style.NORMAL
            else:
                if c == '---':
                    aa = '-'
                else:
                    aa = Bio.Seq.translate(c)
                before_color, after_color = self.aa_painter._color_lookup(aa)
                colored_item = before_color + c
                if after_color:
                    return colored_item + after_color
                else:
                    return colored_item
                
        except Bio.Data.CodonTable.TranslationError:
            return c
        except Exception as e:
            print('Warning:', str(e), file=sys.stderr)
            return c                # Temporarily. Adding colors later
