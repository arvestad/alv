import Bio.Seq
from colorama import init as colorama_init, Fore, Back, Style
import sys

class Painter:
    '''
    Base class for all painters
    '''
    def __init__(self):
        self.restrictions = []
        self.bg_neutral = Back.WHITE
        self.fg_neutral = Fore.BLACK

    def set_options(self, args):
        if args.keep_colors_when_redirecting:
            colorama_init(strip=False)
        else:
            colorama_init()
        if args.majority:
            self.restrictions.append(restrict_to_majority)
        if args.no_indels:
            self.restrictions.append(restrict_to_no_indels)

    def color_mode(self, m):
        '''
        Dark or light background?
        '''
        if m == 'light': # White background
            self.bg_neutral = Back.WHITE
            self.fg_neutral = Fore.BLACK
        elif m == 'dark':   #
            self.bg_neutral = Back.BLACK
            self.fg_neutral = Fore.WHITE

    def color_for_bad_data(self):
        return self.bg_neutral + Fore.RED, self.fg_neutral + self.bg_neutral

    def color_for_stop(self):
        return Back.BLACK + Fore.RED, self.fg_neutral

    def indel_color(self):
        return self.bg_neutral, None

    def eol(self):
        return Style.RESET_ALL

    def sol(self):
        '''
        Execute at start-of-line.
        '''
        return self.fg_neutral


# Restriction functions
#
# These look at columns, given as a Counter object, and decides whether the column
# should be colored or not. Returns True for 'color it', or False for 'dont color it'.
#
def restrict_to_majority(column):
    maj = column.most_common(1)[0] # Single most common
    if maj[0] != '-' and maj[1] >= 0.5 * len(list(column.elements())):
        return True
    else:
        return False

def restrict_to_no_indels(column):
    return column['-'] == 0



class AminoAcidPainter(Painter):
    '''
    Put paint of amino acids.
    '''
    def __init__(self):
        super().__init__()

    def _color_lookup(self, c):
        if c in 'AILMFWVCailmfwvc':
            return Back.BLUE, self.bg_neutral
        elif c in 'KRkr':
            return Back.RED, self.bg_neutral
        elif c in 'EDed':
            return Back.MAGENTA, self.bg_neutral
        elif c in 'NQSTnqst':
            return Back.GREEN, self.bg_neutral
        elif c in 'Gg':
            return Back.YELLOW, self.bg_neutral
        elif c in 'Pp':
            return Back.YELLOW, self.bg_neutral
        elif c in 'HYhy':
            return Back.CYAN, self.bg_neutral
        elif c in 'Xx':
            return self.bg_neutral, self.bg_neutral
        elif c in "!?":
            return self.color_for_bad_data()
        elif c in "*":
            return self.color_for_stop()
        elif c in '-_.:':
            return self.indel_color()
        else:
            return self.bg_neutral, None

    def colorizer(self, c, column):
        if c == '!?':
            before_color, after_color = self.color_for_bad_data()
            return before_color + c + after_color
        elif all(map(lambda r: r(column), self.restrictions)): # True also if the list self.restrictions is empty
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


class AminoAcidTaylorPainter(AminoAcidPainter):
    '''
    Put paint to amino acids, an approximation of the "Taylor" style.
    '''
    def __init__(self):
        super().__init__()

    def _color_lookup(self, c):
        if c in 'AILMFVCailMFVC':
            return Back.GREEN, self.bg_neutral
        elif c in 'KRHkrh':
            return Back.BLUE, self.bg_neutral
        elif c in 'EDSTedst':
            return Back.RED, self.bg_neutral
        elif c in 'NQnq':
            return Back.MAGENTA, self.bg_neutral
        elif c in 'CGPcgp':
            return Back.YELLOW, self.bg_neutral
        elif c in 'Yy':
            return Back.CYAN, self.bg_neutral
        elif c in 'X':
            return self.bg_neutral, self.bg_neutral
        elif c in "!?":
            return self.color_for_bad_data()
        elif c in "*":
            return self.color_for_stop()
        elif c in '-_.:':
            return self.indel_color()
        else:
            return self.bg_neutral, None

class AminoAcidHydrophobicity(AminoAcidPainter):
    '''
    Put paint to amino acids, indicating hydrophobicity.
    '''
    def __init__(self):
        super().__init__()

    def _color_lookup(self, c):
        if c in 'AILMFVPGailmfvpg':
            return Back.RED, self.bg_neutral
        elif c in 'QNHSTYCWqnhstycw':
            return Back.BLUE, self.bg_neutral
        elif c in 'RKDErkde':
            return Back.GREEN, self.bg_neutral
        elif c in 'X':
            return self.bg_neutral, self.bg_neutral
        elif c in "!?":
            return self.color_for_bad_data()
        elif c in "*":
            return self.color_for_stop()
        elif c in '-_.:':
            return self.indel_color()
        else:
            return self.bg_neutral, None


class DnaPainter(Painter):
    '''
    Put paint of nucleotides.
    '''
    def __init__(self):
        super().__init__()

    def colorizer(self, c, column):
        if c in 'TUtu':           # Handles RNA too
            return Back.CYAN + c
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
            return self.bg_neutral + c

class DnaClassPainter(Painter):
    '''
    Put paint of nucleotides, even RNA sequences. (Is this one even used??)
    '''
    def __init__(self):
        super().__init__()

    def colorizer(self, c, column=[]):
        if c in 'TUtuCcYy':           # Handles RNA too
            return Back.Cyan + c
        elif c in 'AaGgRr':
            return Back.Magenta + c
        elif c in '!*':
            before, after = self.color_for_bad_data()
            return before + c + after
        elif c in '-.:':
            before, after = self.indel_color()
            print(before_color, after_color)
            if after:
                return before + c + after
            else:
                return before + c
        else:
            return self.bg_neutral + c


class CodonPainter(Painter):
    '''
    Put paint of codons.
    '''
    def __init__(self, aa_painter):
        super().__init__()
        self.aa_painter = aa_painter

    def colorizer(self, c, column):
        '''
        c is a expected to be a codon, or a single letter.
        Single letters may occur in MACSE alignments, due to frameshifts.
        '''
        try:
            if '!' in c:
                before_color, after_color = self.color_for_bad_data()
                assert after_color # Weird if we don't have a color reset.
                return before_color + c + after_color
            elif all(map(lambda r: r(column), self.restrictions)): # True also if the list self.restrictions is empty
                if len(c) != 3:
                    return self_bg_neutral + Fore.RED + Style.BRIGHT + c + self.fg_neutral + Style.NORMAL
                elif c == '---':
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
            print('Alv warning:', str(e), file=sys.stderr)
            return c                # Temporarily. Adding colors later
