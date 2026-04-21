import Bio.Seq
from colorama import init as colorama_init, Fore, Back, Style
import sys


# --- 24-bit RGB helpers ---

def _rgb_bg(r, g, b):
    return f'\033[48;2;{r};{g};{b}m'

def _rgb_fg(r, g, b):
    return f'\033[38;2;{r};{g};{b}m'


# --- Palette definitions ---
# Maps property/letter keys to (before_color, after_color) tuples.

# Named-terminal fallback palette for AA (used by the 'terminal' scheme).
_AA_TERMINAL_PALETTE = {
    'hydrophobic': (Back.BLUE,    Back.WHITE),
    'positive':    (Back.RED,     Back.WHITE),
    'negative':    (Back.MAGENTA, Back.WHITE),
    'polar':       (Back.GREEN,   Back.WHITE),
    'glycine':     (Back.YELLOW,  Back.WHITE),
    'proline':     (Back.YELLOW,  Back.WHITE),
    'aromatic':    (Back.CYAN,    Back.WHITE),
    'unknown':     (Back.WHITE,   Back.WHITE),
}

# ClustalX canonical amino acid colors (from Jalview reference implementation).
_AA_CLUSTAL_PALETTE = {
    'hydrophobic': (_rgb_bg(128, 160, 240),                          Back.WHITE),             # #80A0F0 blue
    'positive':    (_rgb_bg(240,  21,   5) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# #F01505 red
    'negative':    (_rgb_bg(192,  72, 192) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# #C048C0 magenta
    'polar':       (_rgb_bg( 21, 192,  21),                          Back.WHITE),             # #15C015 green
    'glycine':     (_rgb_bg(240, 144,  72),                          Back.WHITE),             # #F09048 orange
    'proline':     (_rgb_bg(192, 192,   0),                          Back.WHITE),             # #C0C000 yellow
    'aromatic':    (_rgb_bg( 21, 164, 164),                          Back.WHITE),             # #15A4A4 cyan
    'unknown':     (Back.WHITE,                                      Back.WHITE),
}

# Taylor amino acid colors (from Jalview reference implementation).
_AA_TAYLOR_PALETTE = {
    'hydrophobic':   (_rgb_bg( 21, 192,  21),                          Back.WHITE),             # #15C015 green
    'basic':         (_rgb_bg(128, 160, 240),                          Back.WHITE),             # #80A0F0 blue
    'polar_charged': (_rgb_bg(240,  21,   5) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# #F01505 red
    'polar':         (_rgb_bg(192,  72, 192) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# #C048C0 magenta
    'special':       (_rgb_bg(192, 192,   0),                          Back.WHITE),             # #C0C000 yellow
    'aromatic':      (_rgb_bg( 21, 164, 164),                          Back.WHITE),             # #15A4A4 cyan
    'unknown':       (Back.WHITE,                                      Back.WHITE),
}

# CINEMA amino acid colors (Parry-Smith et al. 1997, Jalview reference).
_AA_CINEMA_PALETTE = {
    'positive':   (_rgb_bg(  0,   0, 255) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# blue
    'negative':   (_rgb_bg(255,   0,   0) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# red
    'polar':      (_rgb_bg(  0, 255,   0),                          Back.WHITE),             # green
    'hydrophobic':(_rgb_bg(220, 220, 220),                           Back.WHITE),             # light grey
    'glycine':    (_rgb_bg(255, 128,   0),                          Back.WHITE),             # orange
    'proline':    (_rgb_bg(  0, 255, 255),                          Back.WHITE),             # cyan
    'cysteine':   (_rgb_bg(255, 255,   0),                          Back.WHITE),             # yellow
    'unknown':    (Back.WHITE,                                       Back.WHITE),
}

# Zappo amino acid colors (Jalview reference implementation).
_AA_ZAPPO_PALETTE = {
    'aliphatic':  (_rgb_bg(255, 175, 175),                          Back.WHITE),             # #FFAFAF light pink
    'aromatic':   (_rgb_bg(255, 200,   0),                          Back.WHITE),             # #FFC800 orange
    'positive':   (_rgb_bg(100, 100, 255) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# #6464FF blue
    'negative':   (_rgb_bg(255,   0,   0) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# #FF0000 red
    'hydrophilic':(_rgb_bg(  0, 255,   0),                          Back.WHITE),             # #00FF00 green
    'special':    (_rgb_bg(255,   0, 255) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# #FF00FF magenta
    'cysteine':   (_rgb_bg(255, 255,   0),                          Back.WHITE),             # #FFFF00 yellow
    'unknown':    (Back.WHITE,                                       Back.WHITE),
}

# Okabe-Ito colorblind-safe palette for AA — uses clustal groupings.
# Palette: blue #0072B2, vermillion #D55E00, reddish purple #CC79A7, bluish green #009E73,
#          orange #E69F00, yellow #F0E442, sky blue #56B4E9.
_AA_OKABEITO_PALETTE = {
    'hydrophobic': (_rgb_bg(  0, 114, 178) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# #0072B2 blue
    'positive':    (_rgb_bg(213,  94,   0) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# #D55E00 vermillion
    'negative':    (_rgb_bg(204, 121, 167),                          Back.WHITE),             # #CC79A7 reddish purple
    'polar':       (_rgb_bg(  0, 158, 115) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# #009E73 bluish green
    'glycine':     (_rgb_bg(230, 159,   0),                          Back.WHITE),             # #E69F00 orange
    'proline':     (_rgb_bg(240, 228,  66),                          Back.WHITE),             # #F0E442 yellow
    'aromatic':    (_rgb_bg( 86, 180, 233),                          Back.WHITE),             # #56B4E9 sky blue
    'unknown':     (Back.WHITE,                                      Back.WHITE),
}

_AA_HYDROPHOBICITY_PALETTE = {
    'hydrophobic': (Back.RED,   Back.WHITE),
    'polar':       (Back.BLUE,  Back.WHITE),
    'charged':     (Back.GREEN, Back.WHITE),
    'unknown':     (Back.WHITE, Back.WHITE),
}

_DNA_TERMINAL_PALETTE = {
    'A': (Back.GREEN,  Back.WHITE),
    'C': (Back.YELLOW, Back.WHITE),
    'G': (Back.RED,    Back.WHITE),
    'T': (Back.CYAN,   Back.WHITE),
    'U': (Back.CYAN,   Back.WHITE),
    'N': (Back.WHITE,  Back.WHITE),
}

# ClustalX canonical nucleotide colors (from coldna.xml defaults).
_DNA_CLUSTAL_PALETTE = {
    'A': (_rgb_bg(230, 179,  76),                          Back.WHITE),             # orange
    'C': (_rgb_bg(230,  26,  26),                          Back.WHITE),             # red
    'G': (_rgb_bg( 26, 230,  26),                          Back.WHITE),             # green
    'T': (_rgb_bg( 26,  26, 230) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# blue → white text
    'U': (_rgb_bg( 26,  26, 230) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),
    'N': (Back.WHITE,                                      Back.WHITE),
}


# --- Grouping definitions ---
# Maps individual letters to palette keys. Codon and RNA types reuse AA and DNA
# groupings respectively.

_AA_CLUSTAL_GROUPS = {
    **dict.fromkeys('AILMFWVCailmfwvc', 'hydrophobic'),
    **dict.fromkeys('KRkr',             'positive'),
    **dict.fromkeys('EDed',             'negative'),
    **dict.fromkeys('NQSTnqst',         'polar'),
    **dict.fromkeys('Gg',               'glycine'),
    **dict.fromkeys('Pp',               'proline'),
    **dict.fromkeys('HYhy',             'aromatic'),
    **dict.fromkeys('Xx',               'unknown'),
}

# Note: the original Taylor scheme omits lowercase m, f, v, c — preserved here.
_AA_TAYLOR_GROUPS = {
    **dict.fromkeys('AILMFVCailMFVC', 'hydrophobic'),
    **dict.fromkeys('KRHkrh',         'basic'),
    **dict.fromkeys('EDSTedst',       'polar_charged'),
    **dict.fromkeys('NQnq',           'polar'),
    **dict.fromkeys('CGPcgp',         'special'),
    **dict.fromkeys('Yy',             'aromatic'),
    **dict.fromkeys('X',              'unknown'),
}

_AA_CINEMA_GROUPS = {
    **dict.fromkeys('HKRhkr',           'positive'),
    **dict.fromkeys('DEde',             'negative'),
    **dict.fromkeys('STNQstnq',         'polar'),
    **dict.fromkeys('AVILMFWYavilmfwy', 'hydrophobic'),
    **dict.fromkeys('Gg',               'glycine'),
    **dict.fromkeys('Pp',               'proline'),
    **dict.fromkeys('Cc',               'cysteine'),
    **dict.fromkeys('Xx',               'unknown'),
}

_AA_ZAPPO_GROUPS = {
    **dict.fromkeys('ILVAMilvam',       'aliphatic'),
    **dict.fromkeys('FWYfwy',           'aromatic'),
    **dict.fromkeys('KRHkrh',           'positive'),
    **dict.fromkeys('DEde',             'negative'),
    **dict.fromkeys('STNQstnq',         'hydrophilic'),
    **dict.fromkeys('GPgp',             'special'),
    **dict.fromkeys('Cc',               'cysteine'),
    **dict.fromkeys('Xx',               'unknown'),
}

_AA_HYDROPHOBICITY_GROUPS = {
    **dict.fromkeys('AILMFVPGailmfvpg', 'hydrophobic'),
    **dict.fromkeys('QNHSTYCWqnhstycw', 'polar'),
    **dict.fromkeys('RKDErkde',         'charged'),
    **dict.fromkeys('X',                'unknown'),
}


# UCSC Genome Browser nucleotide colors.
_DNA_UCSC_PALETTE = {
    'A': (_rgb_bg(  0, 150,   0) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# green → white text
    'C': (_rgb_bg(  0,   0, 255) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# blue → white text
    'G': (_rgb_bg(255, 170,   0),                          Back.WHITE),             # orange
    'T': (_rgb_bg(255,   0,   0) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# red → white text
    'U': (_rgb_bg(255,   0,   0) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),
    'N': (Back.WHITE,                                      Back.WHITE),
}

# Okabe-Ito colorblind-safe palette (Okabe & Ito 2008, Color Universal Design).
_DNA_OKABEITO_PALETTE = {
    'A': (_rgb_bg( 86, 180, 233),                          Back.WHITE),             # #56B4E9 sky blue
    'C': (_rgb_bg(240, 228,  66),                          Back.WHITE),             # #F0E442 yellow
    'G': (_rgb_bg(  0, 158, 115) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# #009E73 bluish green
    'T': (_rgb_bg(230, 159,   0),                          Back.WHITE),             # #E69F00 orange
    'U': (_rgb_bg(230, 159,   0),                          Back.WHITE),
    'N': (Back.WHITE,                                      Back.WHITE),
}

# Hello Kitty: hot pink, red bow, gold, white body.
_DNA_HELLOKITTY_PALETTE = {
    'A': (_rgb_bg(255, 105, 180),                          Back.WHITE),             # #FF69B4 hot pink
    'C': (_rgb_bg(255, 255, 255),                          Back.WHITE),             # #FFFFFF white
    'G': (_rgb_bg(255, 215,   0),                          Back.WHITE),             # #FFD700 gold
    'T': (_rgb_bg(255,   0,   0) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),# #FF0000 red → white text
    'U': (_rgb_bg(255,   0,   0) + _rgb_fg(255, 255, 255), Back.WHITE + Fore.BLACK),
    'N': (Back.WHITE,                                      Back.WHITE),
}

_DNA_STOCKHOLM_PALETTE = {
    'A': (_rgb_bg(0, 47, 95)    + _rgb_fg(218, 218, 218), Back.WHITE + Fore.BLACK),
    'C': (_rgb_bg(235, 113, 37),  Back.WHITE),
    'G': (_rgb_bg(163, 168, 107), Back.WHITE),
    'T': (_rgb_bg(172, 222, 230), Back.WHITE),
    'U': (_rgb_bg(172, 222, 230), Back.WHITE),
    'N': (Back.WHITE,             Back.WHITE),
}


# --- Color scheme registry ---
# Maps scheme name -> seqtype -> (grouping, palette).
# DNA palettes use letter keys directly; grouping is None (identity lookup).
# Codon alignments delegate to AA schemes; RNA alignments delegate to DNA schemes.

COLOR_SCHEMES = {
    'clustal': {
        'aa':  (_AA_CLUSTAL_GROUPS, _AA_CLUSTAL_PALETTE),
        'dna': (None,               _DNA_CLUSTAL_PALETTE),
    },
    'terminal': {
        'aa':  (_AA_CLUSTAL_GROUPS, _AA_TERMINAL_PALETTE),
        'dna': (None,               _DNA_TERMINAL_PALETTE),
    },
    'taylor': {
        'aa':  (_AA_TAYLOR_GROUPS,  _AA_TAYLOR_PALETTE),
    },
    'cinema': {
        'aa':  (_AA_CINEMA_GROUPS, _AA_CINEMA_PALETTE),
    },
    'zappo': {
        'aa':  (_AA_ZAPPO_GROUPS, _AA_ZAPPO_PALETTE),
    },
    'hydrophobicity': {
        'aa':  (_AA_HYDROPHOBICITY_GROUPS, _AA_HYDROPHOBICITY_PALETTE),
    },
    'ucsc': {
        'dna': (None, _DNA_UCSC_PALETTE),
    },
    'okabe-ito': {
        'aa':  (_AA_CLUSTAL_GROUPS, _AA_OKABEITO_PALETTE),
        'dna': (None,               _DNA_OKABEITO_PALETTE),
    },
    'hellokitty': {
        'dna': (None, _DNA_HELLOKITTY_PALETTE),
    },
    'stockholm': {
        'dna': (None, _DNA_STOCKHOLM_PALETTE),
    },
}


class Painter:
    '''
    Base class for all painters.
    '''
    def __init__(self):
        self.restrictions = []

    def set_options(self, args):
        if args.keep_colors_when_redirecting:
            colorama_init(strip=False)
        else:
            colorama_init()
        if args.majority:
            self.restrictions.append(restrict_to_majority)
        if args.no_indels:
            self.restrictions.append(restrict_to_no_indels)
        if args.only_variable:
            self.restrictions.append(restrict_to_variable)
        if args.only_variable_excluding_indels:
            self.restrictions.append(restrict_to_variable_excluding_indels)

    def color_for_bad_data(self):
        return Back.WHITE + Fore.RED, Fore.BLACK + Back.WHITE

    def color_for_stop(self):
        return Back.BLACK + Fore.RED, Fore.BLACK

    def indel_color(self):
        return Back.WHITE, Back.WHITE

    def eol(self):
        return Style.RESET_ALL

    def sol(self):
        return Fore.BLACK + Back.WHITE

    def colorizer(self, c, column):
        if c == '!?':
            before_color, after_color = self.color_for_bad_data()
            return before_color + c + after_color
        elif all(map(lambda r: r(column), self.restrictions)):
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

def restrict_to_variable(column):
    return len(column) > 1

def restrict_to_variable_excluding_indels(column):
    if '-' in column:
        return len(column) > 2
    else:
        return len(column) > 1


class AminoAcidPainter(Painter):
    '''
    Painter for amino acid alignments. The color scheme is selected by name.
    '''
    def __init__(self, scheme='clustal'):
        super().__init__()
        self.groups, self.palette = COLOR_SCHEMES[scheme]['aa']

    def _color_lookup(self, c):
        if c in '!?':
            return self.color_for_bad_data()
        elif c == '*':
            return self.color_for_stop()
        elif c in '-_.:':
            return self.indel_color()
        key = self.groups.get(c, c)
        return self.palette.get(key, (Back.WHITE, None))


class DnaPainter(Painter):
    '''
    Painter for DNA/RNA alignments. The color scheme is selected by name.
    '''
    def __init__(self, scheme='clustal'):
        super().__init__()
        _, self.palette = COLOR_SCHEMES[scheme]['dna']

    def _color_lookup(self, c):
        if c in '!*':
            return self.color_for_bad_data()
        elif c in '-.:':
            return self.indel_color()
        return self.palette.get(c.upper(), (Back.WHITE, None))


class CodonPainter(Painter):
    '''
    Painter for codon alignments. Translates each codon to an amino acid
    and delegates to an AminoAcidPainter for coloring.
    '''
    def __init__(self, aa_painter):
        super().__init__()
        self.aa_painter = aa_painter

    def colorizer(self, c, column):
        '''
        c is expected to be a codon, or a single letter.
        Single letters may occur in MACSE alignments, due to frameshifts.
        '''
        try:
            if '!' in c:
                before_color, after_color = self.color_for_bad_data()
                assert after_color # Weird if we don't have a color reset.
                return before_color + c + after_color
            elif all(map(lambda r: r(column), self.restrictions)):
                if len(c) != 3:
                    return Back.WHITE + Fore.RED + Style.BRIGHT + c + Fore.BLACK + Style.NORMAL
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
            return c
