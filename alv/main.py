import alv
import argparse
import os
import sys

from .version import __version__
from .alignmentterminal import AlignmentShellTerminal
import .io as io
from .exceptions import AlvPossibleFormatError, AlvEmptyAlignment

# Plain text citation
citation = 'Arvestad, (2018). alv: a console-based viewer for molecular sequence alignments. Journal of Open Source Software, 3(31), 955, https://doi.org/10.21105/joss.00955'

# A citation for BibTeX
bibitem = '''@Article{alv2018,
  author = 	 {Lars Arvestad},
  title = 	 {alv: a console-based viewer for molecular sequence alignments},
  journal = 	 {Journal of Open Source Software},
  year = 	 2018,
  volume = 	 3,
  number = 	 31,
  pages = 	 955,
  doi =          {https://doi.org/10.21105/joss.00955}
}'''

# For the --method option
method_text = 'Alignments were viewed using alv (github.com/arvestad/alv).'

def list_genetic_codes():
    additional_help_text='''
The genetic codes are implemented using BioPython, which reflect what is
available at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi .

The following is the list, with numbers, of known genetic codes:
   1. The Standard Code
   2. The Vertebrate Mitochondrial Code
   3. The Yeast Mitochondrial Code
   4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
   5. The Invertebrate Mitochondrial Code
   6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
   9. The Echinoderm and Flatworm Mitochondrial Code
   10. The Euplotid Nuclear Code
   11. The Bacterial, Archaeal and Plant Plastid Code
   12. The Alternative Yeast Nuclear Code
   13. The Ascidian Mitochondrial Code
   14. The Alternative Flatworm Mitochondrial Code
   16. Chlorophycean Mitochondrial Code
   21. Trematode Mitochondrial Code
   22. Scenedesmus obliquus Mitochondrial Code
   23. Thraustochytrium Mitochondrial Code
   24. Pterobranchia Mitochondrial Code
   25. Candidate Division SR1 and Gracilibacteria Code
   26. Pachysolen tannophilus Nuclear Code
   27. Karyorelict Nuclear
   28. Condylostoma Nuclear
   29. Mesodinium Nuclear
   30. Peritrich Nuclear
   31. Blastocrithidia Nuclear
'''
    print(additional_help_text)

def setup_argument_parsing():
    '''
    Create an argument parser, parse and return args.
    '''
    ap = argparse.ArgumentParser()
    ap.add_argument('infile', nargs='?', help="The infile is the path to a file, or '-' if reading from stdin.")
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    ap.add_argument('-ai', '--alignment-index', type=int, default=0,
                    help='If reading file with many alignments, choose which one to output with this zero-based index. Default: first alignment in file.')
    ap.add_argument('-f', '--format', choices=['guess', 'fasta', 'clustal', 'nexus', 'phylip', 'stockholm'], default='guess',
                    help="Specify what sequence type to assume. Be specific if the file is not recognized automatically. When reading from stdin, the format is always guessed to be FASTA. Default: %(default)s")
    ap.add_argument('-t', '--type', choices=['aa', 'dna', 'rna', 'codon', 'guess'], default='guess',
                    help="Specify what sequence type to assume. Coding DNA/RNA is assumed with the 'codon' option. Guessing the format only chooses between 'aa' and 'dna', but assumes the standard genetic code.  Default: %(default)s")
    ap.add_argument('-g', '--glimpse', action='store_true',
                    help='Give a glimpse of an alignment. If the alignment fits without any scrolling and without line breaks, then just view the alignment. Otherwise, identify a conserved part of the MSA and show a random sample of the sequences that fits the screen.')
    ap.add_argument('-c', '--color-scheme', choices=['clustal', 'taylor', 'hydrophobicity'], default='clustal',
                    help='Color scheme for AA and coding DNA/RNA. The clustal coloring scheme is an approximation of the original, due to the limited color choices for consoles. The "hydrophobicity" gives red to hydrophobic, blue to polar, and green to charged residues.  Default: %(default)s')
    ap.add_argument('--code', choices=[1,2,3,4,5,6,9,10,11, 12, 13, 14, 16,21,22,23,24,25,26,27,28,29,30,31], type=int, default=1,
                    help="Genetic code to use, based on NCBI's code list, see details below. Show alternatives with the --list-codes option. Default: %(default)s.")
    ap.add_argument('-d', '--dotted', action='store_true',
                    help="Let the first sequence in output alignment be a template and, for other sequences, show identity to template using a period. Useful for alignments with high similarity.")
    ap.add_argument('-lc', '--list-codes', action='store_true',
                    help="List the available genetic codes and exit.")
    ap.add_argument('-w', '--width', type=int, default=0,
                    help='Width of alignment blocks. Defaults to terminal width minus accession width, essentially.')
    ap.add_argument('-k', '--keep-colors-when-redirecting', action='store_true',
                    help="Do not strip colors when redirecting to stdout, or similar. In particular useful with the command 'less -R'.")
    ap.add_argument('-l', '--pipe-to-less', action='store_true',
                    help="Do not break the alignment into blocks. Implies -k. Suitable when piping to commands like 'less -RS'.")
    # General info
    info_args = ap.add_argument_group('General information')
    info_args.add_argument('-i', '--info', action='store_true',      help="Append basic information about the alignment at the end.")
    info_args.add_argument('-j', '--just-info', action='store_true', help="Write basic information about the alignment and exit.")
    info_args.add_argument('--cite', action='store_true',            help="Write citation example: plain text and a BibTeX item.")
    info_args.add_argument('--method', action='store_true',          help="Write a suggested text to add to a methods section.")

    # Options for changing sequence order
    ordering_args = ap.add_argument_group('Sequence selection and ordering')
    ordering_args.add_argument('-r', '--random-accessions', type=int, metavar='N', default=0,
                               help='Only view a random sample of the alignment sequences.')
    ordering_args.add_argument('-s', '--sorting', choices=['infile', 'alpha'], default='infile',
                               help="Sort the sequences as given in the infile or alphabetically (by accession). Default: %(default)s")
    ordering_args.add_argument('-si', '--sort-by-id', metavar='ACCESSION', type=str,
                               help='Sort the output alignment by similarity (percent identity) to named sequence. Overrides -s.')
    ordering_args.add_argument('-so', '--sorting-order', metavar='ACCESSIONS', type=str,
                               help='Comma-separated list of accessions. Sequences will be presented in this order. Also note that one can choose which sequences to present with this opion. Overrides -s and -si.')
    ordering_args.add_argument('-sm', '--select-matching', metavar='ACCESSION_PATTERN', type=str,
                               help='Only show sequences with accessions containing ACCESSION_PATTERN.')
    ordering_args.add_argument('-sa', '--sub-alignment', nargs=2, metavar='INT', type=int,
                               help='Only show alignment columns given by FROM and UPTO indices.')

    # Options for limiting colorization
    restriction_args = ap.add_argument_group('Restricting colorization')
    restriction_args.add_argument('--majority', action='store_true',
                                  help='Only color those column where the most common amino acid is found in 50 percent of sequences.')
    restriction_args.add_argument('--no-indels', action='store_true',
                                  help='Only color column without indels.')
    restriction_args.add_argument('--only-variable', action='store_true',
                                  help='Only color columns that contain variation.')
    restriction_args.add_argument('--only-variable-excluding-indels', action='store_true',
                                  help='Only color columns that contain variation, ignoring indels.')

    # Options for removing parts of the accession
    accession_args = ap.add_argument_group('Accession trimming')
    accession_args.add_argument('-as', '--acc-substring', nargs=2, metavar='INT', type=int,
                                help="Specify what substring of an accession to keep. '-as 10 15' discards all but position 10 to 14 in any accession.")
    accession_args.add_argument('-aa', '--acc-abbreviate', type=int, metavar='N',
                                help="Keep only the first N and last N characters of the accession")

    return ap
# Feature to add:
#    ap.add_argument('-p', '--prefix', type=int, default=0, help='Number of characters to remove from the beggining of the accession.')


def input_and_option_adaption(args):
    '''
    Read data, and handle some of the program options.
    Return a pair of an alv.Alignment and an alv.Painter instance.
    Exits on error.
    '''
    try:
        if args.format == 'guess' and args.infile != '-':
            format = io.guess_format(args.infile)
        elif args.format == 'guess' and args.infile == '-':
            format = 'fasta'    # Hard guess, because a bit complicated when reading from pipe (sys.stdin)
        else:
            format = args.format
        if args.sub_alignment:  # Are we restricting the alignment?
            start_col, end_col = args.sub_alignment
        else:
            start_col, end_col = 0, -1
        alignment, painter = io.read_alignment(args.infile, args.type, format, args.color_scheme, args.code, start_col, end_col, args.alignment_index)
        return alignment, painter

    except KeyboardInterrupt:
        sys.exit()
    except AlvPossibleFormatError:
        print('alv: cannot guess the format of input.', file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("alv: file '" + args.infile + "' not found.")
        sys.exit(4)
    except ValueError as e:
        # Bio.AlignIO uses ValueError for a number of reading problems
        if str(e) == 'No records found in handle':
            msg = 'probably wrong input format in '+ args.infile+'. Try option -f.'
        elif str(e) == 'Sequences must all be the same length':
            msg = 'unequal sequence lengths. Maybe not aligned input in '+args.infile+'?'
        else:
            msg = str(e)
        print('alv error:', msg, file=sys.stderr)
        sys.exit(2)
    except Exception:
        print('alv bug: Unknown error when reading input.', file=sys.stderr)
        sys.exit(3)


def main():
    ap = setup_argument_parsing()
    args = ap.parse_args()


    # Handle option that do not require an infile
    if args.list_codes:
        list_genetic_codes()
        ap.exit()

    if args.method:
        print(method_text)
        ap.exit()
    if args.cite:
        print(citation)
        print()
        print(bibitem)
        ap.exit()

    # From here on, we need data from an infile
    if not args.infile:
        ap.print_usage()
        ap.exit()

    # Read the data
    try:
        alignment, painter = input_and_option_adaption(args)
    except IOError as e:
        print(f'alv: {e}', file=sys.stderr)
        ap.exit(7)
    except ValueError as e:
        print(f'alv: {e}', file=sys.stderr)
        ap.exit(8)

    # In case we just want to know the basics:
    if args.just_info:
        io.output_al_info(alignment)
        ap.exit()

    # Shorten accessions, if requested
    if args.acc_substring:
        start = args.acc_substring[0]
        stop = args.acc_substring[1]
        if start >= stop or start<0:
            print("alv: bad indices for option '-as'!", file=sys.stderr)
            ap.exit(5)
        alignment.trim_accessions(start, stop)
    elif args.acc_abbreviate:
        alignment.abbreviate_accessions(args.acc_abbreviate)

    # Prepare for output
    if args.pipe_to_less:
        args.width = alignment.al_width()
        args.keep_colors_when_redirecting = True

    terminal = AlignmentShellTerminal(args)
    painter.set_options(args)
    try:
        if args.glimpse:
            terminal.output_glimpse(alignment, painter, args.width, args.dotted)
        else:
            terminal.output_alignment(alignment, painter, args.width, args.dotted)
        if args.info:
            io.output_al_info(alignment)
    except KeyboardInterrupt:
        ap.exit()
    except AlvEmptyAlignment:
        print('alv: input contains no sequence data?', file=sys.stderr)
        ap.exit(4)
    except BrokenPipeError:
        # This should not cause any specific error at all: correct behaviour is to end the program.
        sys.stderr.close()  # Without this line, python will sometimes give an error
        ap.exit(6)
    except IOError as e:
        print(f'alv: ', file=sys.stderr)
        ap.exit(7)
    except ValueError as e:
        print('alv:', e, file=sys.stderr)
        ap.exit(3)
    except Exception as e:
        print('Alv bug! Please report!', file=sys.stderr)
        print(e, file=sys.stderr)
        ap.exit(2)

if __name__ == '__main__':
    main()
