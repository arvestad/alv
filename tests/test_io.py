import unittest
import argparse
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import alv.io
from alv.alignment import aaAlignment, dnaAlignment, codonAlignment, BaseAlignment
import alv.colorize
import alv.exceptions

class TestIO_etc(unittest.TestCase):
    def test_guess_seq_type(self):
        
        coding = MultipleSeqAlignment([
             SeqRecord(Seq("ACTGCTAGCTAG", generic_dna), id="Alpha"),
             SeqRecord(Seq("ACT-CTAGCTAG", generic_dna), id="Beta"),
             SeqRecord(Seq("ACTGCTAGDTAG", generic_dna), id="Gamma"),
         ])
        noncoding = MultipleSeqAlignment([
             SeqRecord(Seq("CGTCGTCGTCGT", generic_dna), id="Alpha"),
             SeqRecord(Seq("CGACGACGACGA", generic_dna), id="Beta"),
             SeqRecord(Seq("ATAATAATAATA", generic_dna), id="Gamma"),
         ])
        peptide = MultipleSeqAlignment([
             SeqRecord(Seq("ARNDCEQHIKL*", generic_protein), id="Alpha"),
             SeqRecord(Seq("ARNDCEQHIKL*", generic_protein), id="Beta"),
             SeqRecord(Seq("ARNDCEQHIKLM", generic_protein), id="Gamma"),
         ])
        self.assertEqual(alv.io.guess_seq_type(coding), 'codon')
        self.assertEqual(alv.io.guess_seq_type(noncoding), 'dna')
        self.assertEqual(alv.io.guess_seq_type(peptide), 'aa')


class TestFormats(unittest.TestCase):
    def setUp(self):
        self.dna_filename = 'tests/t1.fa'
        self.aa_filename = 'tests/t4.fa'
        self.sthlm_filename= 'tests/t4.sthlm'
        self.pfam_file = 'tests/PF00005_seed.txt'

    def test_reading_fasta_files(self):
        al, painter = alv.io.read_alignment(self.dna_filename, 'dna', 'fasta', '', 'standard')
        self.assertIsInstance(al, dnaAlignment)
        self.assertIsInstance(painter, alv.colorize.dnaPainter)

        al, painter = alv.io.read_alignment(self.dna_filename, 'aa', 'fasta', 'clustal', 'standard')
        self.assertIsInstance(al, aaAlignment)
        self.assertIsInstance(painter, alv.colorize.aaPainter)

        al, painter = alv.io.read_alignment(self.dna_filename, 'codon', 'fasta', 'clustal', 'standard')
        self.assertIsInstance(al, codonAlignment)
        self.assertIsInstance(painter, alv.colorize.codonPainter)

    def test_reading_stockholm_files(self):
        al, painter = alv.io.read_alignment(self.sthlm_filename, 'aa', 'stockholm', '', 'standard')
        self.assertIsInstance(al, aaAlignment)
        self.assertIsInstance(painter, alv.colorize.aaPainter)
        al, painter = alv.io.read_alignment(self.pfam_file, 'aa', 'stockholm', '', 'standard')
        self.assertIsInstance(al, aaAlignment)
        self.assertIsInstance(painter, alv.colorize.aaPainter)


    def test_reading_wrong_format(self):
        with self.assertRaises(alv.exceptions.AlvPossibleFormatError):
            al, painter = alv.io.read_alignment(self.aa_filename, 'aa', 'clustal', '', 'standard')
        with self.assertRaises(alv.exceptions.AlvPossibleFormatError):
            al, painter = alv.io.read_alignment(self.aa_filename, 'aa', 'phylip', '', 'standard')
        with self.assertRaises(alv.exceptions.AlvPossibleFormatError):
            al, painter = alv.io.read_alignment(self.aa_filename, 'aa', 'stockholm', '', 'standard')
        with self.assertRaises(alv.exceptions.AlvPossibleFormatError):
            al, painter = alv.io.read_alignment(self.sthlm_filename, 'aa', 'fasta', '', 'standard')
        with self.assertRaises(alv.exceptions.AlvPossibleFormatError):
            al, painter = alv.io.read_alignment(self.sthlm_filename, 'aa', 'clustal', '', 'standard')
        with self.assertRaises(alv.exceptions.AlvPossibleFormatError):
            al, painter = alv.io.read_alignment(self.sthlm_filename, 'aa', 'phylip', '', 'standard')
        

    def test_guessing_seq_type(self):
        '''
        Given simple DNA or AA input, is the sequence typ (DNA/AA) guessed correctly?
        '''
        al, painter = alv.io.read_alignment(self.dna_filename, 'guess', 'fasta', '', 'standard')
        self.assertIsInstance(al, dnaAlignment)
        al, painter = alv.io.read_alignment(self.aa_filename, 'guess', 'fasta', '', 'standard')
        self.assertIsInstance(al, aaAlignment)
        
if __name__ == '__main__':
    unittest.main()
    
