import unittest
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import alv.io
from alv.alignment import AminoAcidAlignment, DnaAlignment, CodonAlignment, BaseAlignment
import alv.alignmentterminal as at
import alv.colorize
import alv.exceptions

class TestIOetc(unittest.TestCase):
    def test_guess_seq_type(self):

        coding = MultipleSeqAlignment([
             SeqRecord(Seq("AAGCTGATCAGC"), annotations={'molecule_type':'DNA'}, id="Alpha"),
             SeqRecord(Seq("CTGAAGATCAGC"), annotations={'molecule_type':'DNA'}, id="Beta"),
             SeqRecord(Seq("AAGCTGATCAGG"), annotations={'molecule_type':'DNA'}, id="Gamma"),
         ])
        noncoding = MultipleSeqAlignment([
             SeqRecord(Seq("CGTCGTCGTCGT"), annotations={'molecule_type':'DNA'}, id="Alpha"),
             SeqRecord(Seq("CGACGACGACGA"), annotations={'molecule_type':'DNA'}, id="Beta"),
             SeqRecord(Seq("ATAATAATAATA"), annotations={'molecule_type':'DNA'}, id="Gamma"),
         ])
        peptide = MultipleSeqAlignment([
             SeqRecord(Seq("ARNDCEQHIKL*"), annotations={'molecule_type':'protein'}, id="Alpha"),
             SeqRecord(Seq("ARNDCEQHIKL*"), annotations={'molecule_type':'protein'}, id="Beta"),
             SeqRecord(Seq("ARNDCEQHIKLM"), annotations={'molecule_type':'protein'}, id="Gamma"),
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
        self.nexus_filename = 'tests/test.nex'

    def test_reading_fasta_files(self):
        al, painter = alv.io.read_alignment(self.dna_filename, 'dna', 'fasta', '', 'standard')
        self.assertIsInstance(al, DnaAlignment)
        self.assertIsInstance(painter, alv.colorize.DnaPainter)

        al, painter = alv.io.read_alignment(self.dna_filename, 'aa', 'fasta', 'clustal', 'standard')
        self.assertIsInstance(al, AminoAcidAlignment)
        self.assertIsInstance(painter, alv.colorize.AminoAcidPainter)

        al, painter = alv.io.read_alignment(self.dna_filename, 'codon', 'fasta', 'clustal', 'standard')
        self.assertIsInstance(al, CodonAlignment)
        self.assertIsInstance(painter, alv.colorize.CodonPainter)

    def test_reading_stockholm_files(self):
        al, painter = alv.io.read_alignment(self.sthlm_filename, 'aa', 'stockholm', '', 'standard')
        self.assertIsInstance(al, AminoAcidAlignment)
        self.assertIsInstance(painter, alv.colorize.AminoAcidPainter)
        al, painter = alv.io.read_alignment(self.pfam_file, 'aa', 'stockholm', '', 'standard')
        self.assertIsInstance(al, AminoAcidAlignment)
        self.assertIsInstance(painter, alv.colorize.AminoAcidPainter)

    def test_reading_nexus_files(self):
        al, painter = alv.io.read_alignment(self.nexus_filename, 'aa', 'nexus', '', 'standard')
        self.assertIsInstance(al, AminoAcidAlignment)
        self.assertIsInstance(painter, alv.colorize.AminoAcidPainter)

    def test_reading_wrong_format(self):
        with self.assertRaises(ValueError):
            al, painter = alv.io.read_alignment(self.aa_filename, 'aa', 'clustal', '', 'standard')
        with self.assertRaises(ValueError):
            al, painter = alv.io.read_alignment(self.aa_filename, 'aa', 'phylip', '', 'standard')
        with self.assertRaises(ValueError):
            al, painter = alv.io.read_alignment(self.aa_filename, 'aa', 'stockholm', '', 'standard')
        with self.assertRaises(ValueError):
            al, painter = alv.io.read_alignment(self.aa_filename, 'aa', 'nexus', '', 'standard')
        with self.assertRaises(ValueError):
            al, painter = alv.io.read_alignment(self.sthlm_filename, 'aa', 'fasta', '', 'standard')
        with self.assertRaises(ValueError):
            al, painter = alv.io.read_alignment(self.sthlm_filename, 'aa', 'clustal', '', 'standard')
        with self.assertRaises(ValueError):
            al, painter = alv.io.read_alignment(self.sthlm_filename, 'aa', 'phylip', '', 'standard')
        with self.assertRaises(ValueError):
            al, painter = alv.io.read_alignment(self.sthlm_filename, 'aa', 'nexus', '', 'standard')
        with self.assertRaises(ValueError):
            al, painter = alv.io.read_alignment(self.nexus_filename, 'aa', 'fasta', '', 'standard')
        with self.assertRaises(ValueError):
            al, painter = alv.io.read_alignment(self.nexus_filename, 'aa', 'phylip', '', 'standard')
        with self.assertRaises(ValueError):
            al, painter = alv.io.read_alignment(self.nexus_filename, 'aa', 'clustal', '', 'standard')
        with self.assertRaises(ValueError):
            al, painter = alv.io.read_alignment(self.nexus_filename, 'aa', 'stockholm', '', 'standard')


    def test_guessing_seq_type(self):
        '''
        Given simple DNA or AA input, is the sequence typ (DNA/AA) guessed correctly?
        '''
        al, painter = alv.io.read_alignment(self.dna_filename, 'guess', 'fasta', '', 'standard')
        self.assertIsInstance(al, DnaAlignment)
        al, painter = alv.io.read_alignment(self.aa_filename, 'guess', 'fasta', '', 'standard')
        self.assertIsInstance(al, AminoAcidAlignment)
        al, painter = alv.io.read_alignment(self.nexus_filename, 'guess', 'nexus', '', 'standard')
        self.assertIsInstance(al, CodonAlignment)


class TestBlockCalculation(unittest.TestCase):
    def setUp(self):
        self.al1, _ = alv.io.read_alignment('tests/t2.fa', 'dna', 'fasta','', 'standard')
        self.al2, _ = alv.io.read_alignment('tests/t2_longaccessions.fa', 'dna', 'fasta','', 'standard')


    def test_block_width_calculation(self):
        self.assertEqual(10, self.al1.block_width(100, 0))
        self.assertEqual(10, self.al1.block_width(100, 10))
        self.assertEqual(5, self.al1.block_width(100, 5))
        self.assertEqual(10, self.al1.block_width(10, 10))


class TestIndexBar(unittest.TestCase):
    def setUp(self):
        self.tickmark = '^'
        import sys
        if sys.stdout.encoding == 'UTF-8':
            self.tickmark = '↑'

    def test_make_one_tick(self):
        s = at.make_one_tick(5, 10, '^')
        self.assertEqual(len(s), 10)
        self.assertEqual(s, '        5^')

    def test_calc_tick_indices(self):
        indices = list(at.calc_tick_indices(33, 100, 10, 7))
        self.assertEqual(indices, [40, 50, 60, 70, 80, 90])

        indices = list(at.calc_tick_indices(50, 100, 10,10))
        self.assertEqual(indices, [60, 70, 80, 90])

    def test_make_tick_string(self):
        s = at.make_tick_string(6, 0, 40, 40, 10)
        if self.tickmark == '↑':
            self.assertEqual(s, '     0↑')
        else:
            self.assertEqual(s, '     0^')

        s = at.make_tick_string(6, 0, 40, 20, 10)
        if self.tickmark == '↑':
            self.assertEqual(s, '     0↑                 20↑')
        else:
            self.assertEqual(s, '     0^                 20^')

        s = at.make_tick_string(5, 0, 50, 20, 7)
        if self.tickmark == '↑':
            self.assertEqual(s, '    0↑                 20↑                 40↑')
        else:
            self.assertEqual(s, '    0^                 20^                 40^')

        s = at.make_tick_string(7, 35, 100, 20, 10)
        if self.tickmark == '↑':
            self.assertEqual(s, '     35↑                      60↑                 80↑')
        else:
            self.assertEqual(s, '     35^                      60^                 80^')


if __name__ == '__main__':
    unittest.main()
