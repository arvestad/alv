import unittest
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import alv.io

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

if __name__ == '__main__':
    unittest.main()
    
