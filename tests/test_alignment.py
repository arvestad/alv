import unittest
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from alv.alignment import AminoAcidAlignment

class TestAlignment(unittest.TestCase):
    def setUp(self):
        simple_al = MultipleSeqAlignment([
             SeqRecord(Seq("AR", generic_protein), id="Alpha"),
             SeqRecord(Seq("AR", generic_protein), id="Beta_"),
             SeqRecord(Seq("AS", generic_protein), id="Gamma"),
         ])
        self.al = AminoAcidAlignment(simple_al)

    def test_summarize_columns(self):
        self.assertEqual(self.al._summarize_columns()[0]['A'], 3)
        self.assertEqual(self.al._summarize_columns()[1]['R'], 2)
        self.assertEqual(self.al._summarize_columns()[1]['S'], 1)
        self.assertEqual(self.al._summarize_columns()[1]['A'], 0)

    def test_accession_length(self):
        for acc in self.al.accessions():
            self.assertEqual(len(acc), 5)

    def test_basic_info(self):
        info = list(self.al.get_basic_info())
        self.assertEqual(info[0][1], 3)
        self.assertEqual(info[1][1], 2)

    def test_column_conservation(self):
        conservation = self.al.get_column_conservation()
        self.assertEqual(conservation[0], 1.0)
        self.assertEqual(conservation[1], 2/3)

if __name__ == '__main__':
    unittest.main()
