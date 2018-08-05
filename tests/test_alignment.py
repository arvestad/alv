import unittest
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from alv.alignment import aaAlignment

class TestAlignment(unittest.TestCase):
    def test_summarize_columns(self):
        simple_al = MultipleSeqAlignment([
             SeqRecord(Seq("AR", generic_protein), id="Alpha"),
             SeqRecord(Seq("AR", generic_protein), id="Beta"),
             SeqRecord(Seq("AS", generic_protein), id="Gamma"),
         ])
        al = aaAlignment(simple_al)
        self.assertEqual(al._summarize_columns()[0]['A'], 3)
        self.assertEqual(al._summarize_columns()[1]['R'], 2)
        self.assertEqual(al._summarize_columns()[1]['S'], 1)
        self.assertEqual(al._summarize_columns()[1]['A'], 0)


if __name__ == '__main__':
    unittest.main()
    

