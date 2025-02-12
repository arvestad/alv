import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import alv

class TestInitFunctions(unittest.TestCase):
    def setUp(self):
        # Create test alignment
        self.test_alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ARNDCEQHIKL"), annotations={'molecule_type':'protein'}, id="Alpha"),
            SeqRecord(Seq("ARNDCEQHIKL"), annotations={'molecule_type':'protein'}, id="Beta"),
            SeqRecord(Seq("ARNDCEQHIKM"), annotations={'molecule_type':'protein'}, id="Gamma"),
        ])

    def test_view_basic(self):
        # Test basic functionality without errors
        try:
            alv.view(self.test_alignment)
        except Exception as e:
            self.fail(f"view() raised unexpected exception: {e}")

    def test_view_parameters(self):
        # Test view with different parameters
        test_cases = [
            {'seqtype': 'aa', 'width': 60, 'dotted': True},
            {'seqtype': 'guess', 'width': 80, 'dotted': False},
            {'seqtype': 'aa', 'color_scheme': 'taylor', 'width': 40},
            {'al_start': 0, 'al_end': 5, 'width': 60}
        ]

        for params in test_cases:
            try:
                alv.view(self.test_alignment, **params)
            except Exception as e:
                self.fail(f"view() with params {params} raised unexpected exception: {e}")

    def test_glimpse_basic(self):
        # Test basic functionality without errors
        try:
            alv.glimpse(self.test_alignment)
        except Exception as e:
            self.fail(f"glimpse() raised unexpected exception: {e}")

    def test_glimpse_parameters(self):
        # Test glimpse with different parameters
        test_cases = [
            {'n_seq': 2, 'width': 60, 'dotted': True},
            {'n_seq': 3, 'seqtype': 'aa', 'width': 40},
            {'n_seq': 1, 'color_scheme': 'taylor'},
            {'n_seq': 2, 'al_start': 0, 'al_end': 5}
        ]

        for params in test_cases:
            try:
                alv.glimpse(self.test_alignment, **params)
            except Exception as e:
                self.fail(f"glimpse() with params {params} raised unexpected exception: {e}")

if __name__ == '__main__':
    unittest.main()

