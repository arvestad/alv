import unittest
from alv.alignmentterminal import extract_accession_order

class TestOrdering(unittest.TestCase):
    def test_explicit_ordering(self):
        order = extract_accession_order('a,b,c')
        self.assertEqual(len(order), 3)

        with self.assertRaises(ValueError):
            order = extract_accession_order('')

        with self.assertRaises(ValueError):
            order = extract_accession_order('a,')

        with self.assertRaises(ValueError):
            order = extract_accession_order('a,,b')

        
