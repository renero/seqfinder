from unittest import TestCase
from patterns import replace, generate_patterns


class TestReplace(TestCase):
    def test_three_letters(self):
        self.assertEqual(replace('ABC', pos=0, reps=2, gene_length=1), '?BC')
        self.assertEqual(replace('ABCD', 2, 2, 2), 'AB??')

    def test_generate_patterns(self):
        self.assertEqual(generate_patterns('ABC', max_regex_length=1, gene_length=1),
                         {'?BC', 'A?C', 'AB?', 'ABC'})

        self.assertEqual(generate_patterns('ABCD', max_regex_length=1, gene_length=2),
                         {'ABCD', '??CD', 'AB??'})

        self.assertEqual(generate_patterns('ABCDEF', max_regex_length=2, gene_length=2),
                         {'????EF', '??CDEF', 'AB????', 'AB??EF', 'ABCD??', 'ABCDEF'})
