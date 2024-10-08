#!/usr/bin/env python3

import unittest

from convert_superstring import split_superstring, get_k


class TestConvertSuperstring(unittest.TestCase):

    def test_split_superstring(self):
        tests = [
            # (superstring, k, want_result)
            ("ATttCgg", 3, ["ATTT", "CGG"]),
            ("ATttCgga", 4, ["ATTTC", "CGGA"]),
            ("ATtGCtCg", 2, ["ATT", "GCT", "CG"]),
            ("ATCGCTtGCTATACtcCACgttta", 6, ["ATCGCTTGCTA", "GCTATACTCCAC", "CACGTTTA"]),
        ]

        for t in tests:
            got_result = split_superstring(t[0], t[1])

            self.assertEqual(t[2], got_result)

    def test_get_k(self):
        tests = [
            # (superstring, want_result)
            ("ATttCgg", 3),
            ("ATttCgga", 4),
            ("ATC", 1),
            ("Aaaaaa", 6),
        ]

        for t in tests:
            got_result = get_k(t[0])

            self.assertEqual(t[1], got_result)


if __name__ == '__main__':
    unittest.main()
