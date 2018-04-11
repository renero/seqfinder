"""seq_align
Align two sequences to the point where max match is found.
"""

__version__ = '0.1'
__author__ = 'J.Renero'

import sys


def align(a, b):
    best, best_x = 0, 0
    for x in range(len(a)):
        b_substring = b if x is 0 else b[:-x]
        s = sum(i == j for (i, j) in zip(a[x:], b_substring))
        print(a[x:], ' <?> ', b_substring, ' = ', s, sep='')
        if s > best:
            best, best_x = s, x
    print('--')
    for x in range(len(a)):
        a_substring = a if x is 0 else a[:-x]
        s = sum(i == j for (i, j) in zip(a_substring, b[x:]))
        print(a_substring, ' <?> ', b[x:], ' = ', s, sep='')
        if s > best:
            best, best_x = s, x
    print('--')
    for x in range(len(a)):
        a_substring = a if x is 0 else a[:-x]
        b_substring = b if x is 0 else b[:-x]
        s = sum(i == j for (i, j) in zip(a_substring, b_substring))
        print(a_substring, ' <?> ', b_substring, ' = ', s, sep='')
        if s > best:
            best, best_x = s, x
    print('--')
    return best_x


if len(sys.argv) is not 3:
    print(len(sys.argv))
    print("Usage: {} pattern_A pattern_B".format(sys.argv[0]))
    exit(1)

pat_A, pat_B = sys.argv[1], sys.argv[2]
pos = align(pat_B, pat_A)
print(pat_A)
print(pat_B)
print(pos)
s1 = pat_A + '-' * pos
s2 = '-' * pos + pat_B
print(s1)
print(s2)
