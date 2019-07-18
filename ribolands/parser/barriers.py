#
# ribolands.parser.barriers
#   - copy and/or modify together with tests/test_barriers_parser.py
#
# Written by Stefan Badelt (stefan.badelt@gmail.com)
#
# Distributed under the MIT License, use at your own risk.
#

from collections import namedtuple
from pyparsing import (Dict, Word, Literal, Group, Suppress, Optional, ZeroOrMore,
        Combine, White, OneOrMore, alphas, alphanums, nums, delimitedList,
        StringStart, StringEnd, Forward, LineEnd, pythonStyleComment,
        ParseElementEnhance)


def barriers_grammar():
    """ A grammar to interpret output from the program barriers.

    https://github.com/ViennaRNA/Barriers

    tested version: 0.7.1

    # A random barriers example:
         CCCCGGAAGGAAUGGGUAGUGACAGCAGCUUAGGUCGCUGCAUCAUCCCC
       1 ........((.((((((((((((..........)))))))).)))).)). -15.10    0  15.00
       2 ..((....))...((((.(((.((((.((....)).))))))).)))).. -15.00    1  13.50
       3 ..((....))...(((..((((..(((((.......))))).))))))). -14.90    2  11.90
       4 ........(((.(((((((((((..........)))))))).)))))).. -14.90    1   3.30
       5 ..((....))...(((..((((..(((((.......))))).)))).))) -14.40    3   5.60
       6 ..((....))...(((..(((.((((.((....)).)))))))...))). -14.30    2   6.20
       7 ..((....)).((((((((((((..........)))))))).)))).... -13.70    1   3.30
       8 .((.....))...((((.(((.((((.((....)).))))))).)))).. -13.60    2   3.30
       9 .((((.......))))..((((..(((((.......))))).)))).... -13.50    3   7.20
      10 .((.....))...(((..((((..(((((.......))))).))))))). -13.50    3   3.30
      11 ..((....))...(((..((((((((.((....)).))))..))))))). -13.50    6   4.20
      12 .((.....))...(((..((((..(((((.......))))).)))).))) -13.00    5   3.30
      13 ..((....))...(((..((((((((.((....)).))))..)))).))) -13.00    2   5.60
      14 .((.....))...(((..(((.((((.((....)).)))))))...))). -12.90    6   3.30
    """

    # Remove line breaks from DEFAULT_WHITE_CHARS, they are important.
    DWC = "".join([x for x in ParseElementEnhance.DEFAULT_WHITE_CHARS if x != "\n"])
    ParseElementEnhance.setDefaultWhitespaceChars(DWC)

    def T(x, tag):
        def TPA(tag):
            return lambda s, l, t: [tag] + t.asList()
        return x.setParseAction(TPA(tag))

    W = Word
    G = Group
    S = Suppress
    O = Optional
    C = Combine
    L = Literal
    D = Dict

    iupack = W("ACGTURYSMWKVHDBN")
    struct = W("().&")
    number = W(nums, nums)
    num_flt = C(O(L('-') | L('+')) + number + O(L('.') + number))
    num_sci = C(O(L('-') | L('+')) + number + O(L('.') + number) + L('e') + O(L('-') | L('+')) + W(nums))
    gorf = num_sci | num_flt

    # HACK
    identity = O('i') + number
    father   = O('f') + number
    energy   = O('e') + num_flt
    barrier  = O('b') + num_flt

    seq_line = G(T(iupack, 'sequence') + LineEnd().suppress())
    str_line = G(G(T(identity, 'identity')) + G(T(struct, 'structure')) + G(T(energy, 'energy')) + G(T(father, 'father')) + G(T(barrier, 'barrier')) + LineEnd().suppress())

    stmt = seq_line | str_line 

    document = StringStart() + ZeroOrMore(LineEnd().suppress()) + OneOrMore(stmt) + StringEnd()
    document.ignore(pythonStyleComment)

    return document

def as_tuple(data):
    seq = None
    lms = []
    LocalMinimum = namedtuple('LocalMinimum', 'id structure energy father barrier')
    for e, line in enumerate(data):
        if e == 0:
            assert line[0] == 'sequence'
            seq = line[1]
        else:
            id = None
            ss = None
            en = None
            fa = None
            dG = None
            for [tag, val] in line:
                if tag == 'identity':
                    id = int(val)
                elif tag == 'structure':
                    ss = val
                else:
                    return Exception('unidentified input')


def parse_barriers_file(data):
    document = barriers_grammar()
    return document.parseFile(data).asList()

def parse_barriers_string(data):
    document = barriers_grammar()
    return document.parseString(data).asList()

def parse_barriers(data, is_file = True):
    return parse_barriers_file(data) if is_file else parse_barriers_string(data)

