#
# ribolands.parser.barriers
#

import logging
rlog = logging.getLogger(__name__)

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
    ancestor = O('f') + number
    energy   = O('e') + num_flt
    barrier  = O('b') + num_flt

    # saddle option
    saddle = W("~().&")
    senergy   = O('s') + num_flt
    
    #       Lmin[i].my_pool, Lmin[i].ancestors_pool, mfe - kT * log(lmin[i].Z), Lmin[i].my_GradPool, mfe - kT * log(lmin[i].Zg));
    bsize = G(number + number + senergy + number + senergy)

    seq_line = G(T(iupack, 'sequence') + LineEnd().suppress())
    str_line = G(G(T(identity, 'identity')) +
                 G(T(struct, 'structure')) + 
                 G(T(energy, 'energy')) + 
                 G(T(ancestor, 'ancestor')) + 
                 G(T(barrier, 'barrier')) + 
                O(G(T(saddle, 'saddle'))) + 
                O(G(T(bsize, 'bsize'))) + 
                 LineEnd().suppress())

    stmt = seq_line | str_line | LineEnd().suppress()

    document = StringStart() + ZeroOrMore(LineEnd().suppress()) + OneOrMore(stmt) + StringEnd()
    document.ignore(pythonStyleComment)

    return document

def barriers_pathfile_grammar():
    """
    ....(((.((.(((....))))).))).......((((........)))).......... (-10.76) L0005
    ....(((.(..(((....))).).))).......((((........)))).......... ( -7.65) S
    ....(((.(.((((....))))).))).......((((........)))).......... (-10.76) L0004
    ....(((.(.((((....))))).)))...(...((((........))))...)...... ( -7.27) S
    ....(((.(.((((....))))).)))..((...((((........))))...))..... ( -9.12) I
    ....(((.(.((((....))))).))).(((...((((........))))...))).... ( -9.82) I
    ....(((.(.((((....))))).)))((((...((((........))))...))))... (-10.94) L0001
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

    struct = W("().&")
    number = W(nums, nums)
    num_flt = C(O(L('-') | L('+')) + number + O(L('.') + number))
    num_sci = C(O(L('-') | L('+')) + number + O(L('.') + number) + L('e') + O(L('-') | L('+')) + W(nums))
    gorf = num_sci | num_flt

    energy   = S(L('(')) + num_flt + S(L(')'))
    lmin = C(L('L') + number)
    ident = lmin | L('S') | L('I')

    document = StringStart() + ZeroOrMore(LineEnd().suppress()) + OneOrMore(G(struct + energy + ident + LineEnd().suppress())) + StringEnd()
    document.ignore(pythonStyleComment)

    return document

def as_tuple(data):
    seq = None
    lms = []
    LocalMinimum = namedtuple('LocalMinimum', 'id structure energy ancestor barrier saddle pool_size pool_G grad_size grad_G')
    for e, line in enumerate(data):
        if e == 0:
            assert line[0] == 'sequence'
            seq = line[1]
            lms.append(seq)
        else:
            id = None
            ss = None
            en = None
            fa = None
            dG = None
            saddle = None
            pool_size = None
            fpool_size = None
            pool_G = None
            grad_size = None
            grad_G = None
            for [tag, val] in line:
                if tag == 'identity':
                    id = int(val)
                elif tag == 'structure':
                    ss = val
                elif tag == 'energy':
                    en = float(val)
                elif tag == 'ancestor':
                    fa = int(val)
                elif tag == 'barrier':
                    dG = float(val)
                elif tag == 'saddle':
                    saddle = val
                elif tag == 'bsize':
                    pool_size = int(val[0])
                    fpool_size = int(val[1])
                    pool_G = float(val[2])
                    grad_size = int(val[3])
                    grad_G = float(val[4])
                else:
                    raise NotImplementedError('unidentified input')
            lm = LocalMinimum(id, ss, en, fa, dG, saddle, pool_size, pool_G, grad_size, grad_G)
            lms.append(lm)
    return lms

def parse_barriers_file(data):
    document = barriers_grammar()
    return document.parseFile(data).asList()

def parse_barriers_string(data):
    document = barriers_grammar()
    return document.parseString(data).asList()

def parse_barriers(data, is_file = True, return_tuple = False):
    list_data = parse_barriers_file(data) if is_file else parse_barriers_string(data)
    return as_tuple(list_data) if return_tuple else list_data

def parse_barriers_pathfile(data):
    document = barriers_pathfile_grammar()
    return document.parseFile(data).asList()

def parse_barriers_pathstring(data):
    document = barriers_pathfile_grammar()
    return document.parseString(data).asList()

def parse_barriers_path(data, is_file = True):
    return parse_barriers_pathfile(data) if is_file else parse_barriers_pathstring(data)

