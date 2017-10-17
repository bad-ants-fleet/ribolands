import RNA
import ribolands as ril



def fold_exterior_loop(md, seq, con, ext_moves, exterior_only=True):
  """ Constrained folding 
  
  The default behavior is "exterior_only", which replaces all constrained
  helices with short 'NNN' stretches at the sequence level. This reduces 
  the sequence length (n) and therefore the runtime O(n^3)
  
  :param seq: RNA sequence
  :param con: constraint
  :param exterior_only: only fold the extior loop region

  :return: secondary structure
  """

  if exterior_only :
    spacer = 'NNN'
    pt = ril.make_pair_table(con, base=0)
    ext = ''
    ext_con = ''

    # shrink the sequences
    skip = 0
    for i, j in enumerate(pt):
      if i < skip : continue
      if j == -1 : 
        ext += seq[i]
        ext_con += '.'
      else :
        ext += seq[i]
        ext += seq[i+1]
        ext += spacer
        ext += seq[j-1]
        ext += seq[j]
        ext_con += '(('
        ext_con += 'x' * len(spacer)
        ext_con += '))'
        skip = j+1

    #print
    #print 'el', ext
    #print 'ec', ext_con

    # If we have seen this exterior loop before, then we don't need to 
    # calculate again, and we have to trace back if the parents are connected.
    if ext in ext_moves :
      css = ext_moves[ext][1]
    else :
      fc_tmp = RNA.fold_compound(ext, md)
      fc_tmp.constraints_add(ext_con, RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP)
      css, cfe = fc_tmp.mfe()
      #css, cfe = RNA.fold(ext)
      ext_moves[ext] = [set(), css]
      del fc_tmp
    
    # replace characters in constraint
    c, skip = 0, 0
    for i, j in enumerate(pt):
      if i < skip : continue
      if j == -1 : 
        con = con[:i] + css[c] + con[i+1:]
        c += 1
      else :
        c += len(spacer) + 4
        skip = j+1
    ss = con
    
  else :
    raise DeprecationWarning('This mode is deprecated, it ignores stuff')
    # Force copy of string for ViennaRNA swig interface bug
    tmp = (con + '.')[:-1]
    RNA.cvar.fold_constrained = 1
    ss, mfe = RNA.fold(seq, tmp)
    RNA.cvar.fold_constrained = 0

  return ss, ext

