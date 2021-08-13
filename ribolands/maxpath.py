
import findpath as maxpath

def init_findpath_max(sequence, md = None, mp = False):
    assert md is None
    model_details = {
        "noLP": 0,
        "logML": 0,
        "temperature": 37.0,
        "dangles": 2,
        "special_hp": 1,
        "noGU": 0,
        "noGUclosure": 0,
    }
    return maxpath.findpath_class(sequence, mp=mp, model_details=model_details)

def findpath_max(fp, s1, s2, w = 4):
    fp.init(s1, s2, w)
    saddle_energy = fp.get_en()
    path = fp.get_path()
    sspath = []
    ss = s1
    for (i, j, en) in path:
        if 0 < i:
            assert i < j 
            assert ss[i-1] == '.'
            assert ss[j-1] == '.'
            ss = ss[:i-1] + '(' + ss[i:j-1] + ')' + ss[j:]
        if 0 > i:
            assert i > j 
            assert ss[-i-1] == '('
            assert ss[-j-1] == ')'
            ss = ss[:-i-1] + '.' + ss[-i:-j-1] + '.' + ss[-j:]
        sspath.append((ss, en))
    return saddle_energy, sspath

