#author Emanuele Bosi

from pyparsing import nestedExpr
from sys import argv
 
def listToBoolean(l,boolDict={'False':'False','True':'True','or':'or','and':'and'}):
    l_=map(lambda x: boolDict.get(x,'True'),l)
    return ' '.join(l_)
 
def str2Boolean(s,boolDict={'False':'False','True':'True','or':'or','and':'and'}):
    return boolDict.get(s.group(0),'True')
 
def traverse(o, tree_types=(list, tuple)):
    if isinstance(o, tree_types):
        str_,l = '',[]
        for value in o:
            for subvalue in traverse(value, tree_types):
                str_ = '%s %s' %(str_,subvalue)
                l.append(subvalue)
        #print gprSderenator(l)
        yield gprSderenator(l)
    else: yield o
 
def fix(l):
    false_idx = l.index('False')
    if false_idx == 0: return l[2:]
    else: or_idx = false_idx - 1
    return l[:or_idx] + l[false_idx + 1:]
 
def gprSderenator(l,boolDict={'False':'False','True':'True','or':'or','and':'and'}):
            """ main function. If the whole expression is False, then return False.
            If the expression is True, but there are false inside, fix the boolean and return the whole expression."""
            l_=map(lambda x: boolDict.get(x,'True'),l)
            if eval(listToBoolean(l_)):
                        while 'False' in l: l=fix(l)
                        return '(%s)' % ' '.join(l)
            else: return 'False'
 
 
if __name__ == '__main__':
    s = ' '.join(argv[1:])
    for value in traverse(nestedExpr().parseString(s).asList()[0]): pass
    print value
