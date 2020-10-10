#!/usr/bin/python


def sort_by_chromosome(chr_list):
    a = [x.split('chr')[-1] for x in chr_list]
    if 'X' in a: 
        X_present = True
        a.remove('X')
    else:
        X_present = False
        
    if 'Y' in a: 
        Y_present = True
        a.remove('Y')
    else:
        Y_present = False
    if 'M' in a: 
        M_present = True
        a.remove('M')
    else:
        M_present = False
    if 'MT' in a:
        MT_present = True
        a.remove('MT')
    else:
        MT_present = False
        
    a = list(map(int,a))
    b = sorted(a)
    c = ['chr'+str(x) for x in b]
    if X_present:
        c.append('chrX')
    if Y_present:
        c.append('chrY')
    if M_present:
        c.append('chrM')
    elif MT_present:
        c.append('chrMT')
        
    return(c)