from scipy import stats


def KS_test(x,y):
    if x and y:
        less = stats.mstats.ks_2samp(x, y,alternative = 'less')[1] 
        greater = stats.mstats.ks_2samp(x, y,alternative = 'greater')[1]
        if (x == y) or (sum(x) == 0 and sum(y) == 0): #樣本相同 或 兩個樣本皆為0
            two_sided = 1.0
        else:
            two_sided = stats.mstats.ks_2samp(x, y,alternative = 'two-sided')[1]
    else:
        less = greater = two_sided = 0

    return {'two_sided':two_sided, 'greater':less, 'less':greater}
    
def T_test(x,y):
    if x and y:
        d, two_sided = stats.ttest_ind(x, y, equal_var=False)
        if d < 0:
            greater = 1 - two_sided/2 #"greater" is the alternative that x has a larger mean than y
            less = two_sided/2 #"less" is the alternative that x has a larger mean than y
        elif d >0:
            greater = two_sided/2
            less = 1 - two_sided/2
        else:
            greater = less = two_sided/2
    else:
        less = greater = two_sided = 0

    return {'two_sided':two_sided, 'greater':greater, 'less':less}

def U_test(x,y):
    if x and y:
        d, two_sided = stats.ranksums(x, y)
        if d < 0:
            greater = 1 - two_sided/2
            less = two_sided/2
        elif d >0:
            greater = two_sided/2
            less = 1 - two_sided/2
        else:
            greater = less = two_sided/2
    else:
        less = greater = two_sided = 0
        
    return {'two_sided':two_sided, 'greater':greater, 'less':less}

