# See: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PunziFom
# and https://arxiv.org/pdf/physics/0308063.pdf

import math 

def compute_S_min (a,b, B):
    first_term = 0.5 *(b**2)
    second_term = a * math.sqrt(B)
    third_term = 0.5 * b * math.sqrt(b**2 + (4*a*math.sqrt(B)) + (4 * B))
    S_min = first_term + second_term + third_term
    return S_min
    
#print compute_S_min(1,1,1) #returns 3, as it should

def compute_S_min_imp(a,b,B):
    first_term = (1./8.)*(a**2)
    second_term = (9./13.)*(b**2)
    third_term = a * math.sqrt(B)
    fourth_term = 0.5 * b * math.sqrt(b**2 + (4*a*math.sqrt(B)) + (4 * B))
    S_min_imp = first_term + second_term + third_term + fourth_term
#    print first_term, second_term, third_term, fourth_term
    return S_min_imp

#print compute_S_min_imp(1,1,1) #returns 3.31730769231, as it should

#More testing
# a = 1.
# b = 1.
# B = 1.
# 
# print compute_S_min (a,b,B)


