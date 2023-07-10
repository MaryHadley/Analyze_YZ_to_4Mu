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


#print compute_S_min(a=1, B=1, b=2) #returns the same thing as compute_S_min(a=1, b=2, B=2), as it should

#Stat Comm suggests a=2, b=5, so we will use those
#norm_factor_<N> == avg_eff_from_YZ_ZY_weighted_avgs_<N>, comes from
# /Users/maryhadley/analysisNotes/compute_weighted_average_using_Run2_lumi.py

#Upsi mu pT cut of 2 case
#a =2 
#b = 5
#B = 159 + 68. 4 = 227.4
#B = 159.353 + 68.4048  = 227.7578
#norm_factor_2 = 0.473110716857

print "Upsi mu pT cut of 2 case:"
print "########################"
print "S_min:  ", compute_S_min(a=2, b=5, B=227.7578)
print "S_min_imp:  ", compute_S_min_imp(a=2, b=5, B=227.7578)
S_min_2 = compute_S_min(a=2, b=5, B=227.7578)
S_min_imp_2 = compute_S_min_imp(a=2, b=5, B=227.7578)
norm_factor_2 = 0.473110716857
print "Normalized S_min:  ", S_min_2 * (1./norm_factor_2)
print "Normalized S_min_imp:  ", S_min_imp_2 * (1./norm_factor_2)
print "########################"

#Upsi mu pT cut of 3 case
#a =2
#b = 5
#B = 88.9 + 48.1 = 137
#B = 88.9234 +48.104 = 137.0274
#norm_factor_3 = 0.44361221096

print "Upsi mu pT cut of 3 case:"
print "########################"
print "S_min:  ", compute_S_min(a=2, b=5, B=137.0274)
print "S_min_imp:  ", compute_S_min_imp(a=2, b=5, B=137.0274)
norm_factor_3 = 0.44361221096
S_min_3 = compute_S_min(a=2, b=5, B=137.0274)
S_min_imp_3 = compute_S_min_imp(a=2, b=5, B=137.0274)
print "Normalized S_min:  ", S_min_3 * (1./norm_factor_3)
print "Normalized S_min_imp:  ", S_min_imp_3 * (1./norm_factor_3)
print "#######################"

#Upsi mu pT cut of 4 case
#a = 2
#b = 5
#B = 42.4 + 25.3 = 67.7
#B = 42.4245 + 25.3029 = 67.7274
#norm_factor_4 = 0.290611678025

print "Upsi mu pT cut of 4 case:"
print "########################"
print "S_min:  ", compute_S_min(a=2, b=5, B = 67.7274)
print "S_min_imp:  ", compute_S_min_imp(a=2, b=5, B = 67.7274)
S_min_4 = compute_S_min(a=2, b=5, B = 67.7274)
S_min_imp_4 = compute_S_min_imp(a=2, b=5, B = 67.7274)
norm_factor_4 = 0.290611678025
print "Normalized S_min:  ", S_min_4 * (1./norm_factor_4)
print "Normalized S_min_imp:  ", S_min_imp_4 * (1./norm_factor_4)
print "########################"

print "Minimum normalized S_min score:  ", min(S_min_2 * (1./norm_factor_2), S_min_3 * (1./norm_factor_3), S_min_4 * (1./norm_factor_4))
print "Minimum normalized S_min_imp score:  ", min(S_min_imp_2 * (1./norm_factor_2), S_min_imp_3 * (1./norm_factor_3), S_min_imp_4 * (1./norm_factor_4))