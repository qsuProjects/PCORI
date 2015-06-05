
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Project: PCORI Missing Data
# 
# -This program creates covariates for simulated data based on existing VA data 
#    (summary stats provided by Vilija at VA, see Q:\Datasets\PCORI\data simulation\real data estimates)
#
# -Using jointly_generate_binary_normal (binNor package) function by Demirtas, simulate clinical characterisitics and drug exposure for
#  a large number of patients over time.
#
# -This file initializes needed variables. 
#
# population means and standard deviations
# probability of being on drug currently given ever-user
# 1. pd_abac = 0.275
# 2. pd_ataz = 0.222
# 3. pd_dida = 0.198
# 4. pd_efav = 0.5
# 5. pd_emtr = 0.413
# 6. pd_indi = 0.206
# 7. pd_lami = 0.688
# 8. pd_lopi = 0.253 
# 9. pd_nelf = 0.24
# 10. pd_nevi = 0.184
# 11. pd_rito = 0.180
# 12. pd_saqu = 0.103
# 13. pd_stav = 0.332
# 14. pd_teno = 0.526
# 15. pd_zido = 0.544
# 16. p1 (gender) = 0.97

# 17. age ~ N(46, 10.12)
# 18. bmi ~ N(25, 4.2)
# 19. cd4 ~ N(308.1, 50)
# 20. vln ~ N(105337, 10000)
# 21. bps ~ N(126.85, 15.48)
# 22. bpd ~ N(76.7, 9.96)
# 23. ldl ~ N(99.1, 23)
# 24. hdl ~ N(39.02, 7)
# 25. trig ~ N(180.87, 30)

# mean proportion of time on drug
# 26. d_abac ~ N(0.12, sqrt(p(1-p)/n))
# 27. d_ataz ~ N(0.08, sqrt(p(1-p)/n))
# 28. d_dida ~ N(0.09, sqrt(p(1-p)/n))
# 29. d_efav ~ N(0.16, sqrt(p(1-p)/n))
# 30. d_emtr ~ N(0.11, sqrt(p(1-p)/n))
# 31. d_indi ~ N(0.07, sqrt(p(1-p)/n))
# 32. d_lami ~ N(0.40, sqrt(p(1-p)/n))
# 33. d_lopi ~ N(0.12, sqrt(p(1-p)/n))
# 34. d_nelf ~ N(0.08, sqrt(p(1-p)/n))
# 35. d_nevi ~ N(0.07, sqrt(p(1-p)/n))
# 36. d_rito ~ N(0.05, sqrt(p(1-p)/n))
# 37. d_saqu ~ N(0.04, sqrt(p(1-p)/n))
# 38. d_stav ~ N(0.17, sqrt(p(1-p)/n))
# 39. d_teno ~ N(0.22, sqrt(p(1-p)/n))
# 40. d_zido ~ N(0.24, sqrt(p(1-p)/n))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



######################### SET GLOBAL VARIABLES #########################

# these will be constant for all simulations
#n.Drugs = 15  # number of drugs
#n.OtherBins = 1
#n.OtherNorms = 9
race.names = c("white", "black", "other")

zero = 0.000000000000000000000000001
one = 0.99



######################### LOAD MATRICES AND PARAMETER DATAFRAME #########################

setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/GENCOV/vignettes")

# within-subject correlation matrix
wcorin = read.csv("ex1_wcor.csv", header=FALSE)[-1,-1]
wcor = as.numeric(t(wcorin)[lower.tri(wcorin)]) #it gets read in by rows

# population correlation matrix
pcorin = read.csv("ex1_pcor.csv", header=FALSE)[-1,-1]
pcor = as.numeric(t(pcorin)[lower.tri(pcorin, diag=F)]) # need to transpose and read in the lower half to convert the a matrix into vector by row

# read in and complete parameters dataframe
parameters = complete_parameters( read.csv("ex1_parameters.csv"), .n.Subj )

# parameters for categorical variable
cat.parameters = read.csv("ex1_categorical_parameters.csv")


######################### RACE MODEL COEFFICIENTS #########################

# relationship between covariates and race (0=white, 1=black, 2=other)
# ln(P(race=1)/P(race=0)) = bb0 + bb1*abacavir + ... + bb37*trig
# ln(P(race=2)/P(race=0)) = bo0 + bo1*abacavir + ... + bo37*trig
# estimates for black race
bb0 = -1.1125  #intercept
bb1 = -0.118   #abacavir
bb2 = -0.1128 #atazanavir
bb3 = -0.1242  #didanosine
bb4 = -0.0917  #efavirenz
bb5 = 0.0255   #emtricitabine
bb6 = -0.171   #indinavir
bb7 = -0.1092  #lamivudine
bb8 = 0.00497  #lopinavir
bb9 = -0.0815  #nalfinavir
bb10 = -0.2936 #nevirapine
bb11 = -0.0885 #ritonavir
bb12 = -0.3105 #saquinavir
bb13 = -0.0398 #stavudine
bb14 = -0.2657 #temofovir
bb15 = 0.1424  #zidovudine
bb16 = -0.3929 #male
bb17 = 0.1421  #age 35-44
bb18 = 0.3382  #age 45-54
bb19 = 0.1168  #age 55-64
bb20 = -0.1656 #age 65+
bb21 = 0.2922  #vl 400-<3599
bb22 = 0.325   #vl 3500-<10K
bb23 = 0.5946  #10k-<50k
bb24 = 0.5069  #>=50k
bb25 = 0.8204  #cd4 <50
bb26 = 0.4093  #cd4 50-<100
bb27 = 0.2606  #cd4 100-<200
bb28 = 0.2337  #cd4 200-<350
bb29 = 0.0707  #cd4 350-<500
bb30 = 0.1884  #bmi <20
bb31 = 0.0369  #bmi 25-<30
bb32 = 0.3407  #bmi >=30
bb33 = 0.019    #bpd
bb34 = -0.00228 #bps
bb35 = 0.0242   #hdl
bb36 = -0.00411 #ldl
bb37 = -0.00174 #trig

# estimates for other race
bo0 = -3.2279 #intercept
bo1 = 0.0475
bo2 = -0.1799
bo3 = 0.025
bo4 = -0.1504
bo5 = 0.1564
bo6 = 0.2813
bo7 = -0.1601
bo8 = 0.1243
bo9 = 0.1171
bo10 = -0.143
bo11 = 0.0641
bo12 = -0.5204
bo13 = -0.1603
bo14 = -0.184
bo15 = 0.0909
bo16 = 0.5073
bo17 = -0.4517
bo18 = -0.7359
bo19 = -0.8512
bo20 = -0.6642
bo21 = -0.0853  #vl 400-<3599
bo22 = 0.1042   #vl 3500-<10K
bo23 = 0.1807  #10k-<5
bo24 = -0.1015  #>=5
bo25 = 0.2626  #cd4 <50
bo26 = 0.5522  #cd4 50-<100
bo27 = 0.1098  #cd4 100-<200
bo28 = 0.1602  #cd4 200-<350
bo29 = 0.0636  #cd4 >=350
bo30 = 0.2924  #bmi <20
bo31 = 0.308  #bmi 25-<30
bo32 = 0.5335  #bmi >=30
bo33 = 0.0387    #bpd
bo34 = -0.0225 #bps
bo35 = 0.00996   #hdl
bo36 = -0.00186 #ldl
bo37 = -0.00056 #trig



