
##############################################################################
# x_bar_var is a function calculating weighted sample mean and variance      #
# for the covariates (results are consistent with svymean in survey package) #
# INPUT                                                                      #
#  dat: a dataframe including all covariates                                 #
#  wt:  sample weights                                                       #
#  covars:  a vector of covariate names                                      #
# OUPUT                                                                      #
#  X_bar:    sample estimate of population mean                              #
#  X_sigma2: sample estimate of population variance                          #
##############################################################################

x_bar_var = function(dat, wt, covars){
  # Split covariates
  catvars = names(dat[, covars])[sapply(dat[, covars], is.factor)]
  numvars = names(dat[, covars])[!sapply(dat[, covars], is.factor)]

  N_hat = sum(wt)
  n = nrow(dat)
  Xbar_hat_co = NULL
  Xbar_hat_cat = NULL
  var_Xhat_co = NULL
  var_Xhat_cat = NULL

  if (!identical(numvars, character(0))) {
    # Create data frame for numeric covars
    x_matrix_co = dat[, numvars, drop = F]
    # Means and vars for numeric covars
    p = ncol(x_matrix_co)
    wt_matrix = matrix(rep(wt, p), n, p)
    X_hat_co = apply(x_matrix_co*wt_matrix, 2, sum)
    Xbar_hat_co = X_hat_co/N_hat
    var_Xhat_co = c(apply(x_matrix_co, 2, var))
  }

  if (!identical(catvars, character(0))) {
    # Categorize categorical covars and create data frame
    formula = as.formula(paste("trt~", paste0(catvars, collapse = "+")))
    x_matrix_cat = as.data.frame(model.matrix(formula,
                                              data = dat,
                                              contrasts.arg = lapply(dat[sapply(dat, is.factor)],
                                                                     contrasts, contrasts = F))[,-1])
    # Means and vars for categorical covars
    p = ncol(x_matrix_cat)
    wt_matrix = matrix(rep(wt, p), n, p)
    X_hat_cat = apply(x_matrix_cat*wt_matrix, 2, sum)
    Xbar_hat_cat = X_hat_cat/N_hat
    Xbar_cat = c(apply(x_matrix_cat, 2, mean))
    var_Xhat_cat = c(Xbar_cat*(1-Xbar_cat))
  }

  Xbar_hat = c(Xbar_hat_co, Xbar_hat_cat)
  var_Xhat = c(var_Xhat_co, var_Xhat_cat)

  return(list(X_bar = Xbar_hat, X_sigma2 = var_Xhat))
}


#############################################################################################
# smd is a function to calculate SMD for two samples                                        #
# INPUT                                                                                     #
#  psa_dat: a dataframe including covariates to be compared in the two samples              #
#           NOTE: categorical variables should have a format of factor!                     #
#  wt:      a vector of weights for nonprobability sample and the survey sample.            #
#           The order should be consistent with psa_dat                                     #
#  rsp_name: a vector of an indicator for sample membership in the dataframe of psa_dat     #
#           1 for nonprobabilty sample, 0 for survey sample                                 #
#  covars:  a vector of covariate names                                                     #
# OUTPUT                                                                                    #
#  smd:     a vector of smd for all compared covariates                                     #
#############################################################################################

smd = function(psa_dat, wt, rsp_name, covars){
  # separate datasets for nonprobability and survey samples
  np_dat  = psa_dat[psa_dat[, rsp_name]==1,]
  svy_dat = psa_dat[psa_dat[, rsp_name]==0,]
  np_dat$trt = 1
  svy_dat$trt = 0
  # vectors of the weights for the two samples
  pwt   = psa_dat[psa_dat[, rsp_name]==1, wt]
  svywt = psa_dat[psa_dat[, rsp_name]==0, wt]
  # Get pseudo-weighted nonprobability sample and sample weighted survey estimates of pop means & vars
  np_x_est  = x_bar_var(dat = np_dat, wt = pwt, covars = covars)
  svy_x_est = x_bar_var(dat = svy_dat, wt = svywt, covars = covars)
  # Calculate smd for each covariate
  smd = (np_x_est$X_bar-svy_x_est$X_bar)/sqrt((np_x_est$X_sigma2+ svy_x_est$X_sigma2)/2)
  smd
}


################################################################################
# FUNCTION pop.gen is to generate a finite population                          #
# INPUT:                                                                       #
#     seed:       random seed that enables replication                         #
#     N:          population size                                              #
#     gamma:      regression coefficients for generating outcome y             #
#     beta_star:  scaling factor for two-way interaction terms of covariates   #
#     beta_sstar: scaling factor for more complex interactions of covariates   #
# OUTPUT:                                                                      #
#     pop: a data frame pop including                                          #
#       - covariates X1... X10, selected quadratic terms, interactions, and    #
#         higher order terms;                                                  #
#       - covariates X1*...X4*, and X1**...X4** generated from X1...X4 ;       #
#       - outcome variable y;                                                  #
################################################################################

pop.gen = function(seed = 9762342,
                   N = 10000,
                   gamma = c(-2.5, 1, 1, 1, 1, .71, -.19, .26),
                   beta_star = c(0.5, 0.7, 0.5, 0.5, 0.5, 0.7, 0.7),
                   beta_sstar = rep(0.6, 7)){
    # Create the base covariates (V1...V10)
	v_w = as.data.frame(matrix(rnorm(N*10), N, 10))
	# Create covariates X1...X10
	names(v_w) = c(paste0(rep("v", 8), c(1:6, 8, 9), sep = ""), "w7", "w10")
	v_w$w5 = (v_w$v1*0.16+v_w$v5*0.84)*1.171
	v_w$w6 = (v_w$v2*0.67+v_w$v6*0.33)*1.353
	v_w$w8 = (v_w$v3*0.18+v_w$v8*0.82)*1.177
	v_w$w9 = (v_w$v4*0.67+v_w$v9*0.33)*1.317
	pop=v_w[,c(1:4, 11:14, 9, 10)]
	pop = pop[,order(as.numeric(substr(names(pop), 2, 3)))]
	names(pop) = paste0(rep("w", 10), c(1:10), sep = "")
	# Generate outcome variable y
	n.gamma = length(gamma)
	odds_y = exp(as.matrix(cbind(1, pop[,c(1:4, 8:10)]))%*%matrix(gamma, n.gamma, 1))
	py = odds_y/(1+odds_y)
	pop$y=as.numeric((runif(N)<py)); mean(pop$y)
	# Generate non-linear and non-additive terms of X1...X10
	pop$w2_2 = pop$w2^2
	pop$w4_2 = pop$w4^2
	pop$w7_2 = pop$w7^2
	pop$w1_w3 = beta_star[1]*pop$w1*pop$w3
	pop$w2_w4 = beta_star[2]*pop$w2*pop$w4
	pop$w4_w5 = beta_star[4]*pop$w4*pop$w5
	pop$w5_w6 = beta_star[5]*pop$w5*pop$w6
	pop$w3_w5 = beta_star[3]*pop$w3*pop$w5
	pop$w4_w6 = beta_star[6]*pop$w4*pop$w6
	pop$w5_w7 = beta_star[5]*pop$w5*pop$w7
	pop$w1_w6 = beta_star[1]*pop$w1*pop$w6
	pop$w2_w3 = beta_star[2]*pop$w2*pop$w3
	pop$w3_w4 = beta_star[3]*pop$w3*pop$w4
	pop$w3_2_w5_2 = beta_sstar[3]*pop$w3^2*pop$w5^2
	pop$w1_w2_w3  = beta_sstar[1]*pop$w1*pop$w2*pop$w3
	pop$w4_w5_w7  = beta_sstar[4]*pop$w4*pop$w5*pop$w7

	# Generate X1*...X4* (categorized X1...X4)
	pop$w1_c = cut(pop$w1, c(-5, -1, 0, 1, 6))
	pop$w2_c = cut(pop$w2, c(-5, -1, 0, 1, 6))
	pop$w3_c = cut(pop$w3, c(-5, -1, 0, 1, 6))
	pop$w4_c = cut(pop$w4, c(-5, -1, 0, 1, 6))
	# Generate X1**...X4**
	pop$w1_n = pop$w1+rnorm(N)/10+pop$w8/10
	pop$w2_n = pop$w2+rnorm(N)+pop$w9/10
	pop$w3_n = pop$w3+rnorm(N)+pop$w10/10
	pop$w4_n = pop$w4+rnorm(N)+pop$y/2
	pop
} # End of pop.gen

##########################################################################################
# FUNCTION mos.gen is to generate measure of size (MOS) of                               #
#                  Probability Proportional to Size (PPS) sampling for cohort and survey #
# INPUT:                                                                                 #
#     pop:   dataframe of population including all covariates and outcome variable       #
#     beta:  Coefficients of X1...X7 for calculating the MOS                             #
#     alpha: Coefficients of X1*...X4* for calculating the MOS                           #
# OUTPUT:                                                                                #
#     q_c: MOS for cohort selection                                                      #
#     q_s: MOS for survey sample selection                                               #
##########################################################################################

mos.gen = function(pop, beta, alpha){
	N = dim(pop)[1]
	n.beta = length(beta)
	# convert X1*...X4* to a matrix of dummy variables
	ds.m_9 = model.matrix(as.formula(paste("y~",paste(names(pop)[28:31],collapse="+"))), data = pop)
	odds=matrix(0, N, 10)
    odds[,1]  = exp(as.matrix(cbind(1, pop[,c(1:7)]))%*%matrix(beta, n.beta, 1))
    odds[,2]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12)]))%*%matrix(beta[c(1:8, 3)], n.beta+1, 1))
    odds[,3]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14)]))%*%matrix(beta[c(1:8, 3, 5, 8)], n.beta+3, 1))
    odds[,4]  = exp(as.matrix(cbind(1, pop[,c(1:7, 15:18)]))%*%matrix(beta[c(1:8, 2, 3, 5, 6)], n.beta+4, 1))
    odds[,5]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12, 15:18)]))%*%matrix(beta[c(1:8, 3, 2, 3, 5, 6)], n.beta+5, 1))
    odds[,6]  = exp(as.matrix(cbind(1, pop[,c(1:7, 15, 16, 19:24, 17, 18)]))%*%matrix(beta[c(1:8, 2:6, 2:6)], n.beta+10, 1))
    odds[,7]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14, 15, 16, 19:24, 17, 18)]))%*%matrix(beta[c(1:8, 3, 5, 8, 2:6, 2:6)], n.beta+13, 1))
    odds[,8]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14, 15:18, 25:27)]))%*%matrix(beta[c(1:8, 3, 5, 8, 2, 3, 5, 6, 4, 2, 3)], n.beta+10, 1))
    odds[,9] = exp(ds.m_9%*%matrix(alpha, length(alpha), 1))
    odds[,10] = exp(as.matrix(cbind(1, pop[,c(32:35)]))%*%matrix(beta[1:5], 5, 1))
    q_c = cbind(odds[,1]^0.3,   odds[,2]^0.25,  odds[,3]^0.25,  odds[,4]^0.27,  odds[,5]^0.25,
                odds[,6]^0.22,  odds[,7]^0.17,  odds[,8]^0.18,  odds[,9]^0.4,   odds[,10]^0.23)
    q_s = cbind(odds[,1]^-0.2,  odds[,2]^-0.2,  odds[,3]^-0.2,  odds[,4]^-0.17, odds[,5]^-0.17,
                odds[,6]^-0.15, odds[,7]^-0.13, odds[,8]^-0.11, odds[,9]^-0.26, odds[,10]^-0.15)
    return(list(q_c = q_c,
                q_s = q_s))

}# End of mos.gen


#####################################################################################################
# FUNCTION samp.slct to select a sample                                                             #
# The sampling design is described in simulation section                                            #
# INPUT:                                                                                            #
#     seed:     random seed that enables replication                                                #
#     fnt.pop:  the population including response and covariates                                    #
#               including, psu, gender, age, urban/rural, hh income, race/eth,                      #
#               Env, and y                                                                          #
#     n:        sample size                                                                         #
#     Cluster:  total number of clusters in the finite population                                   #
#     Clt.samp: number of clusters to be selected in the sample                                     #
#     dsgn:     sampling designs                                                                    #
#                                                                                                   #
# OUTPUT:                                                                                           #
#     returns a dataframe of selected samp, including variables of                                  #
#      (ID = ID, psu = psu, gender = gender, age = age, urb= urb,                                   #
#       hh_inc = hh_inc, race_eth = race_eth, y = fnt.y, weights)                                   #
#####################################################################################################

samp.slct = function(seed, fnt.pop, n, Cluster=NULL, Clt.samp=NULL, dsgn, size = NULL, size.I = NULL){
	set.seed(seed)
	N = nrow(fnt.pop)
	size.Cluster = N/Cluster
	# one-ste sample design
	if(dsgn=="pps"){
	  fnt.pop$x=size
	  samp = sam.pps(fnt.pop,size, n)
	}
	# two-stage cluster sampling design (informative design at the second stage)
	if(dsgn == "srs-pps"){

	  # calcualte the size for the second stage
	  fnt.pop$x = size

	  #-- first stage: select clusters by srs
	  # Clt.samp = 25
	  # n = 1000
	  index.psuI = sample(1:Cluster, Clt.samp, replac = F)
	  index.psuI = sort(index.psuI)  #sort selected psus
	  sample.I = fnt.pop[fnt.pop$psu %in% index.psuI,]
	  sample.I$wt.I = Cluster/Clt.samp

	  #-- second stage: select subj.samp within selected psus by pps to p^a (or (1-p)^a)
	  samp=NULL
	  for (i in 1: Clt.samp){
	    popn.psu.i= sample.I[sample.I$psu==index.psuI[i, 1],]
	    size.II.i = sample.I[sample.I$psu==index.psuI[i, 1],"x"]
	    samp.i = sam.pps(popn.psu.i,size.II.i, n/Clt.samp)
	    samp.i$wt = samp.i$wt*samp.i$wt.I
	    samp = rbind(samp,samp.i)
	  }#sum(samp.cntl$wt);nrow(fnt.cntl)
	}
	if(dsgn == "pps-pps"){

	  # calcualte the size for the second stage
	  fnt.pop$x = size

	  #-- first stage: select clusters by pps to size aggrated level p^a(or (1-p)^a)
	  # Clt.samp = 25
	  # n = 1000
	  index.psuI = sam.pps(matrix(1:Cluster,,1),size.I, Clt.samp)
	  index.psuI = index.psuI[order(index.psuI[,1]),]  #sort selected psus
	  sample.I = fnt.pop[fnt.pop$psu %in% index.psuI[,1],]
	  sample.I$wt.I = rep(index.psuI[,'wt'],each=size.Cluster)

	  #-- second stage: select subj.samp within selected psus by pps to p^a (or (1-p)^a)
	  samp=NULL
	  for (i in 1: Clt.samp){
	    popn.psu.i= sample.I[sample.I$psu==index.psuI[i,1],]
	    size.II.i = sample.I[sample.I$psu==index.psuI[i,1],"x"]
	    samp.i = sam.pps(popn.psu.i,size.II.i, n/Clt.samp)
	    samp.i$wt = samp.i$wt*samp.i$wt.I
	    samp = rbind(samp,samp.i)
	  }#sum(samp.cntl$wt);nrow(fnt.cntl)
	}
	rownames(samp) = as.character(1:dim(samp)[1])
  return(samp)
}


#################################################################################
#-FUNCTION sam.pps is to select a sample with pps sampling                      #
# INPUT:  popul - the population including response and covariates              #
#         MSize - the size of measure to compute the selection probabilities(>0)#
#         n -     the sample size                                               #
# OUTPUT: sample selected with pps sampling of size n,                          #
#         including response, covariates, and sample weights                    #
#################################################################################
sam.pps<-function(popul,Msize, n){
  N=nrow(popul)
  pps.samID=sample(N,n,replace=F,prob=Msize)   #F but make sure N is large!
  if (dim(popul)[2] == 1){
  	sam.pps.data=as.data.frame(popul[pps.samID,])
  	names(sam.pps.data) = names(popul)
  }else{sam.pps.data=popul[pps.samID,]}
  sam.pps.data$wt=sum(Msize)/n/Msize[pps.samID]
  return(sam.pps = sam.pps.data)
}


#####################################################################################################################
## cmb_dat is a function for data preparation including removing missing values in cohort and survey sample,       ##
## combining the two samples                                                                                       ##
## INPUT:  chtsamp - cohort (a data frame including covariates of the propensity score model)                      ##
##         svysamp - survey sample (a data frame including weight variable and                                     ##
##                                  covariates of the propensity score model)                                      ##
##         svy_wt  - variable name of the survey weight (should be included in svysamp)                            ##
##         Formula - propensity score model, with the format of trt~covariate1+covariate2+...                      ##
## OUTPUT: cmb_dat - combined sample of non-missing cohort and survey sample units including varibles              ##
##                   covariates in propensity score model                                                          ##
##                   "trt" indicator of data source (1 for cohort units, 0 for survey sample units)                ##
##                   "wt" weight varaible (1 for cohort units, survey sample weights)                              ##
##         chtsamp - complete  sample (in terms of covariates)                                                     ##
##         svysamp - complete survey sample (in terms of covariates)                                               ##
## WARNINGS:                                                                                                       ##
##         1, missing values in covariates are not allowed. Records with missing values in the cohort are removed. ##
##         2, missing values in covariates are not allowed. Records with missing values in the survey sample       ##
##         are removed. The complete cases are reweighted. Missing completely at random is assumed.                ##
#####################################################################################################################

cmb_dat = function(chtsamp,svysamp,svy_wt, Formula){
  # Get names of the response variable and predictors for the propensity score estimation model
  Fml_names = all.vars(Formula)
  # Name of the response variable
  rsp_name = Fml_names[1]
  # Name of the predictors
  mtch_var = Fml_names[-1]

  # Remove incomplete records in the cohort, if there are any.
  chtsamp_sub = as.data.frame(chtsamp[, mtch_var])
  if(sum(is.na(chtsamp_sub))>0){
    cmplt.indx = complete.cases(chtsamp_sub)
    chtsamp_sub = chtsamp_sub[cmplt.indx, ]
    chtsamp = chtsamp[cmplt.indx,]
    warning("Missing values in covariates are not allowed. Records with missing values in the cohort are removed.")
  }
  # Remove incomplete survey sample, if there are any.
  svysamp_sub = as.data.frame(svysamp[, mtch_var])
  svy_wt.vec = c(svysamp[, svy_wt])
  if(sum(is.na(svysamp_sub))>0){
    cmplt.indx = complete.cases(svysamp_sub)
    svysamp_sub = svysamp_sub[cmplt.indx, ]
    svy_wt.vec = sum(svysamp[, svy_wt])/sum(svy_wt.vec[cmplt.indx])*svy_wt.vec[cmplt.indx]
    warning("Missing values in covariates are not allowed. Records with missing values in the survey sample are removed.
            The complete cases are reweighted. Missing completely at random is assumed.")
  }# end removing incomplete records
  # size of cohort (complete)
  m = dim(chtsamp_sub)[1]
  # size of survey sample (complete)
  n = dim(svysamp_sub)[1]
  # Combine the two complete samples
  # Set outcome variable for propensity score model
  chtsamp_sub[,rsp_name] = 1                 # z=1 for cohort sample
  svysamp_sub[,rsp_name] = 0                 # z=0 for survey sample
  names(chtsamp_sub) = c(mtch_var, rsp_name) # unify the variable names in cohort and survey sample
  names(svysamp_sub) = c(mtch_var, rsp_name)
  cmb_dat = rbind(chtsamp_sub, svysamp_sub)  # combine the two samples
  cmb_dat$wt = c(rep(1, m), svysamp[, svy_wt])
  return(list(cmb_dat = cmb_dat,
              chtsamp = chtsamp,
              svysamp = svysamp
  )
  )
  }


#####################################################################################################################
## cmb_dat1 is a function combining the cohort and survey sample                                                   ##
## INPUT:  chtsamp - cohort (a data frame including covariates of the propensity score model)                      ##
##         svysamp - survey sample (a data frame including weight variable and                                     ##
##                                  covariates of the propensity score model)                                      ##
##         svy_wt  - variable name of the survey weight (should be included in svysamp)                            ##
##         Formula - propensity score model, with the format of trt~covariate1+covariate2+...                      ##
## OUTPUT: cmb_dat - combined sample of non-missing cohort and survey sample units including varibles              ##
##                   covariates in propensity score model                                                          ##
##                   "trt" indicator of data source (1 for cohort units, 0 for survey sample units)                ##
##                   "wt" weight varaible (1 for cohort units, survey sample weights)                              ##
#####################################################################################################################

cmb_dat1 = function(chtsamp,svysamp,svy_wt, Formula){
  # Get names of the response variable and predictors for the propensity score estimation model
  Fml_names = all.vars(Formula)
  # Name of the response variable
  rsp_name = Fml_names[1]
  # Name of the predictors
  mtch_var = Fml_names[-1]

  chtsamp_sub = as.data.frame(chtsamp[, mtch_var])
  svysamp_sub = as.data.frame(svysamp[, mtch_var])
  svy_wt.vec = c(svysamp[, svy_wt])
  m = dim(chtsamp_sub)[1]
  # size of survey sample (complete)
  n = dim(svysamp_sub)[1]
  # Combine the two complete samples
  # Set outcome variable for propensity score model
  chtsamp_sub[,rsp_name] = 1                 # z=1 for cohort sample
  svysamp_sub[,rsp_name] = 0                 # z=0 for survey sample
  names(chtsamp_sub) = c(mtch_var, rsp_name) # unify the variable names in cohort and survey sample
  names(svysamp_sub) = c(mtch_var, rsp_name)
  cmb_dat = rbind(chtsamp_sub, svysamp_sub)  # combine the two samples
  cmb_dat$wt = c(rep(1, m), svy_wt.vec)
  cmb_dat
}

