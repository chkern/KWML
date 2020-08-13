## code to prepare `simu_dat` dataset goes here

library(survey)
library(cobalt)

# Load R functions for pseudo weights calculation
source("weighting_functions.R")
source("subfunctions.R")
seeds = read.table("seed.txt", header = T)
seed1 = seeds[, 1]
seed2 = seeds[, 2]

# Population generation
N =100000
set.seed(9762342)
v_w = as.data.frame(matrix(rnorm(N*10), N, 10))
names(v_w) = c(paste0(rep("v", 8), c(1:6, 8, 9), sep = ""), "w7", "w10")
v_w$w5 = (v_w$v1*0.16+v_w$v5*0.84)*1.171
v_w$w6 = (v_w$v2*0.67+v_w$v6*0.33)*1.353
v_w$w8 = (v_w$v3*0.18+v_w$v8*0.82)*1.177
v_w$w9 = (v_w$v4*0.67+v_w$v9*0.33)*1.317
pop=v_w[,c(1:4, 11:14, 9, 10)]
pop = pop[,order(as.numeric(substr(names(pop), 2, 3)))]
names(pop) = paste0(rep("w", 10), c(1:10), sep = "")
#round(cor(pop), 1)
beta  = c(0, 1, 1, 1.5, 1.5, -.8, -.5, .7)
n.beta=length(beta)
alpha = c(-2.5, 1, 1, 1, 1, .71, -.19, .26)
n.alpha = length(alpha)
odds_y = exp(as.matrix(cbind(1, pop[,c(1:4, 8:10)]))%*%matrix(alpha, n.alpha, 1))
py = odds_y/(1+odds_y)
pop$y=as.numeric((runif(N)<py)); mean(pop$y)
pop$w2_2 = pop$w2^2
pop$w4_2 = pop$w4^2
pop$w7_2 = pop$w7^2
pop$w1_w3 = .5*pop$w1*pop$w3
pop$w2_w4 = .7*pop$w2*pop$w4
pop$w4_w5 = .5*pop$w4*pop$w5
pop$w5_w6 = .5*pop$w5*pop$w6
pop$w3_w5 = .5*pop$w3*pop$w5
pop$w4_w6 = .7*pop$w4*pop$w6
pop$w5_w7 = .5*pop$w5*pop$w7
pop$w1_w6 = .5*pop$w1*pop$w6
pop$w2_w3 = .7*pop$w2*pop$w3
pop$w3_w4 = .5*pop$w3*pop$w4
pop$w3_2_w5_2 = 0.6*pop$w3^2*pop$w5^2
pop$w1_w2_w3  = 0.6*pop$w1*pop$w2*pop$w3
pop$w4_w5_w7  = 0.6*pop$w4*pop$w5*pop$w7

pop$w1_c = cut(pop$w1, c(-5, -1, 0, 1, 6))
pop$w2_c = cut(pop$w2, c(-5, -1, 0, 1, 6))
pop$w3_c = cut(pop$w3, c(-5, -1, 0, 1, 6))
pop$w4_c = cut(pop$w4, c(-5, -1, 0, 1, 6))

pop$w1_s = 2*sqrt(abs(pop$w1+pop$w2+pop$w3))
pop$w2_s = 2*sqrt(abs(pop$w2+pop$w3+pop$w4))
pop$w3_s = 2*sqrt(abs(pop$w3+pop$w4+pop$w5))
pop$w4_s = 2*sqrt(abs(pop$w4+pop$w5+pop$w7))
#round(cor(pop[,c(1:24,32:35)]),2)

pop$w1_n = pop$w1+rnorm(N)/10+pop$w8/10
pop$w2_n = pop$w2+rnorm(N)+pop$w9/10
pop$w3_n = pop$w3+rnorm(N)+pop$w10/10
pop$w4_n = pop$w4+rnorm(N)+pop$y/2
#round(cor(pop[,c(1:24,36:39)]),2)


ds.m_9 = model.matrix(as.formula(paste("y~",paste(names(pop)[28:31],collapse="+"))), data = pop)
beta_9 = c(0, -0.8, -0.5, 0.5, 0.2, 2, 0.5, 3, 1, 0.5, -1, 2, 0.5)

covars   = names(pop)[1:7]
c_covars = names(pop)[1:7]
f_covars = NULL
if(sum(!(covars%in%c_covars))) f_covars=covars[!(covars%in%c_covars)]

# odds=matrix(0, N, 10)
# odds[,1]  = exp(as.matrix(cbind(1, pop[,c(1:7)]))%*%matrix(beta, n.beta, 1))
# odds[,2]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12)]))%*%matrix(beta[c(1:8, 3)], n.beta+1, 1))
# odds[,3]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14)]))%*%matrix(beta[c(1:8, 3, 5, 8)], n.beta+3, 1))
# odds[,4]  = exp(as.matrix(cbind(1, pop[,c(1:7, 15:18)]))%*%matrix(beta[c(1:8, 2, 3, 5, 6)], n.beta+4, 1))
# odds[,5]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12, 15:18)]))%*%matrix(beta[c(1:8, 3, 2, 3, 5, 6)], n.beta+5, 1))
# odds[,6]  = exp(as.matrix(cbind(1, pop[,c(1:7, 15, 16, 19:24, 17, 18)]))%*%matrix(beta[c(1:8, 2:6, 2:6)], n.beta+10, 1))
# odds[,7]  = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14, 15, 16, 19:24, 17, 18)]))%*%matrix(beta[c(1:8, 3, 5, 8, 2:6, 2:6)], n.beta+13, 1))
odds  = exp(as.matrix(cbind(1, pop[,c(1:7, 12:14, 15:18, 25:27)]))%*%matrix(beta[c(1:8, 3, 5, 8, 2, 3, 5, 6, 4, 2, 3)], n.beta+10, 1))

q_c = odds^0.18
q_s = odds^-0.11

#q = as.data.frame(odds/(1+odds))
#q_c = cbind(odds[,1]^0.3,   odds[,2]^0.25,  odds[,3]^0.25,  odds[,4]^0.27,  odds[,5]^0.25,
#            odds[,6]^0.22,  odds[,7]^0.17,  odds[,8]^0.18,  odds[,9]^0.4,   odds[,10]^0.23)
#q_s = cbind(odds[,1]^-0.2,  odds[,2]^-0.2,  odds[,3]^-0.2,  odds[,4]^-0.17, odds[,5]^-0.17,
#            odds[,6]^-0.15, odds[,7]^-0.13, odds[,8]^-0.11, odds[,9]^-0.26, odds[,10]^-0.15)

n_c = 2000
n_s = 2000

samp.c = samp.slct(seed = seed1[simu],
                   fnt.pop = pop,
                   n = n_c,
                   dsgn = "pps",
                   size = q_c)

samp.s = samp.slct(seed = seed1[simu],
                   fnt.pop = pop,
                   n = n_s,
                   dsgn = "pps",
                   size = q_s)

# Combine survey and cohort data
psa_dat = rbind(samp.c, samp.s)
psa_dat$elig_wt = c(rep(1, n_c), samp.s$wt)
psa_dat$trt_n = c(rep(1, n_c), rep(0, n_s))
psa_dat$trt   = as.factor(psa_dat$trt_n)
# Name of data source indicator in the combined sample
rsp_name = "trt" # 1 for AARP, 0 for NHIS
