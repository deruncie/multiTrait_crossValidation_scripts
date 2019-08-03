# full simulation script
# using Wheat K
# tests:
# 1. known G,R:
#    a. simulate data, hold out 10%, calculate cor(u,u_hat), cor(y,u_hat) analytically, for single, CV1, CV2
#    b. simulate data, hold out 10%, simulate cor(u,u_hat), cor(y,u_hat) for single, CV1, CV2
# 2. baseline
#    a. simulate data, hold out 10%, estimate cor(u,u_hat), cor(y,u_hat) for single, CV1, CV2
# 3. parametric
#    a. simulate data, hold out 50%, estimate cor_g(y,y_hat)/sqrt(h2_y_hat)
# 4. semi-parametric
#    a. simulate data, hold out 10%, calculate cor(u,u_hat), cor(y,u_hat) for single, CV1, CV2 with correction
# 5. non-parametric
#    a. simulate data, hold out 10%, split 2 clones, estimate cor(u,u_hat), cor(y,u_hat) for single, CV1, CV2, each with 2nd clone. Adjust residual variance
#    b. simulate data, hold out 20% (10% validation + nearest relatives), estimate cor(u,u_hat), cor(y,u_hat) for single, CV1, CV2, each with relative

library(lme4)
library(lme4qtl)
library(Matrix)
library(rrBLUP)
library(foreach)
library(doParallel)
library(MASS)
library(BGLR)
source('my_relmatLmer.R')

# set up run
runID = as.numeric(commandArgs(t=T)[1])
if(is.na(runID)) {
  runID = 1
} else{
  set.seed(runID)
}
registerDoParallel(1)
nReps = 1

# get K matrix
load('FileS3.rdata')
K = G
rm(list=c('G','Y'))

n = nrow(K)
rownames(K) = colnames(K) = 1:nrow(K)

#------------------------------------------------------------#
# Test 2.a., 4.a.
#------------------------------------------------------------#
try(dir.create('Results'))
try(dir.create('Results/Raw'))
outfile = sprintf('Results/Raw/test_2a4a_run_%03d.csv',runID)
if(file.exists(outfile)) quit('no')

# Select training, validation sets
n_validation = round(.1*n)
validation = sample(1:n,n_validation)
training = 1:n
training = training[training %in% validation == F]
n_validation = length(validation)
n_train = length(training)

# create data.frames
wide_data_full = data.frame(Line = factor(1:n),ID=factor(1:n))
wide_data_validation = wide_data_full[validation,,drop=F]
wide_data_train = wide_data_full[training,,drop=F]
tall_data_full = rbind(data.frame(wide_data_full,Trait=1),data.frame(wide_data_full,Trait=2))
tall_data_full$Trait = factor(tall_data_full$Trait)
tall_data_validation = tall_data_full[c(validation,n+validation),]
tall_data_train = tall_data_full[c(training,n+training),]

#------------------------------------------------------------#
# Set up data.frames, matrices
#------------------------------------------------------------#
{
  # get rotation matrix, pre-rotate K and ~Trait
  chol_K = chol(K)
  svd_K = svd(K)
  Q = svd_K$u
  D_K = Diagonal(n,svd_K$d)
  rownames(D_K) = colnames(D_K) = wide_data_full$Line
  X = kronecker(diag(1,2),t(Q)) %*% model.matrix(~0+Trait,tall_data_full)

  K_nn = K[validation,validation]
  svd_K_nn = svd(K_nn)
  Q_nn = svd_K_nn$u
  D_K_nn = Diagonal(n_validation,svd_K_nn$d)
  rownames(D_K_nn) = colnames(D_K_nn) = wide_data_validation$Line
  X_validation = kronecker(diag(1,2),t(Q_nn)) %*% model.matrix(~0+Trait,tall_data_validation)

  K_oo = K[training,training]
  svd_K_oo = svd(K_oo)
  Q_oo = svd_K_oo$u
  D_K_oo = Diagonal(n_train,svd_K_oo$d)
  rownames(D_K_oo) = colnames(D_K_oo) = wide_data_train$Line
  X_train = kronecker(diag(1,2),t(Q_oo)) %*% model.matrix(~0+Trait,tall_data_train)

  K_no = K[validation,training]
  K_c = K_nn - K_no %*% solve(K_oo,t(K_no)) # (K^{-1})^{-1}_nn
  K_oo_inv = Q_oo %*% diag(1/diag(D_K_oo)) %*% t(Q_oo)
  K_inv = Q %*% diag(1/diag(D_K)) %*% t(Q)

  # other matrices
  I = Diagonal(n,1)
  I_oo = Diagonal(n_train,1)
  I_nn = Diagonal(n_validation,1)
  I = diag(1,n)
  I_oo = diag(1,n_train)
  I_nn = diag(1,n_validation)
  S = diag(1,n_validation) - 1/n_validation
}

#------------------------------------------------------------#
# Loop through parameter sets
#------------------------------------------------------------#

# set base parameter values
cor_G = 0.8
cor_R = -0.8
cor_R2 = 0
H21 = 0.2
H22 = 0.2
Rep=1


# generate random numbers
U_randn = matrix(rnorm(n*2),n,2)
E_randn = matrix(rnorm(n*2),n,2)


# registerDoParallel(detectCores())
# loop through parameter values
results = foreach(cor_G = c(0,0.3,0.6),.combine = 'rbind',.errorhandling = 'remove') %do% { #0.2,0.4,
  foreach(cor_R = c(-0.9,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.9),.combine = 'rbind',.errorhandling = 'remove') %dopar% {
    foreach(H21 = c(0.2,.6),.combine = 'rbind',.errorhandling = 'remove') %do% { #0.1,0.4,0.7
      foreach(H22 = c(0.2,0.6),.combine = 'rbind',.errorhandling = 'remove') %do% { #0.4,

        # set up genetic parameters
        {
          H2s = c(H21,H22)
          Gcor = diag(1,2);
          Gcor[Gcor==0] = cor_G
          G = diag(sqrt(H2s)) %*% Gcor %*% diag(sqrt(H2s))
          Gsqrt = chol(G)
          Rcor = diag(1,2);
          Rcor[Rcor==0] = cor_R2
          Rcor[1,-1] = Rcor[-1,1] = cor_R
          R = diag(sqrt(1-H2s)) %*% Rcor %*% diag(sqrt(1-H2s))
          Rsqrt = chol(R)
        }

        # create traits
        print(c(cor_G,cor_R,H21,H22))
        # set.seed(Rep)
        # U = t(chol_K) %*% U_randn %*% Gsqrt
        U = Q %*% (sqrt(diag(D_K))*U_randn) %*% Gsqrt
        E = E_randn %*% Rsqrt     # environmental error
        Y = as.matrix(U + E)

        Y_o = Y[training,]
        U_n = U[validation,]
        Y_n = Y[validation,]

        # fit single-trait model to training data
        # m1b = relmatLmer(Y_o[,1]~(1|Line),wide_data_full,relmat = list(Line = K))
        m1 = relmatLmer(t(Q_oo) %*% Y_o[,1]~0 + (t(Q_oo) %*% rep(1,n_train)) + (1|Line),wide_data_train,relmat = list(Line = D_K_oo))
        var1 = as.data.frame(VarCorr(m1))
        g11 = var1$vcov[1]
        r11 = var1$vcov[2]

        # single trait predictions
        V_o1 = g11*K_oo + r11*I_oo
        uhat_single = as.matrix(g11*K_no %*% solve(V_o1,Y_o[,1] - fixef(m1)))
        cor(U[validation,1],uhat_single)

        # fit joint model to training data
        m2 = relmatLmer(c(t(Q_oo) %*% Y_o)~0+X_train+(0+Trait|Line)+(0+Trait|ID),tall_data_train,relmat = list(Line = D_K_oo))
        var2 = as.data.frame(VarCorr(m2))
        Rhat = diag(var2$vcov[7],2)
        Rhat[1,1] = Rhat[1,1] + var2$vcov[4]
        Rhat[2,2] = Rhat[2,2] + var2$vcov[5]
        Rhat[1,2]=Rhat[2,1] = var2$vcov[6]
        Ghat = diag(1,2)
        Ghat[1,1] = var2$vcov[1]
        Ghat[2,2] = var2$vcov[2]
        Ghat[1,2]=Ghat[2,1] = var2$vcov[3]

        # CV1-style predictions
        V_o = kronecker(Ghat,K_oo) + kronecker(Rhat,I_oo)
        uhat_o = kronecker(Ghat,K_oo) %*% solve(V_o,c(sweep(Y_o,2,fixef(m2),'-')))
        uhat_o = matrix(uhat_o,nc=2)
        uhat_CV1 = K_no %*% K_oo_inv %*% uhat_o

        # CV2-style predictions
        V_c = Ghat[2,2]*K_c + Rhat[2,2]*I_nn
        uhat_CV2 = as.matrix(uhat_CV1[,1] + (Ghat[1,2]*K_c) %*% solve(V_c,Y_n[,2]-fixef(m2)[2] - uhat_CV1[,2]))

        # fit joint model to full data
        m2b = relmatLmer(c(t(Q) %*% Y)~0+X+(0+Trait|Line)+(0+Trait|ID),tall_data_full,relmat = list(Line = D_K))
        # predictions from full data - to use as candidate validation data
        U_n_uhatFull = Q %*% as.matrix(t(m2b@optinfo$relmat$relfac$Line) %*% as.matrix(ranef(m2b)$'Line'))
        U_n_uhatFull = U_n_uhatFull[validation,]

        # calculate empirical correction factor
        cov_correction = sum(diag(S %*% (Ghat[2,1]*K_c %*% solve(V_c,Rhat[2,1]*I_nn)))) / (n_validation-1)

        # calculate correlations
        results = data.frame(
          Single_u = cor(uhat_single,U_n[,1]),
          Single_y = cor(uhat_single,Y_n[,1]),
          Single_uhatFull = cor(uhat_single,U_n_uhatFull[,1]),
          CV1_u = cor(uhat_CV1[,1],U_n[,1]),
          CV1_y = cor(uhat_CV1[,1],Y_n[,1]),
          CV1_uhatFull = cor(uhat_CV1[,1],U_n_uhatFull[,1]),
          CV2_u = cor(uhat_CV2,U_n[,1]),
          CV2_y = cor(uhat_CV2,Y_n[,1]),
          CV2_uhatFull = cor(uhat_CV2,U_n_uhatFull[,1]),
          CV2_y_corrected = (cov(uhat_CV2,Y_n[,1]) - cov_correction)/sqrt(var(uhat_CV2) * var(Y_n[,1]))
        )

        results = data.frame(runID,Rep,
                             n_train,n_validation,
                             cor_G, cor_R, cor_R2,
                             H2s=matrix(H2s[1:2],nr=1), # true H2s for calculation corrected estimates
                             results)

      }
    }
  }
}


write.csv(results,file = outfile,row.names = F)

print('done')


