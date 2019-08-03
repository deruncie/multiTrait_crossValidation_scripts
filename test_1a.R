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

# prepare K matrix
load('FileS3.rdata')


library(lme4)
library(lme4qtl)
library(Matrix)
library(rrBLUP)
library(foreach)
library(doParallel)
library(MASS)
library(BGLR)

# set up run
runID = as.numeric(commandArgs(t=T)[1])
if(is.na(runID)) runID = 1
registerDoParallel(1)
nReps = 1

# get K matrix
load('FileS3.rdata')
K = G
rm(list=c('G','Y'))

n = nrow(K)
rownames(K) = colnames(K) = 1:nrow(K)

#------------------------------------------------------------#
# Test 1.a.
#------------------------------------------------------------#
try(dir.create('Results'))
try(dir.create('Results/Raw'))
outfile = sprintf('Results/test_1a.csv',runID)
# if(file.exists(outfile)) quit('no')

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
# X = kronecker(diag(1,2),t(Q)) %*% model.matrix(~0+Trait,tall_data_full)

K_nn = K[validation,validation]
svd_K_nn = svd(K_nn)
Q_nn = svd_K_nn$u
D_K_nn = Diagonal(n_validation,svd_K_nn$d)
rownames(D_K_nn) = colnames(D_K_nn) = wide_data_validation$Line
# X_validation = kronecker(diag(1,2),t(Q_nn)) %*% model.matrix(~0+Trait,tall_data_validation)

K_oo = K[training,training]
svd_K_oo = svd(K_oo)
Q_oo = svd_K_oo$u
D_K_oo = Diagonal(n_train,svd_K_oo$d)
rownames(D_K_oo) = colnames(D_K_oo) = wide_data_train$Line
# X_train = kronecker(diag(1,2),t(Q_oo)) %*% model.matrix(~0+Trait,tall_data_train)

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


registerDoParallel(detectCores())
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

        # analytical calculations
        V_o1 = G[1,1]*K_oo + R[1,1]*I_oo
        V_o1_inv = solve(V_o1)
        E_num_single = G[1,1]^2*sum(diag(S %*% K_no %*% V_o1_inv %*% t(K_no)))
        E_cov_single = E_num_single / (n_validation-1)
        E_var_single = G[1,1]^2*sum(diag(S %*% K_no %*% solve(V_o1,G[1,1]*K_oo+R[1,1]*I_oo) %*% V_o1_inv %*% t(K_no))) / (n_validation-1)

        V_o = kronecker(G,K_oo) + kronecker(R,I_oo)
        V_o_inv = solve(V_o)
        E_num_CV1 = sum(diag(S %*% kronecker(G[1,,drop=F],K_no) %*% V_o_inv %*% kronecker(G[,1,drop=F],t(K_no))))
        E_cov_CV1 = E_num_CV1 / (n_validation-1)
        # E_var_CV1 = sum(diag(S %*% kronecker(G[1,,drop=F],K_no) %*% solve(V_o,kronecker(G,K_oo) + kronecker(R,I_oo)) %*% solve(V_o) %*% kronecker(G[,1,drop=F],t(K_no)))) / (n_validation-1)
        E_var_CV1 = sum(diag(S %*% kronecker(G[1,,drop=F],K_no) %*% V_o_inv %*% kronecker(G[,1,drop=F],t(K_no)))) / (n_validation-1)

        V_c = G[2,2]*K_c + R[2,2]*I_nn
        V_c_inv = solve(V_c)
        E_num_CV2_u = E_num_CV1 -
          sum(diag(S %*% kronecker(G[1,2,drop=F],K_c) %*% V_c_inv %*% kronecker(G[2,,drop=F],K_no) %*% V_o_inv %*% kronecker(G[,1,drop=F],t(K_no)))) +
          sum(diag(S %*% kronecker(G[1,2,drop=F],K_c) %*% V_c_inv %*% kronecker(G[2,1,drop=F],K_nn)))
        E_cov_CV2_u = E_num_CV2_u / (n_validation-1)

        A = kronecker(G[1,,drop=F],K_no) %*% V_o_inv
        C = kronecker(G[1,2,drop=F],K_c) %*% V_c_inv
        B = C %*% kronecker(G[2,,drop=F],K_no) %*% V_o_inv

        E_var_CV2 = sum(diag(S %*% (A-B) %*% (kronecker(G,K_oo) + kronecker(R,I_oo)) %*% t(A-B) +
                               C %*% (G[2,2]*K_nn + R[2,2]*I_nn) %*% t(C) +
                               2*(A-B) %*% (kronecker(G[,2,drop=F],t(K_no))) %*% t(C))) / (n_validation-1)

        E_num_CV2_y = E_num_CV2_u + sum(diag(S %*% kronecker(G[1,2,drop=F],K_c) %*% V_c_inv %*% kronecker(R[2,1,drop=F],I_nn)))
        E_cov_CV2_y = E_num_CV2_y / (n_validation-1)

        E_var_u = sum(diag(S %*% kronecker(G[1,1],K_nn))) / (n_validation-1)
        E_var_y = E_var_u + sum(diag(S %*% kronecker(R[1,1],I_nn))) / (n_validation-1)

        results = data.frame(
          E_var_u,
          E_var_y,
          E_var_single,
          E_var_CV1,
          E_var_CV2,
          E_cov_single,
          E_cov_CV1,
          E_cov_CV2_u,
          E_cov_CV2_y,
          E_cor_single_u = E_cov_single/sqrt(E_var_u * E_var_single),
          E_cor_single_y = E_cov_single/sqrt(E_var_y * E_var_single),
          E_cor_CV1_u = E_cov_CV1 / sqrt(E_var_u*E_var_CV1),
          E_cor_CV1_y = E_cov_CV1 / sqrt(E_var_y*E_var_CV1),
          E_cor_CV2_u = E_cov_CV2_u / sqrt(E_var_u*E_var_CV2),
          E_cor_CV2_y = E_cov_CV2_y / sqrt(E_var_y*E_var_CV2)
        )

        results = data.frame(runID,Rep,
                             n_train,n_validation,
                             cor_G, cor_R, cor_R2, H2s=matrix(H2s[1:2],nr=1),
                             # score_function,
                             results)

      }
    }
  }
}

write.csv(results,file = outfile,row.names = F)

print('done')





