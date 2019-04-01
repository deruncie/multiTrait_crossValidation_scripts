# in this simulation, we have a simple family structure to create genetic relatedness
# we measure nT traits on a training population
# we then measure nT traits on a validation population
# The validation population has pairs of individuals.
#   We observe trait 2...nT on 1/2 the indivuals.
#   We use these traits to predict trait 1 on the other 1/2 of the individuals

# In simulations, we modify:
#  n_train: number of training individuals
#  cor_Family: genetic relatedness within families
#  cor_Sib: percent of within-family genetic variation shared by sibs (0-1)
#  cor_G: Genetic correlation of traits
#  cor_R: Residual correlation of traits
#  H21: H2 of trait 1
#  H22: H2 of trait 2


library(lme4)
library(lme4qtl)
library(Matrix)
library(rrBLUP)
library(foreach)
library(doParallel)
library(MASS)
source('my_relmatLmer.R')


runID = as.numeric(commandArgs(t=T)[1])
if(is.na(runID)) runID = 1
registerDoParallel(1)
nReps = 1

cor_Family = .5
cor_Sib = 1
n_Line = as.numeric(commandArgs(t=T)[2])
if(is.na(n_Line)) n_Line = 10
n_Family = 5000/n_Line
n_Sib = 1
nTrait = 2
percent_masked_Lines = .1

outfile = sprintf('Results_1level_big/run_%dLines_%03d.csv',n_Line,runID)
if(file.exists(outfile)) quit('no')


score_function = 'cor'
score = function(preds,actual) {
  if(is.null(preds)) return(NA)
  if(is.null(actual)) return(rep(NA,ncol(preds)))
  preds = as.matrix(preds)
  data.frame(cov = cov(preds,actual), beta = sapply(seq_len(ncol(preds)),function(i) coef(lm(actual ~ preds[,i]))[2]),
             var_uhat = apply(preds,2,var),var_pred = var(actual),
             row.names = NULL)
}

# we have families, and Lines in families, and Sibs in Lines

wide_data = expand.grid(Family = sprintf('F.%03d',1:n_Family),Sib = 1:n_Sib,Line = 1:n_Line)
wide_data$Line = factor(sprintf('%s.L.%02d',wide_data$Family,wide_data$Line))
wide_data$ID = sprintf('%3d',1:nrow(wide_data))
wide_data$y = rnorm(nrow(wide_data))
n = nrow(wide_data)

Z_L = Matrix(model.matrix(~0+Line,wide_data))
colnames(Z_L) = sub('Line','',colnames(Z_L))
Z_L = Z_L[,wide_data$Line]

Line_data = aggregate(y~Line+Family,wide_data,FUN = mean)
Line_data = Line_data[match(wide_data$Line,Line_data$Line),]
Lines = Line_data$Line

Z_family = model.matrix(~0+Family,Line_data)

tall_data = do.call(rbind,lapply(1:nTrait,function(t) data.frame(wide_data,Trait = t)))
tall_data$Trait = factor(tall_data$Trait)

Z_tall = Matrix(model.matrix(~0+Line:Trait,tall_data))
Z_tall = Z_tall[,paste0('Line',c(paste0(Lines,":Trait1"),c(paste0(Lines,":Trait2"))))]

# we mask line 1 from a portion of the families
masked_Lines = as.character(Lines[grep('L.01',Lines)[1:(percent_masked_Lines * length(Lines))]])
length(masked_Lines)
mask = match(masked_Lines,wide_data$Line)
training_Lines = as.character(Lines[as.character(Lines) %in% masked_Lines == F])
Z_training = Z_tall[-mask,]
Z_star = Z_tall[tall_data$Trait == 2 & tall_data$Line %in% masked_Lines,paste0('Line',c(paste0(masked_Lines,":Trait1"),c(paste0(masked_Lines,":Trait2"))))]

# for validation data, we begin with the masked lines.
# we then make a relationship matrix between the validation lines and the observed + masked lines
# we will modify this Kinship matrix by changing the relationship between these lines and the masked_lines
validation_data = wide_data[mask,]

# centering matrix
S = diag(1,length(masked_Lines)) - 1/length(masked_Lines)

# set up genetic matrices
cor_G = 0.6
cor_R = 0.6
H2s = c(.3,rep(.5,nTrait-1))  # broad-sense heritabilities
cor_Sibs = c(0.25,0.5,1)


cor_G = 0.6
cor_R = -0.4
cor_R2 = 0
H21 = 0.4
H22 = 0.7
Rep=1


U_randn = matrix(rnorm(length(Lines)*nTrait),length(Lines),nTrait)
U_validation_randn = rnorm(length(masked_Lines))
E_randn = matrix(rnorm(n*nTrait),n,nTrait)
E_validation_randn = matrix(rnorm(length(masked_Lines)*nTrait),nc=nTrait)

In = Diagonal(n,1)
I_train = Diagonal(n - length(mask),1)


# results = foreach(cor_Family = c(0.1,0.25,0.5,1),.combine = 'rbind',.errorhandling = 'remove') %do% {
results = foreach(cor_Family = c(0.25),.combine = 'rbind',.errorhandling = 'remove') %do% {

  K = cor_Family * tcrossprod(Z_family) + diag((1-cor_Family)*cor_Sib,length(Lines))
  K = K + diag(1e-10,nrow(K))
  rownames(K) = colnames(K) = Lines
  chol_K = chol(K)
  inv_K = Matrix(chol2inv(chol_K))
  rownames(inv_K) = colnames(inv_K) = Lines
  K = Matrix(K,sparse=T)

  K_train = K[training_Lines,training_Lines]
  K_tilde = K[masked_Lines,masked_Lines]
  K_star = K[masked_Lines,training_Lines]
  K_validation = cbind(K_tilde,K_star)

  inv_K_train = solve(K_train)
  K_cond = K_tilde - K_star %*% inv_K_train %*% t(K_star)

  foreach(cor_G = c(0,0.3,0.6),.combine = 'rbind',.errorhandling = 'remove') %do% { #0.2,0.4,
    foreach(cor_R = c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6),.combine = 'rbind',.errorhandling = 'remove') %do% {
      foreach(H21 = c(0.2,.6),.combine = 'rbind',.errorhandling = 'remove') %do% { #0.1,0.4,0.7
        foreach(H22 = c(0.2,0.6),.combine = 'rbind',.errorhandling = 'remove') %do% { #0.4,


          H2s[1] = H21
          H2s[-1] = H22
          Gcor = diag(1,nTrait);
          Gcor[Gcor==0] = cor_G
          G = diag(sqrt(H2s)) %*% Gcor %*% diag(sqrt(H2s))
          Gsqrt = chol(G)
          Rcor = diag(1,nTrait);
          Rcor[Rcor==0] = cor_R2
          Rcor[1,-1] = Rcor[-1,1] = cor_R
          R = diag(sqrt(1-H2s)) %*% Rcor %*% diag(sqrt(1-H2s))
          Rsqrt = chol(R)

          print(c(cor_Family,cor_G,cor_R,H21,H22))
          # set.seed(Rep)
          U = t(chol_K) %*% U_randn %*% Gsqrt
          E = E_randn %*% Rsqrt     # environmental error
          Y = as.matrix(Z_L %*% U + E)

          # fit model to full data
          m0 = lmer(c(Y)~0+Trait+(0+Trait|Family)+(0+Trait|ID),tall_data,control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))

          var1 = as.data.frame(VarCorr(m0))
          Rhat = diag(var1$vcov[7],2)
          Rhat[1,1] = Rhat[1,1] + var1$vcov[1]
          Rhat[2,2] = Rhat[2,2] + var1$vcov[2]
          Rhat[1,2]=Rhat[2,1] = var1$vcov[3]
          Ghat = diag(1,2)
          Ghat[1,1] = var1$vcov[4]
          Ghat[2,2] = var1$vcov[5]
          Ghat[1,2]=Ghat[2,1] = var1$vcov[6]
          Ghat = Ghat/cor_Family
          Rhat = Rhat - Ghat*(1-cor_Family)


          Rhat_full = kronecker(Rhat,In)
          Ghat_full = kronecker(Ghat,K)

          Vhat_full = Rhat_full + Z_tall %*% Ghat_full %*% t(Z_tall)
          uhat = matrix(Ghat_full %*% t(Z_tall) %*% solve(Vhat_full,c(sweep(Y,2,fixef(m0),'-'))),nc = nTrait)
          rownames(uhat) = Lines
          Uhat_full_validation = uhat[masked_Lines,1]


          Y_train = Y
          Y_train[mask,] = NA
          m1 = lmer(c(Y_train)~0+Trait+(0+Trait|Family)+(0+Trait|ID),tall_data,control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))

          var1 = as.data.frame(VarCorr(m1))
          Rhat = diag(var1$vcov[7],2)
          Rhat[1,1] = Rhat[1,1] + var1$vcov[1]
          Rhat[2,2] = Rhat[2,2] + var1$vcov[2]
          Rhat[1,2]=Rhat[2,1] = var1$vcov[3]
          Ghat = diag(1,2)
          Ghat[1,1] = var1$vcov[4]
          Ghat[2,2] = var1$vcov[5]
          Ghat[1,2]=Ghat[2,1] = var1$vcov[6]
          Ghat = Ghat/cor_Family
          Rhat = Rhat - Ghat*(1-cor_Family)

          h2_hat_joint = Ghat[1,1]/(Ghat[1,1]+Rhat[1,1])

          Rhat_full = kronecker(Rhat,I_train)
          Ghat_full = kronecker(Ghat,K_train)
          Vhat_full = Rhat_full + Ghat_full
          uhat = matrix(Ghat_full %*% solve(Vhat_full,c(sweep(Y_train[-mask,],2,fixef(m1),'-'))),nc = nTrait)

          # now validation lines conditional on Training lines
          ustar_mean = K_star %*% inv_K_train %*% uhat
          Sigma_star = Ghat[2,2]*K_cond + diag(Rhat[2,2],length(mask))
          Sigma_star_inv = solve(Sigma_star)
          What_star = Ghat[1,2]*K_cond %*% Sigma_star_inv
          Uhat_2step_validation = as.matrix(ustar_mean[,1] + What_star %*% (Y[mask,2] - fixef(m1)[2] - ustar_mean[2]))

          Uhat = matrix(c(Uhat_2step_validation,matrix(uhat,nc = nTrait)[,1]),nc=1)

          cov_correction = Rhat[1,2]*sum(diag(S %*% What_star)) / (length(masked_Lines)-1)

          # make predictions for validation data without y2_validation depending on cor_Sibs
          Uhat_validation_CV1 =  as.matrix(K_validation[,training_Lines] %*% inv_K_train %*% matrix(uhat,nc=nTrait)[,1])

          # make predictions for validation data depending on cor_Sibs
          Uhat_validation = foreach(cor_Sib = cor_Sibs,.combine = cbind) %do% {
            K_validation_i = K_validation
            # diag(K_validation_i[masked_Lines,masked_Lines]) = cor_Family + (1-cor_Family) * cor_Sib
            diag(K_validation_i[masked_Lines,masked_Lines]) = cor_Sib
            as.matrix(K_validation_i %*% inv_K[Lines,Lines] %*% Uhat[,1])
          }

          # make new "True" Validation data based on cor_Sibs
          # U_validation is the conditional distribution given U and the relationship between the validation and training individuals
          cov_A_observed_actual_inv = kronecker(solve(G),inv_K[Lines,Lines])
          U_rot = cov_A_observed_actual_inv %*% Z_tall %*% c(U)

          U_validation = foreach(cor_Sib = cor_Sibs,.combine = cbind) %do% {
            K_validation_i = K_validation
            # diag(K_validation_i[masked_Lines,masked_Lines]) = cor_Family + (1-cor_Family) * cor_Sib
            diag(K_validation_i[masked_Lines,masked_Lines]) = cor_Sib
            cov_A_validation = kronecker(G[1,,drop=FALSE],K_validation_i)
            Sigma = diag(G[1,1],length(masked_Lines)) - cov_A_validation %*% cov_A_observed_actual_inv %*% t(cov_A_validation)
            if(min(diag(Sigma)) < 0) diag(Sigma) = diag(Sigma) + 1e-10 - pmin(0,min(diag(Sigma)))
            chol_Sigma = chol(Sigma,pivot = F)
            U_validation_mean = cov_A_validation %*% U_rot
            U_validation_mean + t(chol_Sigma) %*% U_validation_randn
          }

          # make new "Ys" for validation individuals
          # E_validation is uncorrelated with E, and is the same across cor_Sibs
          E_validation = E_validation_randn %*% Rsqrt[,1]

          Y_validation = U_validation + c(E_validation)

          res = lmer(Y_train[,1]~(1|Family),data = wide_data)
          res_vars = as.data.frame(VarCorr(res))
          s2_g = res_vars$vcov[1]/cor_Family
          s2_r = res_vars$vcov[2] - s2_g*(1-cor_Family)
          h2_hat_single = s2_g/(s2_g+s2_r)

          V = s2_g * K_train + s2_r * I_train
          Uhat_validation_single = s2_g * as.matrix(K_validation[,training_Lines] %*% solve(V,Y_train[-mask,1]-fixef(res)))[,1]


      # Results

      results = rbind(
        data.frame(Method = 'Joint', Data = 'U_Sib', cor_Sib = cor_Sibs, cov_correction = 0, h2_correction = 1,
                   do.call(rbind,lapply(1:length(cor_Sibs),function(i) score(Uhat_validation[,i],U_validation[,i])))),
        data.frame(Method = 'Joint', Data = 'Y_Sib', cor_Sib = cor_Sibs, cov_correction = 0, h2_correction = h2_hat_joint,
                   do.call(rbind,lapply(1:length(cor_Sibs),function(i) score(Uhat_validation[,i],Y_validation[,i])))),
        data.frame(Method = 'Joint', Data = 'U_training', cor_Sib = 1, cov_correction = 0, h2_correction = 1,
                   do.call(rbind,lapply(length(cor_Sibs),function(i) score(Uhat_validation[,i],U[masked_Lines,1])))),
        data.frame(Method = 'Joint', Data = 'Y_training', cor_Sib = 1, cov_correction = cov_correction, h2_correction = h2_hat_joint,
                   do.call(rbind,lapply(length(cor_Sibs),function(i) score(Uhat_validation[,i],Y[mask,1])))),
        data.frame(Method = 'Joint', Data = 'Uhat_joint', cor_Sib = 1, cov_correction = 0, h2_correction = 1,
                   do.call(rbind,lapply(length(cor_Sibs),function(i) score(Uhat_validation[,i],Uhat_full_validation)))),
        data.frame(Method = 'Single', Data = 'U_Sib', cor_Sib = cor_Sibs, cov_correction = 0, h2_correction = 1,
                   do.call(rbind,lapply(1:length(cor_Sibs),function(i) score(Uhat_validation_single,U_validation[,i])))),
        data.frame(Method = 'Single', Data = 'Y_Sib', cor_Sib = cor_Sibs, cov_correction = 0, h2_correction = h2_hat_single,
                   do.call(rbind,lapply(1:length(cor_Sibs),function(i) score(Uhat_validation_single,Y_validation[,i])))),
        data.frame(Method = 'Single', Data = 'U_training', cor_Sib = 1, cov_correction = 0, h2_correction = 1,
                   score(Uhat_validation_single,U[masked_Lines,1])),
        data.frame(Method = 'Single', Data = 'Y_training', cor_Sib = 1, cov_correction = 0, h2_correction = h2_hat_single,
                   score(Uhat_validation_single,Y[mask,1])),
        data.frame(Method = 'Single', Data = 'Uhat_joint', cor_Sib = 1, cov_correction = 0, h2_correction = 1,
                   score(Uhat_validation_single,Uhat_full_validation)),
        data.frame(Method = 'CV1', Data = 'U_training', cor_Sib = 1, cov_correction = 0, h2_correction = 1,
                   score(Uhat_validation_CV1,U[masked_Lines,1])),
        data.frame(Method = 'CV1', Data = 'Y_training', cor_Sib = 1, cov_correction = 0, h2_correction = h2_hat_joint,
                   score(Uhat_validation_CV1,Y[mask,1])),
        data.frame(Method = 'CV1', Data = 'Uhat_joint', cor_Sib = 1, cov_correction = 0,h2_correction = 1,
                   score(Uhat_validation_CV1,Uhat_full_validation))
      )
      results = data.frame(runID,Rep,
                           n_Family, n_Line, nTrait, n_Sib,
                           percent_masked_Lines,
                           cor_Family, cor_Sib,
                           cor_G, cor_R, cor_R2, H2s=matrix(H2s[1:2],nr=1),
                           Ghat = matrix(c(Ghat),nr=1),
                           Rhat = matrix(c(Rhat),nr=1),
                           score_function,
                           results)
      }
      }}}
}

write.csv(results,file = outfile,row.names = F)

print('done')

# collect results
# library(foreach)
# library(data.table)
# library(doParallel)
# files = list.files(path='Results_1level_big',pattern='.csv',full.names=T)
# registerDoParallel(10)
# results = foreach(file = files,.combine = 'rbind') %dopar% {
#   print(match(file,files)/length(files))
#   fread(file,data.table=F)
# }
# saveRDS(results,file = 'Results_1level_big/collected_results_1level_big.rds')

