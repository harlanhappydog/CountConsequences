

# 
myargs= as.numeric(commandArgs(trailingOnly=TRUE))
print(myargs)

iii <- myargs

# all_pkgs <- c(names(sessionInfo()$loadedOnly), names(sessionInfo()$otherPkgs))

# while(length(all_pkgs)){
# lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)

# if(length(names(sessionInfo()$otherPkgs))){
# invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
# }
# all_pkgs <- c(names(sessionInfo()$loadedOnly), names(sessionInfo()$otherPkgs))
# }


# ls()
# names(sessionInfo()$loadedOnly)
# names(sessionInfo()$otherPkgs)

# Determine missing packages and load them:
# required_packages <- c("MASS", "pscl", "lmtest", "stats", "countreg", "MuMIn")
# not_installed <- required_packages[!(required_packages %in%
                   # installed.packages()[, "Package"])]
# if (length(not_installed)) {
    # install.packages(not_installed)
   # }
# suppressWarnings(lapply(required_packages, require, character.only = TRUE))


required_packages <- c("MASS", "pscl", "lmtest", "stats", "countreg", "MuMIn")
suppressWarnings(lapply(required_packages, require, character.only = TRUE))

set.seed(12345)
#################################################
## Custom functions:

rzinb <-function (n, size, mu, rho) 
{
    x <- ifelse(rbinom(n, 1, rho), 0, rnbinom(n, size, mu = mu))
    return(x)
}

#######
## vuong_f test :  compare model1 to model2
#######

vuong_f<-function (m1, m2, digits = getOption("digits")) 
{
    m1y <- m1$y
    m2y <- m2$y
    m1n <- length(m1y)
    m2n <- length(m2y)
    if (m1n == 0 | m2n == 0) 
        stop("Could not extract dependent variables from models.")
    if (m1n != m2n) 
        stop(paste("Models appear to have different numbers of observations.\n", 
            "Model 1 has ", m1n, " observations.\n", "Model 2 has ", 
            m2n, " observations.\n", sep = ""))
    if (any(m1y != m2y)) {
        stop(paste("Models appear to have different values on dependent variables.\n"))
    }
    p1 <- predprob(m1)
    p2 <- predprob(m2)
    if (!all(colnames(p1) == colnames(p2))) {
        stop("Models appear to have different values on dependent variables.\n")
    }
    whichCol <- match(m1y, colnames(p1))
    whichCol2 <- match(m2y, colnames(p2))
    if (!all(whichCol == whichCol2)) {
        stop("Models appear to have different values on dependent variables.\n")
    }
    m1p <- rep(NA, m1n)
    m2p <- rep(NA, m2n)
    for (i in 1:m1n) {
        m1p[i] <- p1[i, whichCol[i]]
        m2p[i] <- p2[i, whichCol[i]]
    }
    k1 <- length(coef(m1))
    k2 <- length(coef(m2))
    lm1p <- log(m1p)
    lm2p <- log(m2p)
    m <- lm1p - lm2p
    bad1 <- is.na(lm1p) | is.nan(lm1p) | is.infinite(lm1p)
    bad2 <- is.na(lm2p) | is.nan(lm2p) | is.infinite(lm2p)
    bad3 <- is.na(m) | is.nan(m) | is.infinite(m)
    bad <- bad1 | bad2 | bad3
    neff <- sum(!bad)
    if (any(bad)) {
        cat("NA or numerical zeros or ones encountered in fitted probabilities\n")
        cat(paste("dropping these", sum(bad), "cases, but proceed with caution\n"))
    }
    aic.factor <- (k1 - k2)/neff
    bic.factor <- (k1 - k2)/(2 * neff) * log(neff)
    v <- rep(NA, 3)
    arg1 <- matrix(m[!bad], nrow = neff, ncol = 3, byrow = FALSE)
    arg2 <- matrix(c(0, aic.factor, bic.factor), nrow = neff, 
        ncol = 3, byrow = TRUE)
    num <- arg1 - arg2
    s <- apply(num, 2, sd)
    numsum <- apply(num, 2, sum)
    v <- numsum/(s * sqrt(neff))
    names(v) <- c("Raw", "AIC-corrected", "BIC-corrected")
    pval <- rep(NA, 3)
    msg <- rep("", 3)
    for (j in 1:3) {
        if (v[j] > 0) {
            pval[j] <- 1 - pnorm(v[j])
            msg[j] <- "model1 > model2"
        }
        else {
            pval[j] <- pnorm(v[j])
            msg[j] <- "model2 > model1"
        }
    }
    out <- data.frame(v, msg,(pval))
    names(out) <- c("Vuong z-statistic", "H_A", "p-value")

    return(out)
}
#################################################

#######
## simstudy:  function to conduct the simulation study
#######


simstudy<-function(phi=1, n=50, beta0=0.4, p_ZI=0, nSim=3000, alpha_level=0.05){

# set-up for the simulation
results_mat <- matrix(0, nSim, 20)

# X is the fixed design matrix of covariates
# X <- cbind(1,sample(c(0,1),n, replace=TRUE))

for(jj in 1:nSim) {
	
X <- cbind(1,rnorm(n,0,10))

if(jj/41==round(jj/41)){print(jj)}

betavec <- c(beta0, 0); xvar <- X[,-1]

lambda <- exp(X%*%betavec)
omega <- p_ZI
nu <- lambda*phi

y <- rzinb(n, size = nu, mu = lambda, rho = omega)

#####################
### Models:


#####################
### poisson model ###
poissonglm <- glm(y ~ xvar, family = poisson)
poisson_pval <- waldtest(poissonglm, test = "Chisq")["Pr(>Chisq)"][2, ]
poisson_aic <- AIC(poissonglm)
poisson_aicc <- AICc(poissonglm)
poisson_bic <- BIC(poissonglm)

c(poisson_pval, poisson_aic, poisson_aicc, poisson_bic)

#####################
### NB model ###
nb_aic <- Inf
nb_aicc <- Inf
nb_bic <- Inf
nb_pval <- 0.99
tryCatch({
    nb_mod <- glm.nb(y ~ xvar)
    nb_pval <- waldtest(nb_mod, test = "Chisq")["Pr(>Chisq)"][2, ]
    nb_aic <- AIC(nb_mod)
    nb_aicc <- AICc(nb_mod)
    nb_bic <- BIC(nb_mod)
  },
  error = function(e){}
)

c(nb_pval, nb_aic, nb_aicc, nb_bic)

#####################
### ZIP model ###
zip_aic <- Inf
zip_aicc <- Inf
zip_bic <- Inf
zip_mod <- 99
zip_pval <- 0.99
zip_modna <- TRUE
tryCatch({
    zip_mod <- zeroinfl(y ~ xvar | xvar, dist = "poiss")
    zip_modna <- sum(is.na(unlist((summary(zip_mod))$coefficients))) > 0
  },
  error = function(e) {"e"}
)

if (is.double(zip_mod) | zip_modna) {
  tryCatch(
    {
      zip_mod <- zeroinfl(y ~ xvar | xvar, dist = "poiss", EM = TRUE)
    },
    error = function(e){"e"}
    )
}

tryCatch(
  {
    zip_pval <- waldtest(zip_mod, test = "Chisq")["Pr(>Chisq)"][2, ]
    zip_aic <- AIC(zip_mod)
    zip_aicc <- AICc(zip_mod)
    zip_bic <- BIC(zip_mod)
  },
  error = function(e){"e"}
)

c(zip_pval, zip_aic, zip_aicc, zip_bic)

#####################
### ZINB model ###
zinb_aic <- Inf
zinb_aicc <- Inf
zinb_bic <- Inf
zinb_mod <- 99
zinb_pval <- 0.99
zinb_modna <- TRUE

tryCatch(
  {
    zinb_mod <- zeroinfl(y ~ xvar | xvar, dist = "negbin")
    zinb_modna <- sum(is.na(unlist((summary(zinb_mod))$coefficients))) > 0
  },
  error = function(e){"e"}
)

if (is.double(zinb_mod) | zinb_modna) {
  tryCatch(
    {
      zinb_mod <- zeroinfl(y ~ xvar | xvar, dist = "negbin", EM = TRUE)
      summary(zinb_mod)
    },
    error = function(e){"e"}
  )
}

tryCatch(
  {
    zinb_pval <- waldtest(zinb_mod, test = "Chisq")["Pr(>Chisq)"][2, ]
    dimK <- dim(summary(zinb_mod)$coefficients$count)[1]
    zinb_aic <- AIC(zinb_mod)
    zinb_aicc <- AICc(zinb_mod)    
    zinb_bic <- BIC(zinb_mod)    
  },
  error = function(e){"e"}
)

c(zinb_pval, zinb_aic, zinb_aicc, zinb_bic)

#####################
### Tests:


###########
### the D&L score test for overdispersion:
lambdahat <- yhat <- predict(poissonglm, type = "response")
T_1 <- sum((y - lambdahat)^2 - y) / sqrt(2 * sum(lambdahat^2))
LRT_pval <- DLtest_pval <- pnorm(T_1, lower.tail = FALSE)

###########
### the Vuong test for zero-inflation:
vuong_P_zip_pval <- vuong_NB_zinb_pval <- 1

if (exists("zip_mod") & sum(y == 0) > 1) {
  tryCatch(
    {
      vv_P_ZIP <- (vuong_f(poissonglm, zip_mod))
      vuong_P_zip_pval <- as.numeric(as.character(vv_P_ZIP[1, 3]))
    },
    error = function(e) {}
  )
}

if (exists("zinb_mod") & sum(y == 0) > 1) {
  tryCatch(
    {
      vv_NB_ZINB <- (vuong_f(nb_mod, zinb_mod))
      vuong_NB_zinb_pval <- as.numeric(as.character(vv_NB_ZINB[1, 3]))
    },
    error = function(e) {}
  )
}

c(LRT_pval, vuong_P_zip_pval, vuong_NB_zinb_pval)


# Step 2: If both score tests fail to reject the null...
if( (LRT_pval > alpha_level) & (vuong_P_zip_pval > alpha_level) ){
pval1 <- poisson_pval	
choice <- 1
}

# Step 3: If the DL score test rejects the null 
# and the vdB score test fails to reject the null...
if( (LRT_pval <=alpha_level) & (vuong_NB_zinb_pval > alpha_level) ){		
pval1 <- nb_pval
choice <- 2
}


# Step 4: If the DL score test fails to reject the null 
# and the vdB score test does reject the null...
if( (LRT_pval > alpha_level) & (vuong_P_zip_pval <= alpha_level) ){
pval1 <- zip_pval
choice <- 3
}

# Step 5: If both the  DL score test and the vdB reject the null...
if( (LRT_pval <=  alpha_level) & (vuong_NB_zinb_pval <= alpha_level) ){
pval1 <- zinb_pval
choice <- 4
}

pvalAIC <- c(poisson_pval, nb_pval, zip_pval, zinb_pval)[
which.min(c(poisson_aic , nb_aic , zip_aic, zinb_aic))]

pvalAICc <- c(poisson_pval, nb_pval, zip_pval, zinb_pval)[
which.min(c(poisson_aicc , nb_aicc , zip_aicc, zinb_aicc))]

pvalBIC <- c(poisson_pval, nb_pval, zip_pval, zinb_pval)[
which.min(c(poisson_bic , nb_bic , zip_bic, zinb_bic))]


choiceAIC <- which.min(c(poisson_aic , nb_aic , zip_aic , zinb_aic))
choiceAICc <- which.min(c(poisson_aicc , nb_aicc , zip_aicc , zinb_aicc))
choiceBIC <- which.min(c(poisson_bic , nb_bic , zip_bic , zinb_bic))

# save results to results matrix:
results_mat[jj,] <- c(n, beta0, phi, p_ZI, nSim, LRT_pval, vuong_P_zip_pval, vuong_NB_zinb_pval, poisson_pval, nb_pval, zip_pval, zinb_pval, pval1, pvalAIC, pvalAICc, pvalBIC, choice, choiceAIC, choiceAICc, choiceBIC)
}

colnames(results_mat)<-c("n", "beta_0", "phi", "p_ZI", "nSim", "LRT_pval", "vuong_P_zip_pval", "vuong_NB_zinb_pval", "poisson_pval", "nb_pval", "zip_pval", "zinb_pval", "pval1", "pvalAIC", "pvalAICc", "pvalBIC", "choice", "choiceAIC", "choiceAICc",  "choiceBIC")


return(results_mat)
} 
 

######### ######### ######### ######### ######### #########

 

######### ######### ######### ######### ######### #########


lvls <- list()

lvls[[1]] <- c(50, 100, 250, 500, 1000, 2000, 5000, 10000)		# n
lvls[[2]] <- c(0.5, 1, 1.5, 2, 2.5)					# beta_0
lvls[[3]] <- c(Inf, 4/2, 1, 1/2, 1/3)	 			# phi
lvls[[4]] <- c(0, 0.05, 0.1, 0.2, 0.5)				# omega

dsgn <- as.matrix(expand.grid(lvls[[1]], lvls[[2]], lvls[[3]], lvls[[4]]))

colnames(dsgn)<-c("n", "beta_0", "phi", "omega")

dim(dsgn)

results_list<-list()

print(c("iii", iii))
 
ls()

results_list <- simstudy(
   phi = dsgn[iii,3], 
   n = dsgn[iii,1], 
   beta0 = dsgn[iii,2],
   p_ZI = dsgn[iii,4], 
   alpha_level = 0.05, 
   nSim = 15000)
   
ls()
   
print(head(results_list))

saveRDS(results_list, file=paste(paste(
 "results_Nov2020_15000", iii, sep="_"),".rds",sep=""))

warnings()