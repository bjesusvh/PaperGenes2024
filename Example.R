############################ Chunk 1 ############################
library(BGLR); library(MPS)
setwd("Put you directory here")
geno <- readRDS("./data/X.rds")
pheno <- readRDS("./data/Y.rds")
set.seed(63) # for reproducibility
id_candidates <- sample(x = 1:50, size = 50, replace = FALSE)
Y <- pheno[-id_candidates,]
X <- geno[-id_candidates,]
Xcandidates <- geno[id_candidates,]
no_iter <- 50e3; no_burn <- 20e3; thin <- 5
id_samples <- which(seq(1, no_iter, thin) > no_burn)
ETA <- list(list(X = X, model='BRR', saveEffects = TRUE))

############################ Chunk 2 ############################
# Fit Multitrait regression model on quantitative traits (traits 1,2,3)
fmQ <- Multitrait(y = Y[,1:3], 
                  ETA = ETA, intercept = TRUE,   
                  resCov = list(type = "UN", saveEffects = TRUE), 
                  nIter = no_iter, burnIn = no_burn, thin = thin, 
                  verbose = FALSE, saveAt = "./out/E1/")

# Fit ordinal regression on categorical trait (trait 4)
fmO1 <- BGLR(Y[,4], response_type = 'ordinal',
            ETA = ETA, nIter = no_iter, burnIn = no_burn, 
            thin = thin, verbose = FALSE, saveAt = "./out/E1/O1")

# Fit ordinal regression on categorical trait (trait 5)
fmO2 <- BGLR(Y[,5], response_type = 'ordinal',
             ETA = ETA, nIter = no_iter, burnIn = no_burn, 
             thin = thin, verbose = FALSE, saveAt = "./out/E1/O2")

############################ Chunk 3 ############################
# Read MCMC samples from Multitrait model
muQ <- as.matrix(read.table(file = "./out/E1/mu.dat", header = FALSE))[id_samples, ]
betasQ <- readBinMatMultitrait("./out/E1/ETA_1_beta.bin")[id_samples, ,]  
covQ <- as.matrix(read.table(file = "./out/E1/R.dat", header = FALSE))[id_samples, ]

# Read MCMC samples from ordinal model 1
betasO1 <- readBinMat(paste0(file.path("./out/E1/O1ETA_1_b.bin")))         

# Read MCMC samples from ordinal model 2
betasO2 <- readBinMat(paste0(file.path("./out/E1/O2ETA_1_b.bin")))     
n <- nrow(Y)            # number of lines
M <- length(id_samples) # number of MCMC after burnIn and thining
k <- ncol(X)   # number of SNPs
t <- ncol(Y)   # total number of traits

# Combine mean of Quantitative and ordinal
mu <- as.matrix(cbind(muQ, 0, 0)) 
# Fill the array of regression coefficients
betas <- array(NA , c(M, k, t))
betas[ , ,c(1,2,3)] <- betasQ
betas[,,4] <-  betasO1
betas[,,5] <-  betasO2
# Create covariance matrix
covMatrix <- as.matrix(cbind(covQ[,1], covQ[,2], covQ[,3], 0, 0, 
                             covQ[,4], covQ[,5], 0, 0, covQ[,6], 0, 0, 1, 0, 1))
out <- FastMPS(Xcand = Xcandidates, B0 = mu, B = betas, R = covMatrix)
