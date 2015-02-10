source("slqr.r")
load("tropcycl.rd")

#install.packages("Rcpp", dependencies = TRUE)
install.packages("RcppArmadillo", dependencies = TRUE)
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("modifiedCode.cpp")
tc.slqr <- slqr(x = Year, y = WmaxST, ry = c(0, 200), nsweep = 1e4)
ss <- 1e3 + seq(9, 9e3, 9)
tau <- seq(.01, .99, .01)
summary(tc.slqr, subset = ss, tau = tau, show.intercept = FALSE)

#To add more sweeps to your MCMC above, you can use
tc.slqr <- update(tc.slqr, nsweep = 1e4)
