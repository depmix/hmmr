\name{stan_2stateHMM_speed}
\alias{stan_2stateHMM_speed}
\docType{data}
\title{
stan 2stateHMM speed
}
\description{
Posterior samples from a two state HMM for the speed data with rstan 
}
\usage{data("stan_2stateHMM_speed")}
\format{
  A \code{stanfit} object
}
\details{
The estimation was done with rstan version 2.17.3 through the following code:

\preformatted{library(depmixS4)
library(hmmr)
library(rstan)

data(speed)

set.seed(4321)

stanmodelcode <- "
data {
  // counters
  int<lower=2> nstates; // number of states
  int<lower=0> T; // total observations
  int<lower=0> nrep; // number of replications 
  int<lower=0,upper=T> bt[nrep]; // start times for replications
  int<lower=0,upper=T> et[nrep]; // end times for replications
  // data
  real rt[T];
  int<lower=0,upper=1> corr[T];
  real<lower=0,upper=1> Pacc[T];
  // hyper-parameters
  vector<lower=0>[nstates] alpha_init;  // init prior
  vector<lower=0>[nstates] alpha_trans;  // transit prior
  real<lower=0> alpha_obs; // emit prior
  real<lower=0> beta_obs; // emit prior
}
parameters {
  simplex[nstates] init; // initial state probabilities
  simplex[nstates] transition[nstates]; // transition matrices
  real<lower=0,upper=1> theta[nstates]; // p(correct)|state
  ordered[nstates] mu; // means
  real<lower=0> sigma[nstates]; // standard deviations
}
model {
  real acc[nstates];
  real gamma[T,nstates];

  // prior for initial state probabilities
  init ~ dirichlet(alpha_init);
  // prior for transition probabilities
  for (s in 1:nstates) 
    transition[s] ~ dirichlet(alpha_trans);
  // prior for probability of a correct response
  for (s in 1:nstates)
    theta[s] ~ beta(alpha_obs,beta_obs);
  // prior for Normal logRT distributions 
  for (s in 1:nstates) {
    mu[s] ~ normal(0,10);
    sigma[s] ~ cauchy(0,10);
  }
  // forward-backward
  for (i in 1:nrep) {
    // first t
    for (s in 1:nstates) {
      gamma[bt[i],s] = log(init[s]) + bernoulli_lpmf(corr[bt[i]] | theta[s]) + normal_lpdf(rt[bt[i]]|mu[s],sigma[s]);
    }
    // remaining t
    for (t in (bt[i]+1):et[i]) {
      for (s in 1:nstates) {
        for (u in 1:nstates) {
          acc[u] = gamma[t-1,u] + log(transition[u,s]) + bernoulli_lpmf(corr[t] | theta[s]) + normal_lpdf(rt[t]|mu[s],sigma[s]);
        }
        gamma[t,s] = log_sum_exp(acc);
      }
    }
    target += log_sum_exp(gamma[et[i]]);
  }
}
"

standata <- list(
  nstates = 2,
  T = nrow(speed),
  nrep = 3,
  bt = c(1,169,303),
  et = c(168,302,439),
  rt = speed$rt,
  corr = as.numeric(speed$corr == "cor"),
  Pacc = speed$Pacc,
  alpha_init = c(1,1),
  alpha_trans = c(1,1),
  alpha_obs = 1,
  beta_obs = 1
)

stan_2stateHMM_speed <- stan(model_name = "two-state-HMM",
                model_code = stanmodelcode,
                data = standata,
                chains = 4,
                iter = 4000,
                warmup = 2000)
}
}
\examples{
data(stan_2stateHMM_speed)
## maybe str(stan_2stateHMM_speed) ; plot(stan_2stateHMM_speed) ...
}
\keyword{datasets}
