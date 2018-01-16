## SETUP--------------------------------------------------------
library(stats)
library(plyr)
library(nloptr)
answers <- list()
agg_econ <- read.csv("data_agg_econ.csv")

## FUNCTIONS----------------------------------------------------
find_if_maj <- function(var_name, curr_year, k) {
  # ERROR CATCHES FOR INPUTS
  if (!("agg_econ" %in% ls(, envir=.GlobalEnv))) {
    stop("agg_econ DATA NOT IN MEMORY!")
  }
  if (!(var_name %in% colnames(agg_econ))) {
    stop("VARIABLE NOT IN DATA!")
  }
  if (!(curr_year %in% agg_econ$year)) {
    stop("YEAR NOT OBSERVED IN DATA!")
  }
  if (!((curr_year-k+1) %in% agg_econ$year)) {
    stop(c("NOT ENOUGH DATA TO LOOK BEHIND ",curr_year, " BY ",k, " YEARS!"))
  }

  # MAIN ROUTINE
  if (var_name=="unemp_chng") {
    agg_econ$unemp_chng <- agg_econ$unemp_chng*-1
    # because unemp has the oppsite "desirable" moves as the rest of the vars
  }
  var_sub <- agg_econ[agg_econ$year%in%((curr_year-k+1):curr_year),colnames(agg_econ)==var_name]
  # grabs subsetted data necessary for calculations
  return(as.numeric((sum(var_sub>0)/k > .5)))
  # determines if ratio of positives is greater than 1/2
}

find_if_maj_vec <- function(k_vec, curr_year) {
  return(c(find_if_maj("exprt_chng", curr_year, k_vec[1]),
           find_if_maj("unemp_chng", curr_year, k_vec[2]),
           find_if_maj("prices_chng", curr_year, k_vec[3])))
}

is.scalar <- function(x) {is.vector(c(x)) & length(c(x))==1}

find_gammas <- function(X, G_obs, guess=c(.33,.33,.33)) {
  # ERROR CATCHES FOR INPUTS
  if (!(is.scalar(G_obs) & is.vector(X) & length(X)==3)) {
    stop("X AND G DO NOT COVER THE SAME YEARS!")
  } else if ((length(G_obs)>1 & is.matrix(X)) && !(nrow(X)==nrow(G_obs))) {
    stop("X AND G DO NOT COVER THE SAME YEARS!")
  }

  # MAIN ROUTINE
  A <- matrix(c(X[1],X[2],X[3],-X[1],-X[2],-X[3]),2,3, byrow=T)
  # constraint matrix A --> Ax >= b where A is (X1,X2,X3;-X1,-X2,-X3)
  b <- c(0,-1) # b is a vector of (0,-1)
  # constraint vector b --> Ax >= b where b is (0,-1) because X*B (dot prod) must
  # be between 0 and 1 and constraints need to be >= (hence why -1 instead of 1)
  sol <- constrOptim(guess, diff_G, NULL, ui=A, ci=b, X=X, G_obs=G_obs)
  return(sol$par)
}

diff_G <- function(gammas, X, G_obs) {
  if (is.matrix(X)) {
    gammas <- matrix(rep(gammas,nrow(X)),nrow(X),ncol(X))
    ans <- (sum((rowSums(X*gammas) - G_obs)**2))/length(G_obs)
  } else {
    ans <- (((X%*%gammas) - G_obs)**2)/length(G_obs)
  }
  return(ans)
}

compute_X_mtx <- function(k_vec, max_year, list_return=0) {
  min_year <- agg_econ$year[nrow(agg_econ)] + max(k_vec) - 1
  if (list_return) {
    return(list(t(sapply(max_year:min_year, find_if_maj_vec, k_vec=k_vec)), min_year, max_year))
  } else {
    return(t(sapply(max_year:min_year, find_if_maj_vec, k_vec=k_vec)))
  }
}

inside_search <- function(k_vec,G_obs,gamma_guess, max_year) {
  min_year_inside <- agg_econ$year[nrow(agg_econ)] + max(k_vec) - 1
  G_obs_sub_inside <- G_obs[agg_econ$year%in%c(((max_year+1):(min_year+1)))]
  X_mtx <- compute_X_mtx(k_vec, max_year)
  A_inside <- rbind(X_mtx,-X_mtx) # same as above just in multiple years case
  b_inside <- c(rep(0,nrow(X_mtx)),rep(-1,nrow(X_mtx)))
  sol_inside <- constrOptim(gamma_guess, diff_G, NULL, ui=A_inside, ci=b_inside, X=X_mtx, G_obs=G_obs_sub_inside)
  return((sum((rowSums(X*sol_inside$par) - G_obs)**2))/length(G_obs))
}

find_gammas_ks <- function(G_obs, k_guess=c(4,2,3), gamma_guess=c(1/3,1/3,1/3), max_year) {
  # MAIN ROUTINE
  min_year <- agg_econ$year[nrow(agg_econ)] + max(k_guess) - 1
  G_obs_sub <- G_obs[agg_econ$year%in%c((max_year+1):(min_year+1))]
  A_outside <- matrix(c(-1,0,0,0,-1,0,0,0,-1),3,3)
  b_outside <- c(-length(G_obs),-length(G_obs),-length(G_obs)) # k cannot be greater than years in data
  sol_outside <- constrOptim(k_guess, inside_search, NULL, ui=A_outside, ci=b_outside,
                             max_year=max_year, G_obs=G_obs_sub, gamma_guess=gamma_guess)
  return(sol_outside$par)
}

#find_gammas_ks <- function(G_obs, k_guess=c(1,1,1), gamma_guess=c(1/3,1/3,1/3)) {
  # MAIN ROUTINE





## ASSIGNMENT PROBLEMS----------------------------------------
# 2.1
answers$chngs_846 <- find_if_maj_vec(c(8,4,6), 2012)

#2.2
G_obs <- agg_econ$gdp[agg_econ$year==2013]
X <- answers$chngs_846
answers$weights_given <- find_gammas(X,G_obs)

#2.3
min_year <- agg_econ$year[nrow(agg_econ)] + 8 - 1
max_year <- 2015
X_mtx <- t(sapply(max_year:min_year, find_if_maj_vec, k_vec=c(8,4,6)))
answers$X_mtx <- X_mtx

#2.4
gamma_vec <- c(1/3,1/3,1/3)

k_vec_A  <- c(1,1,1)
list_A   <- compute_X_mtx(k_vec_A, max_year=2015, list_return=1)
X_A      <- list_A[[1]]
min_year <- list_A[[2]]
max_year <- list_A[[3]]
gamma_A  <- matrix(rep(gamma_vec,(max_year-min_year) + 1),(max_year-min_year) + 1,3)
mn_sq_A  <- diff_G(gamma_A,X_A,agg_econ$gdp[agg_econ$year%in%((max_year+1):(min_year+1))])

k_vec_B <- c(2,2,3)
list_B   <- compute_X_mtx(k_vec_B, max_year=2015, list_return=1)
X_B      <- list_B[[1]]
min_year <- list_B[[2]]
max_year <- list_B[[3]]
gamma_B  <- matrix(rep(gamma_vec,(max_year-min_year) + 1),(max_year-min_year) + 1,3)
mn_sq_B  <- diff_G(gamma_B,X_B,agg_econ$gdp[agg_econ$year%in%((max_year+1):(min_year+1))])

k_vec_C <- c(12,3,4)
list_C   <- compute_X_mtx(k_vec_C, max_year=2015, list_return=1)
X_C      <- list_C[[1]]
min_year <- list_C[[2]]
max_year <- list_C[[3]]
gamma_C  <- matrix(rep(gamma_vec,(max_year-min_year) + 1),(max_year-min_year) + 1,3)
mn_sq_C  <- diff_G(gamma_C,X_C,agg_econ$gdp[agg_econ$year%in%((max_year+1):(min_year+1))])

k_vec_D <- c(5,7,9)
list_D   <- compute_X_mtx(k_vec_D, max_year=2015, list_return=1)
X_D      <- list_D[[1]]
min_year <- list_D[[2]]
max_year <- list_D[[3]]
gamma_D  <- matrix(rep(gamma_vec,(max_year-min_year) + 1),(max_year-min_year) + 1,3)
mn_sq_D  <- diff_G(gamma_D,X_D,agg_econ$gdp[agg_econ$year%in%((max_year+1):(min_year+1))])

k_vec_E <- c(6,3,4)
list_E   <- compute_X_mtx(k_vec_E, max_year=2015, list_return=1)
X_E      <- list_E[[1]]
min_year <- list_E[[2]]
max_year <- list_E[[3]]
gamma_E  <- matrix(rep(gamma_vec,(max_year-min_year) + 1),(max_year-min_year) + 1,3)
mn_sq_E  <- diff_G(gamma_E,X_E,agg_econ$gdp[agg_econ$year%in%((max_year+1):(min_year+1))])

lett_of_best <- LETTERS[which.min(c(mn_sq_A,mn_sq_B,mn_sq_C,mn_sq_D,mn_sq_E))]
answers$k_bestof5 <- get(paste0("k_vec_",lett_of_best))

#2.5
G_obs <- agg_econ$gdp
max_year <- 2015
