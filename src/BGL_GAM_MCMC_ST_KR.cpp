#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

const double log2pi = std::log(2.0 * 3.14159265359);

// MVN density.

// [[Rcpp::export]]
double lk_Chol_arma(arma::colvec x_, arma::colvec mu_,
                    arma::mat Sigma_, bool logd = true){
  int k = x_.n_rows;
  double out;
  arma::mat Chol_Sigma = trimatu(arma::chol(Sigma_));
  arma::vec rooti = arma::solve(Chol_Sigma.t(), x_ - mu_);
  double quads = arma::dot(rooti, rooti);
  out = -(k/2.0)*log2pi - 0.5 * quads -
    sum(log(Chol_Sigma.diag()));
  if (logd == false) {out = exp(out);}
  return(out);
}

// To generate inverse gaussian random variables.

// [[Rcpp::export]]
NumericVector rinvgauss_rcpp(int n, double mu, double lambda){
  
  NumericVector random_vector(n);
  double z,y,x,u;
  
  for(int i=0; i<n; ++i){
    z=R::rnorm(0,1);
    y=z*z;
    x=mu+0.5*mu*mu*y/lambda - 0.5*(mu/lambda)*sqrt(4*mu*lambda*y+mu*mu*y*y);
    u=R::runif(0,1);
    if(u <= mu/(mu+x)){
      random_vector(i)=x;
    }else{
      random_vector(i)=mu*mu/x;
    };
  }
  return(random_vector);
}

// logSumExp (accurate log(sum(exp())) for small values).

// [[Rcpp::export]]
double logSumExp_C(NumericVector lx){
  double a = max(lx);
  double sum_term = 0.0;
  int n = lx.size();
  
  for(int i = 0; i < n; i++){
    sum_term += exp(lx[i] - a);
  }
  return a + log(sum_term);
}

// [[Rcpp::export]]
arma::mat tensor_rotate_X(arma::mat X,
                          arma::mat eigvecS,
                          arma::mat eigvecT,
                          arma::vec scales){
  
  int p = X.n_cols;
  int n = X.n_rows;
  int ns = eigvecS.n_cols;
  int nt = eigvecT.n_cols;
  arma::mat out(arma::size(X));
  
  for(int j = 0; j < p; j++){
    arma::mat X_j = arma::reshape(X.col(j), nt, ns);
    out.col(j) = arma::vectorise(eigvecT.t()*X_j*eigvecS);
  }
  
  for(int i = 0; i < n; i++){
    out.row(i) /= scales[i];
  }
  
  return(out);
}

// [[Rcpp::export]]
double dmvnorm_kronecker_nugget(arma::colvec x, 
                                arma::colvec mu, 
                                arma::mat Sigma_S, //spatial
                                arma::mat Sigma_T, //temporal
                                double tau2, //variance
                                double sigma2, //nugget
                                bool logd = true){
  
  // x is vector of space time, filled by time first.
  // so that reshape(x, nt, ns) has a time series in each col.
  
  int ns = Sigma_S.n_cols;
  int nt = Sigma_T.n_cols;
  int k = ns*nt; //num of observations
  double out;
  
  //get eigen decomposition for both S T covs
  arma::vec eigval_S; arma::vec eigval_T;
  arma::mat eigvec_S; arma::mat eigvec_T;
  
  arma::eig_sym(eigval_S, eigvec_S, tau2*Sigma_S);
  arma::eig_sym(eigval_T, eigvec_T, Sigma_T);  
  
  //reshape and rotate response
  arma::mat y_mat = arma::reshape(x-mu, nt, ns);
  arma::mat y_new = eigvec_T.t()*y_mat*eigvec_S;
  arma::vec y_vec = arma::vectorise(y_new);
  
  //reshaped has diagonal covariance matrix.
  arma::vec eig_kron = arma::kron(eigval_S, eigval_T);
  arma::vec var_vec = eig_kron + sigma2;
  double quads = sum(pow(y_vec,2)/var_vec);
  
  out = -(k/2.0)*log2pi - 0.5 * quads - sum(log(var_vec))/2;     
  
  if (logd == false) {out = exp(out);}
  return(out);
}

// The BGL-SS Gibbs sampler, for separable ST data.

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List BGL_SS_MCMC(arma::mat C_full,
                 arma::colvec Y,
                 IntegerVector partition, //which covars belong to groups
                 arma::mat DistS, //spatial distances
                 arma::mat DistT, //temporal distances
                 double sig2_alpha = 0.1,
                 double sig2_gamma = 0.1,
                 double sig2_beta = 1e6,
                 double pi0_a = 1,
                 double pi0_b = 1,
                 double tune_rho = 1,
                 int n_samples = 1000,
                 int print_step = 100,
                 String rhomethod = "exact",
                 bool lamsep = false){
  
  Rcout << "Setting up" << std::endl;
  int G = (unique(partition)).size();
  NumericVector m = rep(0.0, G);
  for(int g = 0; g < G; g++){
    m[g] = sum(partition == g) - 1; //no linear term
  }
  
  //initial parameters
  // int pc = C_full.n_cols;
  NumericVector lambda = rep(1e-4, G);
  double sig2 = .3; //the nugget (used in BGL)
  double nu2 = .3; //the scale
  // double nugg = 1; //nugget
  NumericVector tau2 = rep(1.0, G);
  double pi0 = 0.5;
  NumericVector Z = rep(1.0, G);
  NumericVector Zeta = rep(0.0, G); //nonlinear
  arma::vec D_vec = rep(0.0, G);
  double rhoS = 0.1;
  double rhoT = 0.1;
  arma::mat CTC = C_full.t()*C_full; CTC.diag() += 1e-9;
  arma::colvec w = arma::solve(CTC, C_full.t()*Y); //all coefficients
  // arma::colvec u = rep(1.0, pc - G); //nonlinear terms
  // arma::colvec b = rep(1.0, G); //linear terms
  int n = Y.n_rows;
  int nt = DistT.n_cols;
  int ns = DistS.n_cols;
  arma::vec partition_ = as<arma::vec>(partition);
  
  //initial variance-covariance matrix
  arma::mat Sigma_S = exp(-rhoS*DistS);
  arma::mat Sigma_T = exp(-rhoT*DistT);
  
  //samples
  IntegerMatrix Z_samples(G, n_samples); //linear terms
  IntegerMatrix Zeta_samples(G, n_samples); //nonlinear terms
  NumericMatrix tau2_samples(G, n_samples);
  std::fill(tau2_samples.begin(), tau2_samples.end(), 0.0);
  arma::mat w_samples(C_full.n_cols, n_samples);
  NumericVector sig2_samples(n_samples);
  NumericVector nu2_samples(n_samples);
  NumericVector pi0_samples(n_samples);
  NumericVector rhoS_samples(n_samples);
  NumericVector rhoT_samples(n_samples);
  NumericMatrix lambda_samples(G, n_samples);
  std::fill(lambda_samples.begin(), lambda_samples.end(), 0.1);
  NumericVector did_accept(n_samples);
  std::fill(did_accept.begin(), did_accept.end(), 0.0);
  
  double log_pi1 = 0.0;
  arma::mat root_Sig_g;
  arma::mat Sigma_g;
  arma::colvec XT_Resids;
  
  double tmp = 0.1;
  double mn = 0.1;
  arma::vec mn_now = C_full*w;
  
  double lk_curr = dmvnorm_kronecker_nugget(Y, 
                                            mn_now, 
                                            Sigma_S, //spatial
                                            Sigma_T, //temporal
                                            nu2, //variance
                                            sig2); //nugget
  
  //begin sampling
  Rcout << "Begin sampling" << std::endl;
  for(int i = 0; i < n_samples; i++){
    if (i % 10 == 0){
      Rcpp::checkUserInterrupt();    // check for interrupt every 10 iterations
    }
    
    //get eigen decomposition for both S T covs
    arma::vec eigval_S; arma::vec eigval_T;
    arma::mat eigvec_S; arma::mat eigvec_T;
    
    arma::eig_sym(eigval_S, eigvec_S, nu2*Sigma_S);
    arma::eig_sym(eigval_T, eigvec_T, Sigma_T);  
    arma::vec new_vars = arma::kron(eigval_S, eigval_T) + sig2; //new ind. variances
    
    for(int g = 1; g < G; g++){ //skip intercept
      // sample b.g
      arma::mat C_g = C_full.cols(find(partition_ == g));
      arma::mat X_g = C_g.col(0); //keep linear part
      arma::mat Z_g = C_g.cols(1, C_g.n_cols - 1); //remove linear part
      
      arma::mat C_notg = C_full.cols(find(partition_ != g));
      arma::colvec w_notg = w.rows(find(partition_ != g));
      arma::colvec u_g = w.rows(find(partition_ == g));
      u_g.shed_row(0); //remove linear part
      
      //reshape and rotate and scale "residuals"
      arma::mat y_mat = arma::reshape(Y - (C_notg*w_notg) - (Z_g*u_g), nt, ns); //reshape resids
      arma::mat y_new = eigvec_T.t()*y_mat*eigvec_S; //rotate resids
      arma::vec y_vec = arma::vectorise(y_new); //stretch back out
      y_vec = y_vec/sqrt(new_vars); //scale by sds
      
      //reshape design
      arma::mat X_lin_star = tensor_rotate_X(X_g, eigvec_S, eigvec_T, sqrt(new_vars));
      arma::mat ctc_g = X_lin_star.t()*X_lin_star;
      arma::mat I_g; I_g.eye(size(ctc_g));
      arma::mat Lambda = I_g/sig2_beta;
      arma::mat chol_Sig_inv = chol(ctc_g + Lambda);
      root_Sig_g = trans(inv(trimatu(chol_Sig_inv)));
      Sigma_g = root_Sig_g.t()*root_Sig_g;
      
      XT_Resids = X_lin_star.t()*y_vec;
      
      log_pi1 = log(1-pi0) - (log(sig2_beta)/2)+arma::sum(log(root_Sig_g.diag()))+pow(norm(root_Sig_g*XT_Resids),2)/(2);
      NumericVector lx = NumericVector::create(log(pi0),log_pi1);
      double log_lg = log(pi0) - logSumExp_C(lx);
      Z[g] = 1 - rbinom(1,1,exp(log_lg))[0];
      
      if(Z[g] == 0){
        arma::uvec which_ones = find(partition_ == g);
        which_ones.shed_rows(1, which_ones.n_elem - 1); //remove nonlinear terms
        w.rows(which_ones).zeros();
      }else{
        arma::uvec which_ones = find(partition_ == g);
        which_ones.shed_rows(1, which_ones.n_elem - 1); //remove nonlinear terms
        arma::colvec mu_g = Sigma_g*XT_Resids;
        w.rows(which_ones) = (root_Sig_g.t())*(as<arma::colvec>(rnorm(1))) + mu_g;
      }
      
      // sample u.g
      X_lin_star = tensor_rotate_X(Z_g, eigvec_S, eigvec_T, sqrt(new_vars));
      ctc_g = X_lin_star.t()*X_lin_star;
      I_g.eye(size(ctc_g));
      Lambda = I_g/tau2[g];
      double m_g = m[g];
      chol_Sig_inv = chol(ctc_g + Lambda);
      root_Sig_g = trans(inv(trimatu(chol_Sig_inv)));
      Sigma_g = root_Sig_g.t()*root_Sig_g;
      
      arma::colvec b_g = w.rows(find(partition_ == g));
      b_g.shed_rows(1, b_g.n_elem - 1); //keep linear part
      
      //reshape and rotate and scale "residuals"
      y_mat = arma::reshape(Y - (C_notg*w_notg) - (X_g*b_g), nt, ns); //reshape resids
      y_new = eigvec_T.t()*y_mat*eigvec_S; //rotate resids
      y_vec = arma::vectorise(y_new); //stretch back out
      y_vec = y_vec/sqrt(new_vars); //scale by sds
      
      XT_Resids = X_lin_star.t()*y_vec;
      
      log_pi1 = log(1-pi0) - ((m_g)*log(tau2[g])/2)+arma::sum(log(root_Sig_g.diag()))+pow(norm(root_Sig_g*XT_Resids),2)/(2);
      lx = NumericVector::create(log(pi0),log_pi1);
      log_lg = log(pi0) - logSumExp_C(lx);
      Zeta[g] = 1 - rbinom(1,1,exp(log_lg))[0];
      
      if(Zeta[g] == 0){
        arma::uvec which_ones = find(partition_ == g);
        which_ones.shed_row(0); //remove linear term
        w.rows(which_ones).zeros();
      }else{
        arma::uvec which_ones = find(partition_ == g);
        which_ones.shed_row(0); //remove linear term
        arma::colvec mu_g = Sigma_g*XT_Resids;
        w.rows(which_ones) = (root_Sig_g.t())*(as<arma::colvec>(rnorm(m_g))) + mu_g;
      }
      
      //sample tau2
      if(Zeta[g] == 0){
        double shp_tau2 = (m_g + 1)/2;
        double rt_tau2 = pow(lambda[g],2)/2;
        tau2[g] = R::rgamma(shp_tau2, 1/rt_tau2);
      } else {
        arma::vec u_vec = w.rows(find(partition_ == g));
        arma::uword ind1 = 1;
        arma::uword ind2 = u_vec.n_elem - 1;
        arma::vec u_vec2 = u_vec.subvec(ind1, ind2); //not the linear term
        tmp = norm(u_vec2);
        mn = lambda[g]/tmp;
        tau2[g] = abs(1/rinvgauss_rcpp(1,mn,pow(lambda[g],2)))[0];
      }
      
      if(tau2[g] > 1e6){tau2[g] = 1e6;} //a check, for some wild tau2s
      
      //need D_vec
      arma::vec w_gtemp = w.rows(find(partition_ == g));
      arma::uword ind_1 = 1;
      arma::uword ind_2 = w_gtemp.n_elem - 1;
      arma::vec u_gtemp = w_gtemp.subvec(ind_1, ind_2);
      D_vec[g] = pow(norm(u_gtemp),2)/tau2[g];
    }
    
    // sample intercept(s)
    arma::mat C_g = C_full.cols(find(partition_ == 0));
    arma::mat X_lin_star = tensor_rotate_X(C_g, eigvec_S, eigvec_T, sqrt(new_vars));
    arma::mat ctc_g = X_lin_star.t()*X_lin_star;
    arma::mat IDint(ctc_g.n_rows, ctc_g.n_cols); IDint.eye();
    arma::mat chol_Sig_inv = chol(ctc_g + IDint/sig2_beta);
    root_Sig_g = trans(inv(trimatu(chol_Sig_inv)));
    Sigma_g = root_Sig_g.t()*root_Sig_g;
    
    arma::mat C_notg = C_full.cols(find(partition_ != 0));
    arma::colvec w_notg = w.rows(find(partition_ != 0));
    
    //reshape and rotate and scale "residuals"
    arma::mat y_mat = arma::reshape(Y - (C_notg*w_notg), nt, ns); //reshape resids
    arma::mat y_new = eigvec_T.t()*y_mat*eigvec_S; //rotate resids
    arma::vec y_vec = arma::vectorise(y_new); //stretch back out
    y_vec = y_vec/sqrt(new_vars); //scale by sds
    
    XT_Resids = X_lin_star.t()*y_vec;
    
    arma::colvec mu_g = Sigma_g*XT_Resids;
    w.rows(find(partition_ == 0)) = (root_Sig_g.t())*(as<arma::colvec>(rnorm(C_g.n_cols))) + mu_g;
    
    mn_now = C_full*w;
    
    // // sample sig2
    // y_mat = arma::reshape(Y - mn_now, nt, ns); //reshape resids
    // y_new = eigvec_T.t()*y_mat*eigvec_S; //rotate resids
    // y_vec = arma::vectorise(y_new); //stretch back out
    // arma::vec sumterms = Y - mn_now - y_vec/new_vars;
    // double shp = n/2;// + (m*Zeta)[0]/2;
    // double rt = dot(sumterms, sumterms)/2;// + sum(D_vec)/2;
    // sig2 = 1/(rgamma(1, shp, 1/rt))[0];
    
    double qcurr = lk_curr + log(sig2);// + R::dgamma(rho, 1, 1/2, 1);
    
    double sig2_prop = exp(log(sig2) + tune_rho*R::rnorm(0,1));
    
    double lk_prop = dmvnorm_kronecker_nugget(Y, 
                                              mn_now, 
                                              Sigma_S, //spatial
                                              Sigma_T, //temporal
                                              nu2, //variance
                                              sig2_prop);
    
    double qprop = lk_prop + log(sig2_prop);// + R::dgamma(rho_prop, 1, 1/2, 1);
    
    if(log(runif(1, 0, 1)[0]) < qprop - qcurr){
      sig2 = sig2_prop;
      lk_curr = lk_prop;
    }
    
    // sample pi0
    // pi0 = rbeta(1, 
    //             pi0_a + sum(1-Z) - (1 - Z[0]) + sum(1-Zeta) - (1 - Zeta[0]), 
    //             pi0_b + 2*(G-1) - sum(1-Z) + (1 - Z[0]) - sum(1-Zeta) + (1 - Zeta[0]))[0];
    
    // pi0 = rbeta(1, 
    //             pi0_a + sum(1-Z) + sum(1-Zeta), 
    //             pi0_b + 2*(G) - sum(1-Z) - sum(1-Zeta))[0];
    
    pi0 = rbeta(1, 
                pi0_a + sum(1-Z) + sum(1-Zeta), 
                pi0_b + (G + sum(m)) - sum(1-Z) - sum(1-Zeta))[0];
    
    // sample rho
    qcurr = lk_curr + log(rhoS) + log(rhoT);// + R::dgamma(rho, 1, 1/2, 1);
    
    double rhoS_prop = exp(log(rhoS) + tune_rho*R::rnorm(0,1));
    double rhoT_prop = exp(log(rhoT) + tune_rho*R::rnorm(0,1));
    arma::mat Sigma_S_prop = exp(-rhoS_prop*DistS);
    arma::mat Sigma_T_prop = exp(-rhoT_prop*DistT);
    
    lk_prop =  dmvnorm_kronecker_nugget(Y, 
                                        mn_now, 
                                        Sigma_S_prop, //spatial
                                        Sigma_T_prop, //temporal
                                        nu2, //variance
                                        sig2);
    qprop = lk_prop + log(rhoS_prop) + log(rhoT_prop);// + R::dgamma(rho_prop, 1, 1/2, 1);
    
    if(log(runif(1, 0, 1)[0]) < qprop - qcurr){
      rhoS = rhoS_prop;
      rhoT = rhoT_prop;
      Sigma_S = Sigma_S_prop;
      Sigma_T = Sigma_T_prop;
      lk_curr = lk_prop;
    }
    
    // sample nu
    qcurr = lk_curr + log(nu2);// + R::dgamma(rho, 1, 1/2, 1);
    
    double nu2_prop = exp(log(nu2) + tune_rho*R::rnorm(0,1));
    
    lk_prop =  dmvnorm_kronecker_nugget(Y, 
                                        mn_now, 
                                        Sigma_S, //spatial
                                        Sigma_T, //temporal
                                        nu2_prop, //variance
                                        sig2);
    qprop = lk_prop + log(nu2_prop);// + R::dgamma(rho_prop, 1, 1/2, 1);
    
    if(log(runif(1, 0, 1)[0]) < qprop - qcurr){
      nu2 = nu2_prop;
      lk_curr = lk_prop;
    }
    
    // estimate lambda
    if(lamsep){ //separate lambdas
      if(rhomethod == "empirical"){
        for(int g = 0; g < G; g++){
          if(i > 10){
            NumericVector currtau2 = tau2_samples.row(g);
            lambda[g] = sqrt(m[g]/mean(currtau2[Range(i-10, i)]));
          }
        }
      } else if (rhomethod == "exact"){
        for(int g = 0; g < G; g++){
          lambda[g] = R::rgamma(m[g]/2 + 1 + .1, 1/(tau2[g]/2 + .1));
        }
      }
    } else {
      if(rhomethod == "empirical"){
        lambda = rep(sqrt((sum(m) + G - 1)/sum(tau2)), G); //no intercept. G doesn't include partition == 0
      } else if (rhomethod == "exact"){
        lambda = rep(sqrt(R::rgamma(sum(m)/2 + 1 + .1, 1/(sum(tau2)/2 + .1))), G);
      }
    }
    
    //
    // if(i > 500){
    //   for(int g = 0; g < G; g++){
    //     double where_to_start = i - 100;
    //     Range how_many(where_to_start, i);
    //     NumericVector tau2gs = tau2_samples.row(g);
    //     double mn_tau2g = mean(tau2gs[how_many]);
    //     lambda[g] = pow(m[g]/mn_tau2g, 0.5);
    //   }
    // }
    
    // collect samples
    Z_samples(_,i) = Z;
    Zeta_samples(_,i) = Zeta;
    w_samples.col(i) = w;
    pi0_samples[i] = pi0;
    sig2_samples[i] = sig2;
    tau2_samples(_,i) = tau2;
    lambda_samples(_,i) = lambda;
    rhoS_samples[i] = rhoS;
    rhoT_samples[i] = rhoT;
    nu2_samples[i] = nu2;
    
    if((i%print_step == 0) & (i != 0)){
      Rcpp::Rcout << "iteration: " << i << "\t lambda: " << lambda[0] << "\t sig2: " << sig2 << std::endl;
    }
  }
  
  List out; //returns Z, w, pi0, sig2, tau2, rho, lam samples
  out["Z.samples"] = Z_samples;
  out["Zeta.samples"] = Zeta_samples;
  out["tau2.samples"] = tau2_samples;
  out["w.samples"] = w_samples;
  out["sig2.samples"] = sig2_samples;
  out["pi0.samples"] = pi0_samples;
  out["lam.samples"] = lambda_samples;
  out["rhoS.samples"] = rhoS_samples;
  out["rhoT.samples"] = rhoT_samples;
  out["nu2.samples"] = nu2_samples;
  return out;
}


