#include <RcppArmadillo.h>
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

// The BGL-SS Gibbs sampler, for AR1 data.

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List M_BGL_SS_AR1_MCMC(arma::mat C_full,
                       arma::colvec Y,
                       IntegerVector partition,
                       arma::mat Dist,
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
  
  int G = (unique(partition)).size();
  NumericVector m = rep(0.0, G);
  for(int g = 0; g < G; g++){
    m[g] = sum(partition == g) - 1; //no linear term
  }
  
  //initial parameters
  int pc = C_full.n_cols;
  NumericVector lambda = rep(1e-4, G);
  double sig2 = 1;
  double nugg = 1; //nugget
  NumericVector tau2 = rep(1.0, G);
  double pi0 = 0.5;
  NumericVector Z = rep(1.0, G);
  NumericVector Zeta = rep(0.0, G); //nonlinear
  arma::vec D_vec = rep(0.0, G);
  double rho = 0.03;
  arma::colvec w = rep(1.0, pc); //all coefficients
  arma::colvec u = rep(1.0, pc - G); //nonlinear terms
  arma::colvec b = rep(1.0, G); //linear terms
  int n = Y.n_rows;
  double sig2beta = sig2_beta; //weakly informative normal prior on linear terms
  arma::mat ID = arma::mat(size(Dist)); ID.eye();
  
  //precision matrix for data
  // arma::mat CorMat = exp(-pow(Dist/rho, 2)) + 1e-9*ID;
  arma::mat CorMat = exp(-abs(Dist/rho));
  arma::mat Hinv = arma::inv_sympd(CorMat);
  
  arma::vec partition_ = as<arma::vec>(partition);
  
  //samples
  IntegerMatrix Z_samples(G, n_samples); //linear terms
  IntegerMatrix Zeta_samples(G, n_samples); //nonlinear terms
  NumericMatrix tau2_samples(G, n_samples);
  std::fill(tau2_samples.begin(), tau2_samples.end(), 0.0);
  arma::mat w_samples(C_full.n_cols, n_samples);
  NumericVector sig2_samples(n_samples);
  NumericVector pi0_samples(n_samples);
  NumericVector rho_samples(n_samples);
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
  
  //begin sampling
  for(int i = 0; i < n_samples; i++){
    if (i % 10 == 0){
      Rcpp::checkUserInterrupt();    // check for interrupt every 10 iterations
    }
    
    for(int g = 1; g < G; g++){ //skip intercept
      // sample b.g
      arma::mat C_g = C_full.cols(find(partition_ == g));
      arma::mat X_g = C_g.col(0); //keep linear part
      arma::mat Z_g = C_g.cols(1, C_g.n_cols - 1); //remove linear part
      
      arma::mat ctc_g = trans(X_g)*Hinv*X_g/sig2;
      arma::mat I_g; I_g.eye(size(ctc_g));
      arma::mat Lambda = I_g/sig2_beta;
      arma::mat chol_Sig_inv = chol(ctc_g + Lambda);
      root_Sig_g = trans(inv(trimatu(chol_Sig_inv)));
      Sigma_g = root_Sig_g.t()*root_Sig_g;
      
      arma::mat C_notg = C_full.cols(find(partition_ != g));
      arma::colvec w_notg = w.rows(find(partition_ != g));
      arma::colvec u_g = w.rows(find(partition_ == g));
      u_g.shed_row(0); //remove linear part

      XT_Resids = trans(X_g)*Hinv*(Y - (C_notg*w_notg) - (Z_g*u_g))/sig2;
      
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
      ctc_g = trans(Z_g)*Hinv*Z_g;
      I_g.eye(size(ctc_g));
      Lambda = I_g/tau2[g];
      double m_g = m[g];
      chol_Sig_inv = chol(ctc_g + Lambda);
      root_Sig_g = trans(inv(trimatu(chol_Sig_inv)));
      Sigma_g = root_Sig_g.t()*root_Sig_g;
      
      arma::colvec b_g = w.rows(find(partition_ == g));
      b_g.shed_rows(1, b_g.n_elem - 1); //keep linear part
      XT_Resids = trans(Z_g)*Hinv*(Y - (C_notg*w_notg) - (X_g*b_g));

      log_pi1 = log(1-pi0) - ((m_g)*log(tau2[g])/2)+arma::sum(log(root_Sig_g.diag()))+pow(norm(root_Sig_g*XT_Resids),2)/(2*sig2);
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
        w.rows(which_ones) = (sqrt(sig2)*root_Sig_g.t())*(as<arma::colvec>(rnorm(m_g))) + mu_g;
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
        mn = lambda[g]*sqrt(sig2)/tmp;
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

    // sample intercept
    arma::mat C_g = C_full.cols(find(partition_ == 0));
    arma::mat ctc_g = trans(C_g)*Hinv*C_g;
    double m_g = 1;
    arma::mat chol_Sig_inv = chol(ctc_g + 1/sig2_beta);
    root_Sig_g = trans(inv(trimatu(chol_Sig_inv)));
    Sigma_g = root_Sig_g.t()*root_Sig_g;
    
    arma::mat C_notg = C_full.cols(find(partition_ != 0));
    arma::colvec w_notg = w.rows(find(partition_ != 0));
    XT_Resids = trans(C_g)*Hinv*(Y - (C_notg*w_notg));
    
    arma::colvec mu_g = Sigma_g*XT_Resids;
    w.rows(find(partition_ == 0)) = (sqrt(sig2)*root_Sig_g.t())*(as<arma::colvec>(rnorm(m_g))) + mu_g;
    
    // sample sig2
    arma::vec sumterms = (Y-C_full*w);
    // double shp = n/2 + (m*Zeta)[0]/2 + sig2_alpha;
    // double rt = dot(sumterms, sumterms)/2 + sum(D_vec)/2 + sig2_gamma;
    double shp = n/2 + (m*Zeta)[0]/2;
    double rt = dot(sumterms, sumterms)/2 + sum(D_vec)/2;
    sig2 = 1/(rgamma(1, shp, 1/rt))[0];
    
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
    double lk_curr = lk_Chol_arma(Y, C_full*w, sig2*CorMat) + log(rho);// + R::dgamma(rho, 1, 1/2, 1);
    
    double rho_prop = exp(log(rho) + tune_rho*R::rnorm(0,1));
    // arma::mat CorMat_prop = exp(-pow(Dist/rho_prop, 2)) + 1e-9*ID;
    arma::mat CorMat_prop = exp(-abs(Dist/rho_prop));
    
    double lk_prop = lk_Chol_arma(Y, C_full*w, sig2*CorMat_prop) + log(rho_prop);// + R::dgamma(rho_prop, 1, 1/2, 1);
    // Rcout << lk_prop << "\t" << lk_curr << std::endl;
    if(log(runif(1, 0, 1)[0]) < lk_prop - lk_curr){
      rho = rho_prop;
      CorMat = CorMat_prop;
      Hinv = arma::inv_sympd(CorMat);
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
    rho_samples[i] = rho;
    
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
  out["rho.samples"] = rho_samples;
  return out;
}


