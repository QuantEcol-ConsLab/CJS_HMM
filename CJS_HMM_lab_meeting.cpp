#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //Data
  DATA_IMATRIX(CH);                //capture histories
  int n= CH.rows();                // # animals
  int n_occasions= CH.cols();      // ## occasions
  DATA_IVECTOR(first);             // time at first capture
  
  
  //Parameters
  PARAMETER(logit_phi);            //logit mean survival
  Type phi = invlogit(logit_phi);  // mean survival
  PARAMETER(logit_p);              // logit p
  Type p = invlogit(logit_p);      // p
  PARAMETER_VECTOR(epsilon);       // random time effects
  PARAMETER(log_sigma);            // log SD of time effects
  Type sigma = exp(log_sigma);     // SD of time effects
  
  //Variables
  matrix<Type> delta(1,2);         // initial disttribution
  delta(0,0)=1;                    // indexing starts at 0 not 1 in TMB
  delta(0,1)=0;
  
  matrix<Type> Gamma(2,2);         // transition matrix
  Gamma(1,0) = 0;
  Gamma(1,1) = 1;
  
  matrix<Type> obs(2,2);           // observatin matrix
  obs(0,0) = (1-p);
  obs(0,1) = p;
  obs(1,0) = 1;
  obs(1,1) = 0;
  
  //Likelihood
  
  // random effects likelihood to be integrated out using Laplace approx
  vector<Type> phi_t(n_occasions-1);                   // declare vector to hold annual survival probs on prop scale
  Type like_epsilon=0;                                 // initialize random effects LL
  for (int j = 0; j<(n_occasions-1); j++){             // loop over time
    phi_t(j)=invlogit(logit_phi+epsilon(j));           // calculate survival at time on probability scale
    like_epsilon+=dnorm(epsilon(j),Type(0),sigma,true);// probability of random effects
  }
  
  Type lscale = 0;                                     // initialize data LL
  matrix<Type> alpha(1,2);                             // declare 1x2 row matrix to hold probabilities within the forward algorithem

  
  for (int i=0; i<n ; i++){                            // loop through animals 
    alpha = delta;                                     // start probs at inital distribution
   for (int j = first(i); j<(n_occasions); j++){       // loop through time
      
      Gamma(0,0) = phi_t(j-1);                         // fill first row of transition matrix
      Gamma(0,1) = 1-phi_t(j-1);
      
      alpha = alpha * Gamma ;                          // matrix multiply by transition matrix
      alpha(0,0)*=obs(0,CH(i,j));                      // elementwise multiply by observation probabilitie
      alpha(0,1)*=obs(1,CH(i,j));                      // I was having trouble with extracting a column from the observatio nmatrix so Ijust did this workaround
      lscale+= log(sum(alpha));                        // accumulate log likelihood
     
      alpha = alpha / sum( alpha);                     // scale
   }
 }
  

  //Reporting
  REPORT(n);
  REPORT(n_occasions);
  REPORT(alpha);
  REPORT(obs);
  REPORT(Gamma);
  REPORT(delta);
  REPORT(lscale);
  ADREPORT(p);                                         // ask for delta method standard deviations of variables inthe "ADreport"
  ADREPORT(phi);
  ADREPORT(phi_t);
  ADREPORT(sigma);
  Type NLL = -like_epsilon-lscale;                     // calculate NLL
 return(NLL);                                          
}