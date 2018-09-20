#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
using namespace std;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec Phi2Theta(arma::vec Phi, arma::vec binf, arma::vec bsup)
{
  return (exp(Phi)%(bsup-binf)+binf);
}

// [[Rcpp::export]]
arma::vec Theta2Phi(arma::vec Theta, arma::vec binf, arma::vec bsup)
{
  return log((Theta-binf)/(bsup-binf));
}


//' C++ implementation of the algorithm for parameter calibration
//'
//' Run a Metropolis Hastings within Gibbs algorithm and a Metropolis Hastings algorithm with the covariance matrix estimated on the
//' the sample set generated in the Metropolis within Gibbs. This algorithm is suitable only for models without discrepancy.
//'
//' @param Ngibbs the number of iteration in the Metropolis within Gibbs
//' @param Nmh the number of iteration in the Metropolis Hastings
//' @param theta_init the starting point
//' @param r regulation percentage in the modification of the k in the Metropolis Hastings
//' @param Yf the vector of recorded data
//' @param binf the lower bound of the parameters to calibrate
//' @param bsup the upper bound of the parameters to calibrate
//' @param LogTest the log posterior density distribution
//' @param stream (default=1) if stream=0 the progress bar is disabled
//' @return list of outputs: \itemize{
//' \item PHIwg the points of the Metropolis within Gibbs algorithm in the transformed space
//' \item PHI the points of the Metropolis Hastings algorithm in the transformed space
//' \item THETAwg the points of the Metropolis within Gibbs algorithm in the real space
//' \item THETA the points of the Metropolis Hastings algorithm in the real space
//' \item AcceptationRatioWg the vector of the acceptance ratio for each parameter in the Metropolis within Gibbs
//' \item AcceptationRatio the acceptance ratio in the Metropolis Hastings
//' \item S the covariance computed after the Metropolis within Gibbs
//' \item LikeliWG the likelihood computed at each iteration of the Metropolis within Gibbs algorithm
//' \item Likeli the likelihood computed at each iteration of the Metropolis Hastings algorithm
//'  }
//' @export
// [[Rcpp::export]]
List MetropolisHastingsCpp(int Ngibbs, int Nmh, arma::vec theta_init, arma::vec r, arma::vec Yf,
                           arma::vec binf, arma::vec bsup, Function LogTest, bool disc, bool stream)
{
  // Dimention definition
  double Dim = theta_init.size();
  // Variables declaration
  vec k = 1e-2*ones<vec>(Dim);
  mat PHIwg=randu<mat>(Ngibbs,Dim), THETAwg=randu<mat>(Ngibbs,Dim);
  mat LikeliWG=randu<mat>(Ngibbs,Dim);
  vec Likeli=zeros(Nmh,1);
  vec theta, thetaD;
  double Verr, alpha, beta, alpha2, beta2;
  if (Nmh==0) {Nmh=10;}
  mat PHI= randu<mat>(Nmh,Dim), THETA=randu<mat>(Nmh,Dim);
  // Set the first row of THETA at the initial value (for the MWG)
  THETA.row(0)=theta_init.t();
  // Space changing of THETA in PHI
  PHI.row(0)= Theta2Phi(theta_init, binf, bsup).t();
  // Declaration of the acceptation ratios
  vec AcceptationRatioWg=zeros(Dim,1);
  double AcceptationRatio=0;
  // Functions from R
  Function DefPos("DefPos");
  // Set the first row of THETA at the initial value (for the MH)
  THETAwg.row(0)=THETA.row(0);
  // Space changing of THETA in PHI (for the MH)
  PHIwg.row(0) = PHI.row(0);
  // Defining Theta and measurement variance error
  if (disc == 0)
  {
    theta=theta_init.subvec(0,Dim-2);
    Verr=theta_init(Dim-1);
    // Compute the first ratio alpha
    alpha = as<double>(LogTest(theta,Verr));
  } else
  {
    theta=theta_init.subvec(0,Dim-4);
    thetaD=theta_init.subvec(Dim-3,Dim-2);
    Verr=theta_init(Dim-1);
    alpha = as<double>(LogTest(theta,thetaD,Verr));
  }
  // Display activated if stream = true
  // if (stream==1)
  // {
  //   int consoleWidth = 50;
  //   Rcout << setw(consoleWidth / 2) << " " << " ----------------------------------- " << endl;
  //   Rcout << setw(consoleWidth / 2) << " " << "| Metropolis within Gibbs algorithm |" << endl;
  //   Rcout << setw(consoleWidth / 2) << " " << " ----------------------------------- " << endl;
  //   Rcout << " ----------------------- " << endl;
  //   Rcout << "| Number of iterations: |"<< Ngibbs << endl;
  //   Rcout << " ----------------------- "<< endl;
  //   Rcout << "| Initial Theta:        |" << theta_init.subvec(0,Dim-2).t() << endl;
  //   Rcout << " ----------------------- "<< endl;
  //   Rcout << "| Initial Variance:     |" << theta_init(Dim-1) << endl;
  //   Rcout << " ----------------------- " << endl;
  // }
  for (int i=0; i < (Ngibbs-1); i++)
  {
    /*
    if (stream==1)
    {
      // beggining of the bar progress
      Rcout.flush();
      if ((i+2)%(Ngibbs/barWidth)==0)
      {
        float progress = (float)(i+2)/(float)(Ngibbs);
        Rcout << "[";
        q = ((float)(i+2)*(float)barWidth)/(float)(Ngibbs);
        if (q<barWidth) Rcout << string(q, '=');
        else if (q==barWidth) Rcout << string(q, '>');
        else Rcout << " ";
        Rcout << "] " << int(progress * 100.0) << " %\r";
      }
      // end bar progress
    }
    */
    // Get the ith point
    vec phi_star = PHIwg.row(i).t();
    // Beggining of the MHWG part
    for (int j=0; j<Dim; j++)
    {
      if (j>0){
        phi_star.subvec(0,j) = PHIwg.row(i).subvec(0,j).t();
      }
      // Proposition of a new point in the Log-normalized space
      phi_star(j) = randn<double>()*k(j)+PHIwg(i,j);
      vec theta_star = Phi2Theta(phi_star,binf,bsup);
      if (disc == 0)
      {
        theta = theta_star.subvec(0,Dim-2);
        Verr = theta_star(Dim-1);
        beta = as<double>(LogTest(theta,Verr));
      } else
      {
        theta=theta_star.subvec(0,Dim-4);
        thetaD=theta_star.subvec(Dim-3,Dim-2);
        Verr=theta_star(Dim-1);
        beta = as<double>(LogTest(theta,thetaD,Verr));
      }
      // Computing the new LogPost for theta star
      // Ratio for the MH
      //double beta = -2500;
      double logR = beta-alpha;
      if (log(randu<double>()) < logR)
      {
        // Acceptation case
        PHIwg(i+1,j)=phi_star(j);
        THETAwg(i+1,j)=theta_star(j);
        alpha = beta;
        LikeliWG(i+1,j)=beta;
        AcceptationRatioWg(j) += 1;
      }
      else
      {
        // Rejection case
        PHIwg(i+1,j)=PHIwg(i,j);
        THETAwg(i+1,j)=THETAwg(i,j);
        LikeliWG(i+1,j)=alpha;
      }
      // Adaptive algorithm
      if (i%10 == 0)
      {
        if (AcceptationRatioWg(j)/(i+1)<0.2)
        {
          k(j) = k(j)*(1-r(0));
        }
        if (AcceptationRatioWg(j)/(i+1)>0.5)
        {
          k(j) = k(j)*(1+r(0));
        }
      }
    }
    if (stream==1)
    {
      if (i%100 == 0)
      {
        int consoleWidth = 50;
        Rcout << setw(consoleWidth/2)  << "       Metropolis within Gibbs       " << endl;
        Rcout << "At iterate  " << i  << endl;
        Rcout << "Theta = " << THETAwg.row(i).subvec(0,Dim-2) <<endl;
        Rcout << "Variance = " << THETAwg.row(i)(Dim-1) << endl;
        Rcout << "Current accepation Rate: "<< AcceptationRatioWg.t()/(i+1) *100 << endl;
        Rcout << "Ratio: " << r(0) <<"%"<<endl;
        Rcout << "Current k: "<< k.t() << endl;
        Rcout << setw(consoleWidth/2)  << " ----------------------------------- " << endl;
      }
    }
}
  // Establishment of the new covariance matrix
  mat S = as<mat>(DefPos(cov(PHIwg.rows(0.5*Ngibbs,(Ngibbs-1)))));
  vec NewPhi = mean(PHIwg.rows(0.5*Ngibbs,(Ngibbs-1)),0).t();
  vec NewTheta = Phi2Theta(NewPhi,binf,bsup);
  THETA.row(0) = NewTheta.t();
  PHI.row(0) = NewPhi.t();
  //vec NewTheta = THETA.row(0).t();
  //vec NewPhi = PHI.row(0).t();
  // Setting a new starting point for the MH algorithm
  // if (stream==1)
  // {
  //   // Rcout << "Covariance matrix estimated" <<endl;
  //   // Rcout << endl;
  //   // int consoleWidth = 50;
  //   // Rcout << setw(consoleWidth / 2) << " " << " ---------------------- " << endl;
  //   // Rcout << setw(consoleWidth / 2) << " " << "| Metropolis algorithm |" << endl;
  //   // Rcout << setw(consoleWidth / 2) << " " << " ---------------------- " << endl;
  //   // Rcout << " ----------------------- " << endl;
  //   // Rcout << "| Number of iterations: |"<< Nmh << endl;
  //   // Rcout << " ----------------------- "<< endl;
  //   // Rcout << "| Initial Theta:        |" << NewTheta.subvec(0,Dim-2).t() << endl;
  //   // Rcout << " ----------------------- "<< endl;
  //   // Rcout << "| Initial Variance:     |" << NewTheta(Dim-1) << endl;
  //   // Rcout << " ----------------------- "<<  endl;
  // }
  // t is the new k for the second part of the algp
  double t=1;
  // Store that ratio in alpha2
  if (disc == 0)
  {
    theta = NewTheta.subvec(0,Dim-2);
    Verr = NewTheta(Dim-1);
    alpha2 = as<double>(LogTest(theta,Verr));
  } else
  {
    theta = NewTheta.subvec(0,Dim-4);
    thetaD=NewTheta.subvec(Dim-3,Dim-2);
    Verr = NewTheta(Dim-1);
    alpha2 = as<double>(LogTest(theta,thetaD,Verr));
  }
  for (int i=0; i<(Nmh-1); i++)
  {
    /*
    if (stream==1)
    {
      // beggining of the bar progress
      Rcout.flush();
      if ((i+2)%(Nmh/barWidth)==0)
      {
        float progress = (float)(i+2)/(float)(Nmh);
        Rcout << "[";
        q = ((float)(i+2)*(float)barWidth)/(float)(Nmh);
        if (q<barWidth) Rcout << string(q, '=');
        else if (q==barWidth) Rcout << string(q, '>');
        else Rcout << " ";
        Rcout << "] " << int(progress * 100.0) << " %\r";
      }
      // end bar progress
    }
    */
    // The new point is found from the new phi computed line 153
    vec phi_star = mvnrnd(PHI.row(i).t(), t*S);
    // In the original space
    vec theta_star = Phi2Theta(phi_star,binf,bsup);
    if (disc == 0)
    {
      theta = theta_star.subvec(0,Dim-2);
      Verr = theta_star(Dim-1);
      beta2 = as<double>(LogTest(theta,Verr));
    } else
    {
      theta = theta_star.subvec(0,Dim-4);
      thetaD = theta_star.subvec(Dim-2,Dim-4);
      Verr = theta_star(Dim-1);
      beta2 = as<double>(LogTest(theta,thetaD,Verr));
    }
    // Computing the logPost and Ratio for the overall new point
    double logR2 = beta2 - alpha2;
    if(log(randu<double>()) < logR2)
    {
      // Acceptation of the new point
      PHI.row(i+1)=phi_star.t();
      THETA.row(i+1)=theta_star.t();
      Likeli(i+1)=beta2;
      alpha2 = beta2;
      AcceptationRatio += 1;
    }
    else
    {
      // Rejection
      Likeli(i+1)=alpha2;
      PHI.row(i+1)=PHI.row(i);
      THETA.row(i+1)=THETA.row(i);
    }
    // Adaptive part on t
    if (i%10==0)
    {
      if (AcceptationRatio/(i+1)<0.2)
      {
        t = t*(1-r(1));
      }
      if (AcceptationRatio/(i+1)>0.5)
      {
        t = t*(1+r(1));
      }
    }
    if (stream==1)
    {
      if (i%100 == 0)
      {
        int consoleWidth = 50;
        Rcout << setw(consoleWidth/2)  << "         Metropolis Algorithm        " << endl;
        Rcout << "At iterate  " << i  << endl;
        Rcout << "Theta = " << THETA.row(i).subvec(0,Dim-2) <<endl;
        Rcout << "Variance = " << THETA.row(i)(Dim-1) << endl;
        Rcout << "Current accepation Rate: "<< AcceptationRatio/(i+1)*100 << endl;
        Rcout << "Ratio: " << r(1) <<"%"<<endl;
        Rcout << "Current k: "<< t << endl;
        Rcout << setw(consoleWidth/2)  << " ----------------------------------- " << endl;
      }
    }
  }
  return List::create(Named("PHIwg")=PHIwg,Named("THETAwg")=THETAwg,Named("PHI")=PHI,Named("THETA")=THETA,
                            Named("AcceptationRatio")=AcceptationRatio, Named("AcceptationRatioWg")=AcceptationRatioWg
                        , Named("S")=S,Named("LikeliWG")=LikeliWG, Named("Likeli")=Likeli);
}

