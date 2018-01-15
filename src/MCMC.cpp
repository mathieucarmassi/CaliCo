#include <RcppArmadillo.h>
// [[Rcpp::depends(Matrix,RcppArmadillo)]]
#include <iostream>
#include <armadillo>
#include <math.h>


using namespace Rcpp;
using namespace std;
using namespace arma;


// [[Rcpp::export]]
List MetropolisHastingsCpp(int Ngibbs, int Nmh, arma::vec theta_init, arma::vec k, arma::mat SIGMA, arma::vec Yf,
                           arma::vec binf, arma::vec bsup, Function LogTest, int stream)
{
  double Dim = theta_init.size();
  int D;
  arma::mat PHIwg=randu<arma::mat>(Ngibbs,Dim), THETAwg=randu<arma::mat>(Ngibbs,Dim);
  if (Nmh!=0) {D=Nmh;} else {D=10;}
  arma::mat PHI= randu<arma::mat>(D,Dim), THETA=randu<arma::mat>(D,Dim);
  THETA.row(0)=theta_init.t();
  PHI.row(0)= log((THETA.row(0).t()-binf)/(bsup-binf)).t();
  double AcceptationRatio=0;
  arma::vec AcceptationRatioWg=zeros(Dim,1);
  Function unscale("unscale"), rnorm("rnorm"), mvrnorm("mvrnorm"), runif("runif"), DefPos("DefPos");
  THETAwg.row(0)=theta_init.t();
  PHIwg.row(0) = log((THETAwg.row(0).t()-binf)/(bsup-binf)).t();
  arma::vec theta=theta_init.rows(0,Dim-2);
  double Verr=THETAwg(0,(Dim-1));
  //Rcpp::List res = as<Rcpp::List>(model(theta,Verr));
  //arma::vec Yg=res["y"];
  // arma::vec Yg=as<arma::vec>(model(theta,Verr));
  double alpha = as<double>(LogTest(theta,Verr));
  double alpha2 = alpha;
  if (stream==1)
  {
    cout << "Begin of the Metropolis within Gibbs algorithm" << endl;
    cout << "Number of iterations "<< Ngibbs << endl;
  }
  int barWidth = 40;
  int q = 0;
  for (int i=0; i<(Ngibbs-1); i++)
  {
    if (stream==1)
    {
      // beggining of the bar progress
      cout.flush();
      if ((i+2)%(Ngibbs/barWidth)==0)
      {
        float progress = (float)(i+2)/(float)(Ngibbs);
        cout << "[";
        q = ((float)(i+2)*(float)barWidth)/(float)(Ngibbs);
        if (q<barWidth) cout << string(q, '=');
        else if (q==barWidth) cout << string(q, '>');
        else cout << " ";
        cout << "] " << int(progress * 100.0) << " %\r";
      }
      // end bar progress
    }
    vec phi_star = PHIwg.row(i).t();
    vec theta_star = THETAwg.row(i).t();
    for (int j=0; j<Dim; j++)
    {
      if (j>0){
        phi_star.rows(0,j) = PHIwg.row(i+1).cols(0,j).t();
      }
      phi_star(j) = as<double>(rnorm(1,PHIwg(i,j),sqrt(k(j)*SIGMA(j,j))));
      theta_star(j) = as<double>(unscale(exp(phi_star(j)),binf(j),bsup(j)));
      if (j == Dim-1)
      {
        while (theta_star(j) < 0)
        {
          phi_star(j) = as<double>(rnorm(1,PHIwg(i,j),sqrt(k(j)*SIGMA(j,j))));
          theta_star(j) = as<double>(unscale(exp(phi_star(j)),binf(j),bsup(j)));
        }
      }
      Verr = theta_star((Dim-1));
      theta = theta_init.rows(0,Dim-2);
      double beta = as<double>(LogTest(theta,Verr));
      double logR = beta-alpha;
      if (log(as<double>(runif(1,0,1))) < logR)
      {
        PHIwg(i+1,j)=phi_star(j);
        THETAwg(i+1,j)=theta_star(j);
        alpha = beta;
        AcceptationRatioWg(j) += 1;
      }
      else
      {
        PHIwg(i+1,j)=PHIwg(i,j);
        THETAwg(i+1,j)=THETAwg(i,j);
      }
    }
}
  mat Stemp = cov(PHIwg.rows(10/100*Ngibbs,(Ngibbs-1)));
  mat S = as<arma::mat>(DefPos(Stemp));
  mat NewPhi=mean(PHIwg.rows(10/100*Ngibbs,(Ngibbs-1)));
  if (stream==1)
  {
    cout << endl;
    cout << endl;
    cout << "Estimation of the covariance matrix...." <<endl;
    cout << "End of the within gibbs algorithm"<< endl;
    cout << endl;
    cout << "Begin of the metropolis hastings algorithm using the covariance computed" << endl;
    cout << "Number of iterations "<< Nmh <<endl;
  }
  q = 0;
  //double t=1;
  if (Nmh!=0)
  {
  for (int i=0; i<(Nmh-1); i++)
  {
    if (stream==1)
    {
      // beggining of the bar progress
      cout.flush();
      if ((i+2)%(Nmh/barWidth)==0)
      {
        float progress = (float)(i+2)/(float)(Nmh);
        cout << "[";
        q = ((float)(i+2)*(float)barWidth)/(float)(Nmh);
        if (q<barWidth) cout << string(q, '=');
        else if (q==barWidth) cout << string(q, '>');
        else cout << " ";
        cout << "] " << int(progress * 100.0) << " %\r";
      }
      // end bar progress
    }
    vec phi_star = as<vec>(mvrnorm(1,NewPhi.t(),S));
    vec theta_star = as<vec>(unscale(exp(phi_star.t()),binf,bsup));
    theta = theta_init.rows(0,Dim-2);
    Verr = theta_star((Dim-1));
    //Rcpp::List res = as<Rcpp::List>(model(theta,Verr));
    //arma::vec Yg=res["y"];
    // Yg = as<vec>(model(theta,Verr));
    double beta2 = as<double>(LogTest(theta,Verr));
    double logR2 = beta2 - alpha2;
    if(log(as<double>(runif(1,0,1))) < logR2)
    {
      PHI.row(i+1)=phi_star.t();
      THETA.row(i+1)=theta_star.t();
      alpha2 = beta2;
      AcceptationRatio += 1;
    }
    else
    {
      PHI.row(i+1)=PHI.row(i);
      THETA.row(i+1)=THETA.row(i);
    }
    // if (i%100==0)
    // {
    //   if (AcceptationRatio/i<0.2)
    //   {
    //     t = t*0.9;
    //   }else
    //   {
    //     if (AcceptationRatio/i>0.5)
    //     {
    //       t = t*1.1;
    //     }
    //   }
    // }
  }
  if (stream==1)
  {
    std::cout << std::endl;
    cout << "End of the Metropolis Hastings algorithm"<< endl;
  }
  return List::create(Named("PHIwg")=PHIwg,Named("THETAwg")=THETAwg,Named("PHI")=PHI,Named("THETA")=THETA,
                            Named("AcceptationRatio")=AcceptationRatio, Named("AcceptationRatioWg")=AcceptationRatioWg
                        , Named("S")=S);
  }
  else
  {
  return List::create(Named("PHIwg")=PHIwg,Named("THETAwg")=THETAwg,Named("AcceptationRatioWg")=AcceptationRatioWg
                        , Named("S")=S);
  }
  return 0;
}


// [[Rcpp::export]]
arma::mat invMat(arma::mat V)
{
  int n=V.n_rows;
  arma::mat Vinv=solve(V, eye(n,n), solve_opts::fast);
  return Vinv;
}

// [[Rcpp::export]]
void FlushCPP()
{
  cout.flush();
}

// [[Rcpp::export]]
List MetropolisHastingsCppD(int Ngibbs, int Nmh, arma::vec theta_init, arma::vec k, arma::mat SIGMA, arma::vec Yf,
                           arma::vec binf, arma::vec bsup, Function LogTest, int stream)
{
  double Dim = theta_init.size();
  int D;
  arma::mat PHIwg=randu<arma::mat>(Ngibbs,Dim), THETAwg=randu<arma::mat>(Ngibbs,Dim);
  if (Nmh!=0) {D=Nmh;} else {D=10;}
  arma::mat PHI= randu<arma::mat>(D,Dim), THETA=randu<arma::mat>(D,Dim);
  THETA.row(0)=theta_init.t();
  PHI.row(0)= log((THETA.row(0).t()-binf)/(bsup-binf)).t();
  double AcceptationRatio=0;
  arma::vec AcceptationRatioWg=zeros(Dim,1);
  Function unscale("unscale"), rnorm("rnorm"), mvrnorm("mvrnorm"), runif("runif"), DefPos("DefPos");
  THETAwg.row(0)=theta_init.t();
  PHIwg.row(0) = log((THETAwg.row(0).t()-binf)/(bsup-binf)).t();
  arma::vec theta=theta_init.rows(0,Dim-4);
  arma::vec thetaD=theta_init.rows(Dim-3,Dim-2);
  double Verr=THETAwg(0,(Dim-1));
  //Rcpp::List res = as<Rcpp::List>(model(theta,thetaD,Verr));
  //arma::vec Yg=res["y"];
  // arma::vec Yg=as<arma::vec>(model(theta,Verr));
  double alpha = as<double>(LogTest(theta,thetaD,Verr));
  double alpha2 = alpha;
  if (stream==1)
  {
    cout << "Begin of the Metropolis within Gibbs algorithm" << endl;
    cout << "Number of iterations "<< Ngibbs << endl;
  }
  int barWidth = 40;
  int q = 0;
  for (int i=0; i<(Ngibbs-1); i++)
  {
    if (stream==1)
    {
      // beggining of the bar progress
      cout.flush();
      if ((i+2)%(Ngibbs/barWidth)==0)
      {
        float progress = (float)(i+2)/(float)(Ngibbs);
        cout << "[";
        q = ((float)(i+2)*(float)barWidth)/(float)(Ngibbs);
        if (q<barWidth) cout << string(q, '=');
        else if (q==barWidth) cout << string(q, '>');
        else cout << " ";
        cout << "] " << int(progress * 100.0) << " %\r";
      }
      // end bar progress
    }
    vec phi_star = PHIwg.row(i).t();
    vec theta_star = THETAwg.row(i).t();
    for (int j=0; j<Dim; j++)
    {
      if (j>0){
        phi_star.rows(0,j) = PHIwg.row(i+1).cols(0,j).t();
      }
      phi_star(j) = as<double>(rnorm(1,PHIwg(i,j),k(j)*SIGMA(j,j)));
      theta_star(j) = as<double>(unscale(exp(phi_star(j)),binf(j),bsup(j)));
      Verr = theta_star((Dim-1));
      thetaD= theta_star.rows((Dim-3),(Dim-2));
      theta = theta_init.rows(0,Dim-4);
      //Rcpp::List res = as<Rcpp::List>(model(theta,thetaD,Verr));
      //arma::vec Yg=res["y"];
      // Yg = as<vec>(model(theta_star.rows(0,(Dim-2)).t(),theta_star(Dim-1)));
      double beta = as<double>(LogTest(theta,thetaD,Verr));
      double logR = beta-alpha;
      if (log(as<double>(runif(1,0,1))) < logR)
      {
        PHIwg(i+1,j)=phi_star(j);
        THETAwg(i+1,j)=theta_star(j);
        alpha = beta;
        AcceptationRatioWg(j) += 1;
      }
      else
      {
        PHIwg(i+1,j)=PHIwg(i,j);
        THETAwg(i+1,j)=THETAwg(i,j);
      }
    }
    /*if (i%100==0)
    {
    cout<< "Iteration number :" << i << "/" << Ngibbs <<endl;
    cout<< "Theta " << THETAwg.row(i) <<endl;
    }*/
  }
  mat Stemp = cov(PHIwg.rows(10/100*Ngibbs,(Ngibbs-1)));
  mat S = as<arma::mat>(DefPos(Stemp));
  mat NewPhi=mean(PHIwg.rows(10/100*Ngibbs,(Ngibbs-1)));
  if (stream==1)
  {
    cout << endl;
    cout << endl;
    cout << "Estimation of the covariance matrix...." <<endl;
    cout << "End of the within gibbs algorithm"<< endl;
    cout << endl;
    /*cout << "The acceptance rate is: " << AcceptationRatioWg/Ngibbs << endl;*/
    cout << "Begin of the metropolis hastings algorithm using the covariance computed" << endl;
    cout << "Number of iterations "<< Nmh <<endl;
  }
  q = 0;
  if (Nmh!=0)
  {
    for (int i=0; i<(Nmh-1); i++)
    {
      if (stream==1)
      {
        // beggining of the bar progress
        cout.flush();
        if ((i+2)%(Nmh/barWidth)==0)
        {
          float progress = (float)(i+2)/(float)(Nmh);
          cout << "[";
          q = ((float)(i+2)*(float)barWidth)/(float)(Nmh);
          if (q<barWidth) cout << string(q, '=');
          else if (q==barWidth) cout << string(q, '>');
          else cout << " ";
          cout << "] " << int(progress * 100.0) << " %\r";
        }
        // end bar progress
      }
      vec phi_star = as<vec>(mvrnorm(1,NewPhi.t(),S));
      vec theta_star = as<vec>(unscale(exp(phi_star.t()),binf,bsup));
      theta = theta_init.rows(0,Dim-4);
      thetaD = theta_init.rows((Dim-3),(Dim-2));
      Verr = theta_star((Dim-1));
      //Rcpp::List res = as<Rcpp::List>(model(theta,thetaD,Verr));
      //arma::vec Yg=res["y"];
      // Yg = as<vec>(model(theta,Verr));
      double beta2 = as<double>(LogTest(theta,thetaD,Verr));
      double logR2 = beta2 - alpha2;
      if(log(as<double>(runif(1,0,1))) < logR2)
      {
        PHI.row(i+1)=phi_star.t();
        THETA.row(i+1)=theta_star.t();
        alpha2 = beta2;
        AcceptationRatio += 1;
      }
      else
      {
        PHI.row(i+1)=PHI.row(i);
        THETA.row(i+1)=THETA.row(i);
      }
    }
    if (stream==1)
    {
      std::cout << std::endl;
      cout << "End of the Metropolis Hastings algorithm"<< endl;
    }
    return List::create(Named("PHIwg")=PHIwg,Named("THETAwg")=THETAwg,Named("PHI")=PHI,Named("THETA")=THETA,
                              Named("AcceptationRatio")=AcceptationRatio, Named("AcceptationRatioWg")=AcceptationRatioWg
                          , Named("S")=S);
  }
  else
  {
    return List::create(Named("PHIwg")=PHIwg,Named("THETAwg")=THETAwg,Named("AcceptationRatioWg")=AcceptationRatioWg
                          , Named("S")=S);
  }
  return 0;
  }
