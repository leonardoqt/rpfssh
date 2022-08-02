#ifndef __POTENTIAL__
#define __POTENTIAL__

#include <armadillo>

class potential; // time independent part

class potential
{
public:
	const int sz_s = 2;
	const int dim = 2;
	const int sz_f = 1<<sz_s; // this is 2^sz_s
	const double dx = 1e-5, dp = 1e-5;
	int if_test=0;
	//
	double kx,chi,ky,delta,v0,w,mass;
	//
	double omega, gap, b, de;
	//
	double mul, mur;
	arma::vec vsbl, vsbr;
	void init_H(double Kx, double Chi, double Ky, double Delta, double V0, double W, double Mass, double Mul, double Mur, arma::vec Gammal, arma::vec Gammar);
	void init_H(double Omega, double Gap, double De, double W, double Mass, double Mul, double Mur, arma::vec Gammal, arma::vec Gammar);
	arma::cx_mat Hs(arma::vec x);
	void gen_Hs_pd(arma::vec x, arma::vec p, arma::mat& Hs_pd, arma::cx_mat& Us_pd);
	void gen_Hf_pd(arma::vec x, arma::vec p, arma::mat& Hf_pd, arma::cx_mat& Uf_pd);
	void ionic(arma::vec x, double& Eion, arma::vec& dEdx);
	void E_Hf_pd(arma::vec x, arma::vec p, arma::vec& Ef, arma::mat& Uf, arma::vec& Gammal, arma::vec& Gammar);
	// E_expect_Hf_pd: energy expectation of only the electronic part in superadiabats
	void E_expect_Hf_pd(arma::vec x, arma::vec p, arma::vec& Ef);
	void dyn_Hf_pd(arma::vec x, arma::vec p, arma::mat& dHdx, arma::mat& dHdp);
	arma::cx_mat ddt_f(arma::vec x1, arma::vec x2, arma::vec p1, arma::vec p2);
	//
};

#endif
