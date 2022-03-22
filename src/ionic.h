#ifndef __IONIC__
#define __IONIC__

#include <armadillo>
#include "potential.h"

class ionic;

class ionic
{
public:
	int istate;
	int nhops;
	double dt;
	double mass, ek, etot;
	double tem, lambda;
	arma::vec x, p, x_t1, p_t1;
	double xendl, xendr;
	//
	void init(double Mass, double Tem, double Lambda, arma::vec X, arma::vec P, int state, double Dt, double Xendl, double Xendr);
	void move(potential& HH);
	void try_hop(potential& HH, arma::cx_mat rho, arma::mat hop_bath);
	//int check_stop();
};

#endif
