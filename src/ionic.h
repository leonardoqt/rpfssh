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
	int hop_period;
	double mass, ek, etot;
	arma::vec x, p, x_t1, x_t2, p_t1, p_t2;
	double xendl, xendr;
	//
	void init(double Mass, arma::vec X, arma::vec P, int state, double Dt, int Period, double Xendl, double Xendr);
	void move(potential& HH, double tem);
	void try_hop(potential& HH, arma::cx_mat rho, arma::mat hop_bath);
	//int check_stop();
};

#endif
