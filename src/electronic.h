#ifndef __ELECTRONIC_
#define __ELECTRONIC_

#include <armadillo>
#include "potential.h"

class electronic;

class electronic
{
public:
	int sz_s, sz_f, n_dt;
	double beta, dt;
	arma::cx_mat rho;
	arma::mat U;
	arma::mat hop_bath;
	//
	void init(arma::mat Rho0, potential& HH, double Beta, double Dt);
	void evolve(potential& HH, arma::vec x_t1, arma::vec x_t0, arma::vec p_t1, arma::vec p_t0);
	//void try_decoherence(ionic& AA);
private:
	void gen_n_dt(arma::mat T, arma::vec E_adiab);
	arma::cx_mat rho_dot(potential& HH, arma::vec E0, arma::mat Dvt, arma::vec Gl, arma::vec Gr, arma::cx_mat rho0);
	arma::cx_mat Lrho(potential& HH, arma::vec E0, arma::vec Gl, arma::vec Gr, arma::cx_mat rho0);
	void Lrho_ij(arma::cx_mat& res, arma::cx_mat& rho0, int i, int j, double coef);
	//void fit_drho_v1(potential& HH); // impose detailed balance, fit only diagonal
	//void fit_drho_v2(potential& HH); // impose detailed balance, fit full matrix
	void fit_drho_v3(); // impose single orbital the same, fit only diagonal
	void fit_drho_v3_1imp(); // impose single orbital the same, fit only diagonal
};
#endif
