#ifndef __MY_COUNTER__
#define __MY_COUNTER__

#include <fstream>
#include <armadillo>

class counter;

class counter
{
public:
	int nbin, sz_f;
	double dt;
	double *ek, *et, **pp, **ppd, *hop; // pp is population on adiabats, ppd is on diabats
	int *count;
	double *ek_all, *et_all, **pp_all, **ppd_all, *hop_all;
	int *count_all;
	//
	void init(int Sz, double Dt, double Tmax);
	void add(int niter, double eki, double eti, int istate, int hopi, arma::cx_mat U, arma::cx_mat rho);
	void average();
	void print_ene(std::ofstream& output);
	void print_pop(std::ofstream& output);
	void print_popd(std::ofstream& output);
};

#endif
