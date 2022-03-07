#ifndef __MY_COUNTER__
#define __MY_COUNTER__

#include <fstream>

class counter;

class counter
{
public:
	int nbin, sz_f;
	double dt;
	double *ek, *et, **pp, *hop;
	int *count;
	double *ek_all, *et_all, **pp_all, *hop_all;
	int *count_all;
	//
	void init(int Sz, double Dt, double Tmax);
	void add(int niter, double eki, double eti, int istate, int hopi);
	void average();
	void print_ene(std::ofstream& output);
	void print_pop(std::ofstream& output);
};

#endif
