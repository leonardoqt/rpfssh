#include "counter.h"

using namespace std;

void counter::init(int Sz, double Dt, double Tmax)
{
	dt = Dt;
	nbin = ( (int)(Tmax/dt) ) + 10;
	sz_f = Sz;
	// TODO deallocate pointers before assign
	ek  = new double[nbin];
	et  = new double[nbin];
	pp  = new double*[sz_f];
	hop = new double[nbin];
	count  = new int[nbin];
	ek_all  = new double[nbin];
	et_all  = new double[nbin];
	pp_all  = new double*[sz_f];
	hop_all = new double[nbin];
	count_all = new int[nbin];
	for(int t1=0; t1<sz_f; t1++)
	{
		pp[t1] = new double[nbin];
		pp_all[t1] = new double[nbin];
	}
	for (int t1=0; t1<nbin; t1++)
	{
		ek[t1] = et[t1] = hop[t1] = 0;
		count[t1] = 0;
		for (int t2 = 0; t2<sz_f; t2++)
			pp[t2][t1] = 0;
	}
}

void counter::add(int niter, double eki, double eti, int istate, int hopi)
{
	ek[niter] += eki;
	et[niter] += eti;
	pp[istate][niter] += 1.0;
	hop[niter] += hopi;
	count[niter]++;
}

void counter::average()
{
	for(int t1=0; t1<nbin; t1++)
		if(count_all[t1]>0)
		{
			ek_all[t1] /= count_all[t1];
			et_all[t1] /= count_all[t1];
			hop_all[t1] /= count_all[t1];
			for(int t2=0; t2<sz_f; t2++)
				pp_all[t2][t1] /= count_all[t1];
		}
}

void counter::print_ene(ofstream& output)
{
	for(int t1=0; t1<nbin; t1++)
		if(count_all[t1]>0)
			output<<t1*dt<<'\t'<<ek_all[t1]<<'\t'<<et_all[t1]<<endl;
		else
			break;
}

void counter::print_pop(ofstream& output)
{
	for(int t1=0; t1<nbin; t1++)
		if(count_all[t1]>0)
		{
			output<<t1*dt<<'\t';
			for(int t2=0; t2<sz_f; t2++)
				output<<pp_all[t2][t1]<<'\t';
			output<<hop_all[t1]<<endl;
		}
		else
			break;
}
