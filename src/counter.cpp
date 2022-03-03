#include "counter.h"

using namespace std;

void counter::init(double Dt, double Tmax)
{
	dt = Dt;
	nbin = ( (int)(Tmax/dt) ) + 10;
	// TODO deallocate pointers before assign
	ek  = new double[nbin];
	et  = new double[nbin];
	p0  = new double[nbin];
	p1  = new double[nbin];
	hop = new double[nbin];
	count  = new int[nbin];
	ek_all  = new double[nbin];
	et_all  = new double[nbin];
	p0_all  = new double[nbin];
	p1_all  = new double[nbin];
	hop_all = new double[nbin];
	count_all = new int[nbin];
	for (int t1=0; t1<nbin; t1++)
	{
		ek[t1] = et[t1] = p0[t1] = p1[t1] = hop[t1] = 0;
		count[t1] = 0;
	}
}

void counter::add(int niter, double eki, double eti, int istate, int hopi)
{
	ek[niter] += eki;
	et[niter] += eti;
	if (istate ==0)
		p0[niter] += 1.0;
	else
		p1[niter] += 1.0;
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
			p0_all[t1] /= count_all[t1];
			p1_all[t1] /= count_all[t1];
			hop_all[t1] /= count_all[t1];
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
			output<<t1*dt<<'\t'<<p0_all[t1]<<'\t'<<p1_all[t1]<<'\t'<<hop_all[t1]<<endl;
		else
			break;
}
