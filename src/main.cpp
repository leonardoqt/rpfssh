#include <iostream>
#include <fstream>
#include "potential.h"
#include "ionic.h"
#include "electronic.h"
#include "counter.h"
#include <mpi.h>
#include <time.h>
#include <chrono>

using namespace std;
using namespace arma;

int main()
{
	typedef chrono::high_resolution_clock clock;
	typedef chrono::high_resolution_clock::time_point timepoint;
	timepoint now = clock::now();
	MPI_Init(NULL,NULL);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	arma_rng::set_seed(now.time_since_epoch().count()+rank*10);
	//
	double gamma = 0.003, Temp = 0.03, lambda = 0.1, omega = 0.003, gap = 0.03, mur = 0.03, de = 0.03, w = 0.1;
	//
	double ek0 = 1e-3, ek1 = 1e-1;
	int nek = 60, state = 1, sample = 10000;
	double dt = 1.0, Tmax = 100000;
	//
	if ( rank == 0 )
		cin>>omega>>gap>>de>>w>>gamma>>mur>>Temp>>lambda>>ek0>>ek1>>nek>>sample>>dt>>Tmax;
	MPI_Bcast(&omega  ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&gap    ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&de     ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&w      ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&gamma  ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&mur    ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Temp   ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&lambda ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&ek0    ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&ek1    ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&nek    ,1,MPI_INT   ,0,MPI_COMM_WORLD);
	MPI_Bcast(&sample ,1,MPI_INT   ,0,MPI_COMM_WORLD);
	MPI_Bcast(&dt     ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Tmax   ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//====================
	potential HH;
	ionic AA;
	electronic EE;
	counter time_evo;
	//
	double mass = 1000.0;
	double beta = 1/Temp;
	double xstart = 0.0;
	double xend = 20.0;
	//double sigma_x = 0.5;
	//
	int sample_myself;
	double tmp,tmp2;
	double *x_final, *p_final, **p_finaln;
	double *x_final_all, *p_final_all, **p_finaln_all;
	//
	vec vv = linspace(sqrt(2*ek0/mass),sqrt(2*ek1/mass),nek);
	vec counter_t(HH.sz_f,fill::zeros), counter_r(HH.sz_f,fill::zeros);
	mat rho0(HH.sz_f,HH.sz_f,fill::zeros);
	rho0(state,state) = 1;
	//
	sample_myself = sample / size;
	x_final = new double[sample_myself];
	p_final = new double[sample_myself];
	x_final_all = new double[sample_myself*size];
	p_final_all = new double[sample_myself*size];
	p_finaln = new double*[4];
	p_finaln_all = new double*[4];
	for(int t1=0; t1<4; t1++)
	{
		p_finaln[t1] = new double[sample_myself];
		p_finaln_all[t1] = new double[sample_myself*size];
		for(int t2=0; t2<sample_myself; t2++)
			p_finaln[t1][t2] = 1e100;
	}
	//
	vec x0(2,fill::zeros), p0(2,fill::zeros);
	vec gammal(2,fill::zeros);
	gammal(0) = gammal(1) = gamma/4;
	HH.init_H(omega,gap,de,w,mass,0,mur,gammal,gammal);
	//-----------
	//vec t_xx = linspace(-20,20,1000);
	//vec t_x(2,fill::zeros), t_p(2,fill::zeros);
	//double t_eion;
	//vec t_eele,tmpv1,tmpv2;
	//mat tmpm1;
	//cout.precision(10);
	//cout.setf(ios::fixed);
	//for(auto m1:t_xx)
	//{
	//	t_x(0) = m1;
	//	HH.ionic(t_x,t_eion,tmpv1);
	//	HH.E_Hf_pd(t_x,t_p,t_eele,tmpm1,tmpv1,tmpv2);
	//	cout<<m1<<' ';
	//	(t_eele+t_eion).t().raw_print();
	//}
	//-----------
	for (int iv = 0; iv<nek; iv++)
	{
		time_evo.init(HH.sz_f,dt,Tmax);
		counter_t = counter_t*0;
		counter_r = counter_r*0;
		for (int isample=0; isample<sample_myself;isample++)
		{
			//TODO: no rand for testing
			//AA.init(HH,mass,vv(iv)+randn()*(0.5/sigma_x)/mass,xstart+randn()*sigma_x,state,-xend,xend);
			x0(0) = xstart;
			p0(0) = vv(iv)*mass;
			AA.init(mass,Temp,lambda,x0,p0,state,dt,-xend,xend);
			EE.init(rho0,HH,beta,dt);
			// run fssh
			for(int iter = 0; iter*EE.dt < Tmax; iter++)
			{
				AA.move(HH);
				EE.evolve(HH,AA.x_t1,AA.x,AA.p_t1,AA.p);
				//EE.try_decoherence(AA);
				AA.try_hop(HH,EE.rho,EE.hop_bath);
				time_evo.add(iter,AA.ek,AA.etot,AA.istate,AA.nhops,EE.U,EE.rho);
				//cout<<iter*time_evo.dt<<'\t'<<abs(EE.rho_fock(0,0))<<'\t'<<abs(EE.rho_fock(1,1))<<endl;
			}
			x_final[isample] = AA.x(0);
			p_final[isample] = AA.p(0);
			p_finaln[AA.istate][isample] = AA.p(0);
			// count rate
			if (AA.x(0) < 0)
				counter_r(AA.istate) += 1.0;
			else
				counter_t(AA.istate) += 1.0;
		}
		// collect counting
		for (int t1=0; t1<HH.sz_f; t1++)
		{
			tmp = counter_t(t1);
			MPI_Allreduce(&tmp,&tmp2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
			counter_t(t1) = tmp2;
			tmp = counter_r(t1);
			MPI_Allreduce(&tmp,&tmp2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
			counter_r(t1) = tmp2;
		}
		MPI_Gather(x_final,sample_myself,MPI_DOUBLE,x_final_all,sample_myself,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(p_final,sample_myself,MPI_DOUBLE,p_final_all,sample_myself,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(p_finaln[0],sample_myself,MPI_DOUBLE,p_finaln_all[0],sample_myself,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(p_finaln[1],sample_myself,MPI_DOUBLE,p_finaln_all[1],sample_myself,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(p_finaln[2],sample_myself,MPI_DOUBLE,p_finaln_all[2],sample_myself,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(p_finaln[3],sample_myself,MPI_DOUBLE,p_finaln_all[3],sample_myself,MPI_DOUBLE,0,MPI_COMM_WORLD);
		//
		MPI_Allreduce(time_evo.ek,time_evo.ek_all,time_evo.nbin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(time_evo.et,time_evo.et_all,time_evo.nbin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		for(int t1=0; t1<time_evo.sz_f; t1++)
		{
			MPI_Allreduce(time_evo.pp[t1],time_evo.pp_all[t1],time_evo.nbin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
			MPI_Allreduce(time_evo.ppd[t1],time_evo.ppd_all[t1],time_evo.nbin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		}
		MPI_Allreduce(time_evo.hop,time_evo.hop_all,time_evo.nbin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(time_evo.count,time_evo.count_all,time_evo.nbin,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		// print
		if (rank == 0)
		{
			time_evo.average();
			ofstream ff;
			//
			ff.open("ene_v-"+to_string(iv)+".dat");
			time_evo.print_ene(ff);
			ff.close();
			//
			ff.open("traj_v-"+to_string(iv)+".dat");
			time_evo.print_pop(ff);
			ff.close();
			//
			ff.open("dtraj_v-"+to_string(iv)+".dat");
			time_evo.print_popd(ff);
			ff.close();
			//
			ff.open("xx_v-"+to_string(iv)+".dat");
			for(int t1=0; t1<sample_myself*size; t1++)
				ff<<x_final_all[t1]<<endl;
			ff.close();
			//
			ff.open("px_v-"+to_string(iv)+".dat");
			for(int t1=0; t1<sample_myself*size; t1++)
				ff<<p_final_all[t1]<<endl;
			ff.close();
			//
			ff.open("npx0_v-"+to_string(iv)+".dat");
			for(int t1=0; t1<sample_myself*size; t1++)
				if(p_finaln_all[0][t1] < 1e99)
					ff<<p_finaln_all[0][t1]<<endl;
			ff.close();
			//
			ff.open("npx1_v-"+to_string(iv)+".dat");
			for(int t1=0; t1<sample_myself*size; t1++)
				if(p_finaln_all[1][t1] < 1e99)
					ff<<p_finaln_all[1][t1]<<endl;
			ff.close();
			//
			ff.open("npx2_v-"+to_string(iv)+".dat");
			for(int t1=0; t1<sample_myself*size; t1++)
				if(p_finaln_all[2][t1] < 1e99)
					ff<<p_finaln_all[2][t1]<<endl;
			ff.close();
			//
			ff.open("npx3_v-"+to_string(iv)+".dat");
			for(int t1=0; t1<sample_myself*size; t1++)
				if(p_finaln_all[3][t1] < 1e99)
					ff<<p_finaln_all[3][t1]<<endl;
			ff.close();
			//
			cout<<vv(iv)*vv(iv)*mass/2;
			for (int t1=0; t1<HH.sz_f; t1++)
				cout<<'\t'<<counter_t(t1)/sample_myself/size;
			for (int t1=0; t1<HH.sz_f; t1++)
				cout<<'\t'<<counter_r(t1)/sample_myself/size;
			cout<<endl;
		}
	}
	MPI_Finalize();
	return 0;
}
