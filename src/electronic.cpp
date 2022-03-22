#include "electronic.h"
#include <iostream>

using namespace arma;

void electronic::init(mat Rho0, potential& HH, double Beta, double Dt)
{
	sz_s = HH.sz_s;
	sz_f = HH.sz_f;
	beta = Beta;
	dt = Dt;
	rho = cx_mat(Rho0,zeros<mat>(sz_f,sz_f));
	hop_bath = zeros<mat>(sz_f,sz_f);
}

void electronic::evolve(potential &HH, vec x_t1, vec x_t0, vec p_t1, vec p_t0)
{
	vec E_t1, E_t0;
	mat U_t1, U_t0;
	mat TT;
	vec Gl_t1, Gl_t0;
	vec Gr_t1, Gr_t0;
	//
	HH.E_Hf_pd(x_t1,p_t1,E_t1,U_t1,Gl_t1,Gr_t1);
	HH.E_Hf_pd(x_t0,p_t0,E_t0,U_t0,Gl_t0,Gr_t0);
	TT = real( HH.ddt_f(x_t1,x_t0,p_t1,p_t0) )/dt;
	TT = (TT - TT.t())/2;
	//p_t1.t().print();
	//p_t0.t().print();
	//TT.print();
	//
	gen_n_dt(TT,(E_t1+E_t0)/2);
	//
	cx_mat k1,k2,k3,k4;
	vec E2,E1,E0;
	vec Gl2,Gl1,Gl0;
	vec Gr2,Gr1,Gr0;
	for(int t1=0; t1<n_dt; t1++)
	{
		// use linear interpolation
		E2  = (E_t0 - E_t1)/n_dt * t1 + E_t1;
		E0  = E2 + (E_t0 - E_t1)/n_dt;
		E1  = (E2+E0)/2;
		Gl2 = (Gl_t0 - Gl_t1)/n_dt * t1 + Gl_t1;
		Gl0 = Gl2 + (Gl_t0 - Gl_t1)/n_dt;
		Gl1 = (Gl2+Gl0)/2;
		Gr2 = (Gr_t0 - Gr_t1)/n_dt * t1 + Gr_t1;
		Gr0 = Gr2 + (Gr_t0 - Gr_t1)/n_dt;
		Gr1 = (Gr2+Gr0)/2;
		//
		// Use RK4
		// rhodot = -i[H(x) rho] - [D(x) rho] + Lrho
		k1 = rho_dot(HH,E2,TT,Gl2,Gr2,rho);
		k2 = rho_dot(HH,E1,TT,Gl1,Gr1,rho + dt/n_dt/2*k1);
		k3 = rho_dot(HH,E1,TT,Gl1,Gr1,rho + dt/n_dt/2*k2);
		k4 = rho_dot(HH,E0,TT,Gl0,Gr0,rho + dt/n_dt  *k3);
		rho += (k1+2*k2+2*k3+k4)*(dt/n_dt/6);
		//
		// try hop
	}
	//
	U = U_t0;
}

void electronic::gen_n_dt(mat T, vec E_adiab)
{
	vec dt_e(3,fill::zeros);
	dt_e(0) = dt;
	dt_e(1) = 0.02 / ( abs(T).max() );
	dt_e(2) = 0.02 / ( abs(E_adiab - mean(E_adiab)).max() );
	n_dt = ceil( dt / dt_e.min() );
	//cout<<n_dt<<endl;
}

cx_mat electronic::rho_dot(potential& HH, vec E0, mat TT, vec Gl, vec Gr, cx_mat rho0)
{
	cx_mat res, iHD;
	iHD = cx_double(0,-1)*diagmat(E0) - TT;
	res = iHD*rho0 - rho0*iHD;
	if ( dot(Gl,Gl) + dot(Gr,Gr) > 1e-20 )
		res += Lrho(HH,E0,Gl,Gr,rho0);
	else
		hop_bath = zeros<mat>(sz_f,sz_f);
	return res;
}

cx_mat electronic::Lrho(potential& HH, vec E0, vec Gl, vec Gr, cx_mat rho0)
{
	vec fl(sz_s,fill::zeros), fr(sz_s,fill::zeros);
	fl = 1 / ( 1 + exp(beta*(E0.rows(1,sz_s)-E0(0)-HH.mul)) );
	fr = 1 / ( 1 + exp(beta*(E0.rows(1,sz_s)-E0(0)-HH.mur)) );
	//
	cube LL(sz_f,sz_f,8);
	LL(1,0,0) = LL(3,2,1) = LL(0,1,2) = LL(2,3,3) = 1;
	LL(2,0,4) = LL(3,1,5) = LL(0,2,6) = LL(1,3,7) = 1;
	//
	cx_mat res(sz_f,sz_f,fill::zeros);
	mat L;
	//
	L = LL.slice(0);
	res += (Gl(0)*fl(0)+Gr(0)*fr(0)) * ( L*rho0*L.t() - (L.t()*L*rho0+rho0*L.t()*L)/2 );
	L = LL.slice(1);
	res += (Gl(0)*fl(0)+Gr(0)*fr(0)) * ( L*rho0*L.t() - (L.t()*L*rho0+rho0*L.t()*L)/2 );
	L = LL.slice(2);
	res += (Gl(0)*(1-fl(0))+Gr(0)*(1-fr(0))) * ( L*rho0*L.t() - (L.t()*L*rho0+rho0*L.t()*L)/2 );
	L = LL.slice(3);
	res += (Gl(0)*(1-fl(0))+Gr(0)*(1-fr(0))) * ( L*rho0*L.t() - (L.t()*L*rho0+rho0*L.t()*L)/2 );
	//
	L = LL.slice(4);
	res += (Gl(1)*fl(1)+Gr(1)*fr(1)) * ( L*rho0*L.t() - (L.t()*L*rho0+rho0*L.t()*L)/2 );
	L = LL.slice(5);
	res += (Gl(1)*fl(1)+Gr(1)*fr(1)) * ( L*rho0*L.t() - (L.t()*L*rho0+rho0*L.t()*L)/2 );
	L = LL.slice(6);
	res += (Gl(1)*(1-fl(1))+Gr(1)*(1-fr(1))) * ( L*rho0*L.t() - (L.t()*L*rho0+rho0*L.t()*L)/2 );
	L = LL.slice(7);
	res += (Gl(1)*(1-fl(1))+Gr(1)*(1-fr(1))) * ( L*rho0*L.t() - (L.t()*L*rho0+rho0*L.t()*L)/2 );
	//
	hop_bath(1,0) = hop_bath(3,2) = Gl(0)*fl(0)+Gr(0)*fr(0);
	hop_bath(2,0) = hop_bath(3,1) = Gl(1)*fl(1)+Gr(1)*fr(1);
	hop_bath(0,1) = hop_bath(2,3) = Gl(0)*(1-fl(0))+Gr(0)*(1-fr(0));
	hop_bath(0,2) = hop_bath(1,3) = Gl(1)*(1-fl(1))+Gr(1)*(1-fr(1));
	return res;
}
//void electronic::try_decoherence(ionic& AA)
//{
//	if (AA.v_pre*AA.v_new <0 && AA.istate >0)
//	{
//		psi = psi * 0;
//		psi(AA.istate) = cx_double(1,0);
//		rho = psi * psi.t();
//	}
//}
