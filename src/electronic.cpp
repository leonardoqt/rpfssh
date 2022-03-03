#include "electronic.h"
#include <iostream>

using namespace arma;

void electronic::init(mat Rho0, potential& HH, double Beta, double Dt)
{
	first_evo = 1;
	sz_s = HH.sz_s;
	sz_f = HH.sz_f;
	beta = Beta;
	dt = Dt;
	rho = cx_mat(Rho0,zeros<mat>(sz_f,sz_f));
	hop_bath = zeros<mat>(sz_f,sz_f);
}

void electronic::evolve(potential &HH, vec x_t2, vec x_t1, vec x_t0, vec p_t2, vec p_t1, vec p_t0)
{
	// Use RK4
	// Ndot = -i[H rho] - damping*(rho_b - rho_eq)
	// rhodot = -i[H(x) rho] - [D(x) rho] + 
	if (first_evo)
	{
		first_evo = 0;
		x_last = 2*x_t2 - x_t1;
		p_last = 2*p_t2 - p_t1;
	}
	vec E_t2, E_t1, E_t0;
	mat U_t3, U_t2, U_t1, U_t0;
	cx_mat Dvt_t2, Dvt_t1, Dvt_t0;
	vec Gl_t2, Gl_t1, Gl_t0;
	vec Gr_t2, Gr_t1, Gr_t0;
	//
	HH.E_Hf_pd(x_last,p_last,E_t2,U_t3,Gl_t2,Gr_t2);
	HH.E_Hf_pd(x_t2,p_t2,E_t2,U_t2,Gl_t2,Gr_t2);
	HH.E_Hf_pd(x_t1,p_t1,E_t1,U_t1,Gl_t1,Gr_t1);
	HH.E_Hf_pd(x_t0,p_t0,E_t0,U_t0,Gl_t0,Gr_t0);
	Dvt_t2 = logmat(U_t3.t()*U_t2);
	Dvt_t1 = logmat(U_t2.t()*U_t1);
	Dvt_t0 = logmat(U_t1.t()*U_t0);
	//
	cx_mat k1,k2,k3,k4;
	k1 = rho_dot(HH,E_t2,Dvt_t2,Gl_t2,Gr_t2,rho);
	k2 = rho_dot(HH,E_t1,Dvt_t1,Gl_t1,Gr_t1,rho + dt/2*k1);
	k3 = rho_dot(HH,E_t1,Dvt_t1,Gl_t1,Gr_t1,rho + dt/2*k2);
	k4 = rho_dot(HH,E_t0,Dvt_t0,Gl_t0,Gr_t0,rho + dt  *k3);
	rho += (k1+2*k2+2*k3+k4)*(dt/6);
	//
	// save last xp for next RK4
	x_last = x_t1;
	p_last = p_t1;
}

cx_mat electronic::rho_dot(potential& HH, vec E0, cx_mat Dvt, vec Gl, vec Gr, cx_mat rho0)
{
	cx_mat res, iHD;
	iHD = cx_double(0,-1)*diagmat(E0) - Dvt*2/dt;
	res = iHD*rho0 - rho0*iHD;
	res += Lrho(HH,E0,Gl,Gr,rho0);
	return res;
}

cx_mat electronic::Lrho(potential& HH, vec E0, vec Gl, vec Gr, cx_mat rho0)
{
	vec fl(sz_s,fill::zeros), fr(sz_s,fill::zeros);
	fl = 1 / ( 1 + exp(beta*(E0.rows(1,sz_s)-HH.mul)) );
	fr = 1 / ( 1 + exp(beta*(E0.rows(1,sz_s)-HH.mur)) );
	//
	cube LL(sz_f,sz_f,8);
	LL(1,0,0) = LL(3,2,1) = LL(0,1,2) = LL(2,3,3) = 1;
	LL(2,0,4) = LL(3,1,5) = LL(0,2,6) = LL(1,3,7) = 1;
	//
	cx_mat res(sz_f,sz_f,fill::zeros);
	mat L;
	//
	L = LL.slice(0);
	res += (Gl(0)*fl(0)+Gr(0)*fr(0)) * ( L*rho0*L.t() - (L.t()*L*rho+rho*L.t()*L)/2 );
	L = LL.slice(1);
	res += (Gl(0)*fl(0)+Gr(0)*fr(0)) * ( L*rho0*L.t() - (L.t()*L*rho+rho*L.t()*L)/2 );
	L = LL.slice(2);
	res += (Gl(0)*(1-fl(0))+Gr(0)*(1-fr(0))) * ( L*rho0*L.t() - (L.t()*L*rho+rho*L.t()*L)/2 );
	L = LL.slice(3);
	res += (Gl(0)*(1-fl(0))+Gr(0)*(1-fr(0))) * ( L*rho0*L.t() - (L.t()*L*rho+rho*L.t()*L)/2 );
	//
	L = LL.slice(4);
	res += (Gl(1)*fl(1)+Gr(1)*fr(1)) * ( L*rho0*L.t() - (L.t()*L*rho+rho*L.t()*L)/2 );
	L = LL.slice(5);
	res += (Gl(1)*fl(1)+Gr(1)*fr(1)) * ( L*rho0*L.t() - (L.t()*L*rho+rho*L.t()*L)/2 );
	L = LL.slice(6);
	res += (Gl(1)*(1-fl(1))+Gr(1)*(1-fr(1))) * ( L*rho0*L.t() - (L.t()*L*rho+rho*L.t()*L)/2 );
	L = LL.slice(7);
	res += (Gl(1)*(1-fl(1))+Gr(1)*(1-fr(1))) * ( L*rho0*L.t() - (L.t()*L*rho+rho*L.t()*L)/2 );
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
