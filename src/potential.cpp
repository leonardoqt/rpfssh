#include "potential.h"

using namespace arma;

void potential::init_H(double Kx, double Chi, double Ky, double Delta, double V0, double W, double Mass, double Mul, double Mur,vec Gammal, vec Gammar)
{
	if_test = 0;
	kx = Kx;
	chi = Chi;
	ky = Ky;
	delta = Delta;
	v0 = V0;
	w = W;
	mass = Mass;
	mul = Mul;
	mur = Mur;
	vsbl = sqrt(Gammal);
	vsbr = sqrt(Gammar);
}

void potential::init_H(double Omega, double Gap, double De, double W, double Mass, double Mul, double Mur, arma::vec Gammal, arma::vec Gammar)
{
	if_test = 1;
	omega  = Omega;
	gap = Gap;
	b = sqrt(4*omega*gap);
	de = De;
	w = W;
	mass = Mass;
	mul = Mul;
	mur = Mur;
	vsbl = sqrt(Gammal);
	vsbr = sqrt(Gammar);
}

cx_mat potential::Hs(vec x)
{
	cx_mat HH(sz_s,sz_s,fill::zeros);
	if(if_test)
	{
		HH(0,0) = cx_double(b*x(0)+de, 0);
		HH(1,1) = cx_double(-b*x(0)+de, 0);
		HH(0,1) = cx_double(cos(w*x(1)), sin(w*x(1)) ) * (1.5*gap*exp(-0.1*x(0)*x(0)));
		HH(1,0) = cx_double(cos(w*x(1)),-sin(w*x(1)) ) * (1.5*gap*exp(-0.1*x(0)*x(0)));
		return HH;
	}
	else
	{
		HH(0,0) = cx_double(x(0)+delta,0);
		HH(1,1) = cx_double(-x(0)-delta,0);
		HH(0,1) = cx_double( v0*cos(w*x(1)), v0*sin(w*x(1)) );
		HH(1,0) = cx_double( v0*cos(w*x(1)),-v0*sin(w*x(1)) );
		return HH;
	}
}

void potential::gen_Hs_pd(vec x, vec p, mat &Hs_pd, cx_mat &Us_pd)
{
	//if (if_test)
	//{
	//	Hs_pd = real(Hs(x)) + dot(p,p)/2/mass*eye(sz_s,sz_s);
	//	Us_pd = cx_mat(eye(sz_s,sz_s),zeros<mat>(sz_s,sz_s));
	//}
	//else
	//{
		Us_pd = zeros<cx_mat>(sz_s,sz_s);
		Us_pd(0,0) = cx_double(1,0);
		Us_pd(1,1) = cx_double(cos(w*x(1)),-sin(w*x(1)));
		Hs_pd = real(Us_pd.t()*Hs(x)*Us_pd);
		Hs_pd += p(0)*p(0)/2/mass*eye(sz_s,sz_s);
		mat p_tmp(sz_s,sz_s,fill::zeros);
		p_tmp(0,0) = p(1);
		p_tmp(1,1) = p(1)-w;
		Hs_pd += p_tmp*p_tmp/2/mass;
	//}
}

void potential::gen_Hf_pd(vec x, vec p, mat &Hf_pd, cx_mat &Uf_pd)
{
	mat Hs_pd;
	cx_mat Us_pd;
	gen_Hs_pd(x,p,Hs_pd,Us_pd);
	Hf_pd = zeros<mat>(sz_f,sz_f);
	Hf_pd(0,0) = dot(p,p)/2/mass;
	Hf_pd(3,3) = Hs_pd(0,0) + Hs_pd(1,1) - Hf_pd(0,0);
	Hf_pd(span(1,sz_s),span(1,sz_s)) = Hs_pd;
	//
	Uf_pd = cx_mat(eye(sz_f,sz_f),zeros<mat>(sz_f,sz_f));
	Uf_pd(span(1,sz_s),span(1,sz_s)) = Us_pd;
	Uf_pd(3,3) = Us_pd(1,1);
}

void potential::ionic(vec x, double& Eion, vec& dEdx)
{
	if (if_test)
	{
		Eion = omega*dot(x,x);
		dEdx = 2*omega*x;
	}
	else
	{
		Eion = x(0)*x(0)/2 + kx*x(0) + x(1)*x(1)*chi/2 + ky*x(1);
		dEdx = zeros<vec>(dim);
		dEdx(0) = x(0) + kx;
		dEdx(1) = chi*x(1) + ky;
	}
}

void potential::E_Hf_pd(vec x, vec p, vec& Ef, mat& Uf, vec& Gammal, vec& Gammar)
{
	vec Es;
	mat Hs_pd, Us;
	cx_mat Us_pd;
	gen_Hs_pd(x,p,Hs_pd,Us_pd);
	eig_sym(Es,Us,Hs_pd);
	Ef = zeros<vec>(sz_f);
	Ef(0) = dot(p,p)/2/mass;
	Ef.rows(1,sz_s) = Es;
	Ef(3) = Ef(1)+Ef(2)-Ef(0);
	Uf = eye(sz_f,sz_f);
	Uf(span(1,sz_s),span(1,sz_s)) = Us;
	//
	cx_vec vsbl_a = Us.t()*Us_pd.t()*vsbl;
	cx_vec vsbr_a = Us.t()*Us_pd.t()*vsbr;
	Gammal = square(abs(vsbl_a));
	Gammar = square(abs(vsbr_a));
}

void potential::dyn_Hf_pd(vec x, vec p, mat &dHdx, mat &dHdp)
{
	vec x_tmp, p_tmp;
	vec Ef1, Ef2, vt1, vt2;
	mat mt1;
	dHdx = zeros<mat>(dim,sz_f);
	dHdp = zeros<mat>(dim,sz_f);
	//
	for(int t1=0; t1<dim; t1++)
	{
		// d/dx
		x_tmp = x; x_tmp(t1) -= dx;
		E_Hf_pd(x_tmp,p,Ef1,mt1,vt1,vt2);
		x_tmp = x; x_tmp(t1) += dx;
		E_Hf_pd(x_tmp,p,Ef2,mt1,vt1,vt2);
		dHdx.row(t1) = (Ef2-Ef1).t()/2/dx;
		// d/dp
		p_tmp = p; p_tmp(t1) -= dp;
		E_Hf_pd(x,p_tmp,Ef1,mt1,vt1,vt2);
		p_tmp = p; p_tmp(t1) += dp;
		E_Hf_pd(x,p_tmp,Ef2,mt1,vt1,vt2);
		dHdp.row(t1) = (Ef2-Ef1).t()/2/dp;
	}
}

cx_mat potential::ddt_f(arma::vec x1, arma::vec x2, arma::vec p1, arma::vec p2)
{
	mat U1, U2;
	vec E,vt1,vt2;
	E_Hf_pd(x1,p1,E,U1,vt1,vt2);
	E_Hf_pd(x2,p2,E,U2,vt1,vt2);
	for(int t1=0; t1<sz_f; t1++)
		if (dot(U1.col(t1),U2.col(t1)) < 0)
			U2.col(t1) *= -1;
	return logmat(U1.t()*U2);
}
