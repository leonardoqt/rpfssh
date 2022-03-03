#include "ionic.h"

using namespace arma;

void ionic::init(double Mass, vec X, vec P, int state, double Dt, int Period, double Xendl, double Xendr)
{
	dt = Dt;
	hop_period = Period;
	mass = Mass;
	x_t2 = x_t1 = x = X;
	p_t2 = p_t1 = p = P;
	ek = dot(p,p)/2/mass;
	istate = state;
	xendl = Xendl;
	xendr = Xendr;
	nhops = 0;
}

void ionic::move(potential& HH)
{
	// Sympletic integrator for separable Hamiltonian, i.e. grad U_pd is a constant
	// See https://en.wikipedia.org/wiki/Symplectic_integrator#:~:text=In%20concrete
	// TODO: change the algorithm for general Hamiltonian
	mat dHdx, dHdp;
	vec dEiondx;
	double Eion;
	//
	x_t2 = x_t1; x_t1 = x;
	p_t2 = p_t1; p_t1 = p;
	//
	HH.dyn_Hf_pd(x,p,dHdx,dHdp);
	HH.ionic(x,Eion,dEiondx);
	p = p - (dt/2)*(dHdx.col(istate)+dEiondx);
	//
	HH.dyn_Hf_pd(x,p,dHdx,dHdp);
	x = x + dHdp * dt;
	//
	HH.dyn_Hf_pd(x,p,dHdx,dHdp);
	HH.ionic(x,Eion,dEiondx);
	p = p - (dt/2)*(dHdx.col(istate)+dEiondx);
	//
	ek = dot(p,p)/2/mass;
}

void ionic::try_hop(potential &HH, arma::cx_mat rho, arma::mat hop_bath)
{
	cx_mat T = HH.ddt_f(x_t2,x,p_t2,p);
	vec rate_s(HH.sz_f), rate_b(HH.sz_f);
	//
	// TODO: rethink about this part
	//*************
	/*
	// rho_ii_dot can be splitted into two terms, the ration decides whether it undergoes the derivative coupling procedue or relaxation procedure
	cx_mat rho_dot1 = (T*rho - rho*T)/dt;
	double rho_ii_dot1 = real(rho_dot1(istate,istate));
	double rho_ii_dot2 = sum(hop_bath.row(istate))*real(rho(istate,istate));
	*/
	//*************
	//
	vec Ef, vt1, vt2, vt3;
	mat mt1;
	HH.E_Hf_pd(x,p,Ef,mt1,vt1,vt2);
	//
	for (int t1=0; t1<HH.sz_f; t1++)
	{
		if (t1==istate)
		{
			rate_s(t1) = 0;
			rate_b(t1) = 0;
		}
		else
		{
			rate_s(t1) = real( T(istate,t1)*rho(t1,istate) ) * 2 / real(rho(istate,istate));
			rate_b(t1) = ( hop_bath(t1,istate)*real(rho(istate,istate))-hop_bath(istate,t1)*real(rho(t1,t1)) )/real(rho(istate,istate)) * dt*hop_period;
		}
		if (rate_s(t1) < 0)
			rate_s(t1) = 0;
		if (rate_b(t1) < 0)
			rate_b(t1) = 0;
		if (ek + Ef(istate) < Ef(t1))
			rate_s(t1) = 0;
	}
	//rate_s.t().print();
	//rate_b.t().print();
	//cout<<endl;
	//
	vec rate = join_vert(rate_s,rate_b);
	for (int t1=1; t1<HH.sz_f*2; t1++)
		rate(t1) += rate(t1-1);
	//
	vec rnd(1,fill::randu);
	int new_state = istate;
	int from_bath = 1;
	for (int t1=0; t1<HH.sz_f*2; t1++)
		if( rnd(0) < rate(t1) )
		{
			new_state = t1 % HH.sz_f;
			from_bath = t1 / HH.sz_f;
			break;
		}
	//
	// adjust velocity
	if (from_bath == 0)
	{
		double dene = Ef(new_state) - Ef(istate);
		double dx = 1e-5;
		vec x_t1, x_t2;
		cx_mat DD;
		vec dp(HH.dim,fill::zeros);
		// for two impurities, rescale only occurs between 1,2
		for (int t1=0; t1<HH.dim; t1++)
		{
			x_t1 = x_t2 = x;
			x_t1(t1) -= dx;
			x_t2(t1) += dx;
			// TODO: check the meaning of DD and p dependence of DD
			DD = HH.ddt_f(x_t1,x_t2,p,p);
			dp(t1) = real(DD(1,2))/2/dx;
		}
		if (dot(p,dp)<0) dp*=-1;
		// p^2 = (p+a*dp)^2 + 2m*deme
		p += ( sqrt(4*dot(p,dp)*dot(p,dp)-8*mass*dene*dot(dp,dp)) - 2*dot(p,dp) ) / (2*dot(dp,dp)) * dp;
		ek = dot(p,p)/2/mass;
	}
	if (istate != new_state) nhops++;
	istate = new_state;
	double Eion;
	HH.ionic(x,Eion,vt1);
	etot = ek + Ef(istate) + Eion;
}

//int ionic::check_stop()
//{
//	if (x(0) < xendl)
//		return -1;
//	else if (x(0) > xendr)
//		return 1;
//	else
//		return 0;
//}
