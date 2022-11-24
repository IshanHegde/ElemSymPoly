#include <random>
#include <string>
#include "Eigen/Dense"
#include <future>
#include<iostream>

#ifndef NEWTON_H

#define NEWTON_H



struct LUdcmp
{
	int n;
	Eigen::MatrixXd lu;
	Eigen::VectorXd indx;
	double d;
	LUdcmp(Eigen::MatrixXd &a);
	void solve(Eigen::VectorXd &b, Eigen::VectorXd &x);
	void solve(Eigen::MatrixXd &b, Eigen::MatrixXd &x);
	void inverse(Eigen::MatrixXd &ainv);
	double det();
	void mprove(Eigen::VectorXd &b, Eigen::VectorXd &x);
	Eigen::MatrixXd &aref;
};
LUdcmp::LUdcmp(Eigen::MatrixXd &a):aref(a){

	n=a.rows();
	lu = a;

	indx = Eigen::VectorXd::Zero(n);
	const double TINY=1.0e-40;
	int i,imax,j,k;
	double big,temp;
	Eigen::VectorXd vv(n);
	d=1.0;

	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=abs(lu.coeff(i,j)) > big)) big=temp;
		if (big == 0.0) throw("Singular matrix in LUdcmp");
		vv(i)=1.0/big;
	}
	
	for (k=0;k<n;k++) {
		big=0.0;
		imax=k;
		for (i=k;i<n;i++) {
			temp=vv.coeff(i)*abs(lu.coeff(i,j));
			if (temp > big) {
				big=temp;
				imax=i;
			}
		}
		
		if (k != imax) {
			for (j=0;j<n;j++) {
				temp=lu(imax,j);
				lu(imax,j)=lu.coeff(k,j);
				lu(k,j)=temp;
			}
			d = -d;
			vv(imax)=vv.coeff(k);
		}
		
		indx(k)=imax;
		
		if (lu.coeff(k,k) == 0.0) lu(k,k)=TINY;
		
		for (i=k+1;i<n;i++) {
			
			temp=lu(i,k) /= lu(k,k);
			
			for (j=k+1;j<n;j++)
				{
					
				
					lu(i,j) -= temp*lu(k,j);
				}
				
		}
	}
}
void LUdcmp::solve(Eigen::VectorXd &b, Eigen::VectorXd &x)
{
	int i,ii=0,ip,j;
	double sum;
	
	if (b.size() != n || x.size() != n)
		throw("LUdcmp::solve bad sizes");
	for (i=0;i<n;i++) 
		{
			
			x(i) = b(i);
		}
	for (i=0;i<n;i++) 
	{
		
		ip=indx(i);
		sum=x(ip);
		x(ip)=x(i);
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= lu.coeff(i,j)*x(j);
		else if (sum != 0.0)
			ii=i+1;
		x(i)=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=x(i);
		for (j=i+1;j<n;j++) sum -= lu.coeff(i,j)*x(j);
		x(i)=sum/lu.coeff(i,i);
	}
	
}
void LUdcmp::solve(Eigen::MatrixXd &b, Eigen::MatrixXd &x)
{
	int i,j,m=b.cols();
	if (b.rows() != n || x.rows() != n || b.cols() != x.cols())
		throw("LUdcmp::solve bad sizes");
	Eigen::VectorXd xx(n);
	for (j=0;j<m;j++) {
		for (i=0;i<n;i++) xx(i) = b.coeff(i,j);
		solve(xx,xx);
		for (i=0;i<n;i++) x(i,j) = xx.coeff(i);
	}
}
void LUdcmp::inverse(Eigen::MatrixXd &ainv)
{
	int i,j;
	ainv.resize(n,n);
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) ainv(i,j) = 0.;
		ainv(i,i) = 1.;
	}
	solve(ainv,ainv);
}
double LUdcmp::det()
{
	double dd = d;
	for (int i=0;i<n;i++) dd *= lu.coeff(i,i);
	return dd;
}
void LUdcmp::mprove(Eigen::VectorXd &b, Eigen::VectorXd &x)
{
	int i,j;
	Eigen::VectorXd r(n);
	for (i=0;i<n;i++) {
		double sdp = -b.coeff(i);
		for (j=0;j<n;j++)
			sdp += aref.coeff(i,j) * x.coeff(j);
		r(i)=sdp;
	}
	solve(r,r);
	for (i=0;i<n;i++) x(i) -= r.coeff(i);
}

template <class T>
struct Vfdjac {
	const double EPS;
	T &func;
	Vfdjac(T &funcc) : EPS(1.0e-8),func(funcc) {}
	Eigen::MatrixXd  operator() (const Eigen::VectorXd &x, const Eigen::VectorXd &fvec) {
		int n=x.size();
		Eigen::MatrixXd  df(n,n);
		Eigen::VectorXd xh=x;
		for (int j=0;j<n;j++) {
			double temp=xh(j);
			double h=EPS*abs(temp);
			if (h == 0.0) h=EPS;
			xh(j)=temp+h;
			h=xh(j)-temp;
			Eigen::VectorXd f=func(xh);
			xh(j)=temp;
			
			for (int i=0;i<n;i++)
				{
					df(j,i)=(f[i]-fvec[i])/h;
				}
			
		}
		return df;
	}
};
template <class T>
struct Vfmin {
	Eigen::VectorXd fvec;
	T &func;
	int n;
	Vfmin(T &funcc) : func(funcc){}
	double operator() (Eigen::VectorXd &x) {
		n=x.size();
		double sum=0;
		fvec=func(x);
		for (int i=0;i<n;i++) sum += fvec(i)*fvec(i);
		return 0.5*sum;
	}
};

template <class T>
void lnsrch(const Eigen::VectorXd &xold, const double fold, const Eigen::VectorXd &g, Eigen::VectorXd &p,
Eigen::VectorXd &x, double &f, const double stpmax, bool &check, T &func) {
	const double ALF=1.0e-4, TOLX=std::numeric_limits<double>::epsilon();
	double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
	double rhs1,rhs2,slope=0.0,sum=0.0,temp,test,tmplam;
	int i,n=xold.size();
	check=false;
	for (i=0;i<n;i++) sum += p.coeff(i)*p.coeff(i);
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=0;i<n;i++)
			p(i) *= stpmax/sum;
	for (i=0;i<n;i++)
		slope += g.coeff(i)*p.coeff(i);
	if (slope >= 0.0) throw("Roundoff problem in lnsrch.");
	test=0.0;
	for (i=0;i<n;i++) {
		temp=abs(p.coeff(i))/std::max(abs(xold.coeff(i)),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (;;) {
		for (i=0;i<n;i++) x(i)=xold.coeff(i)+alam*p(i);
		f=func(x);
		if (alam < alamin) {
			for (i=0;i<n;i++) x(i)=xold.coeff(i);
			check=true;
			return;
		} else if (f <= fold+ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(f-fold-slope));
			else {
				rhs1=f-fold-alam*slope;
				rhs2=f2-fold-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc < 0.0) tmplam=0.5*alam;
					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
					else tmplam=-slope/(b+sqrt(disc));
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = f;
		alam=std::max(tmplam,0.1*alam);
	}
}

template <class T>
void newt(Eigen::VectorXd &x, bool &check, T &vecfunc) {
	const int MAXITS=200;
	const double TOLF=1.0e-8,TOLMIN=1.0e-12,STPMX=100.0;
	const double TOLX=std::numeric_limits<double>::epsilon();
	int i,j,its,n=x.rows();

	double den,f,fold,stpmax,sum,temp,test;
	Eigen::VectorXd g(n),p(n),xold(n);
	Eigen::MatrixXd fjac(n,n);

	Vfmin<T> Vfmin(vecfunc);
	Vfdjac<T> Vfdjac(vecfunc);

	Eigen::VectorXd &fvec=Vfmin.fvec;
	f=Vfmin(x);
	test=0.0;

	for (i=0;i<n;i++)
		if (abs(fvec.coeff(i)) > test) test=abs(fvec.coeff(i));
	if (test < 0.01*TOLF) {
		check=false;
		return;
	}
	sum=0.0;

	for (i=0;i<n;i++) sum += x.coeff(i)*x.coeff(i);
	stpmax=STPMX*std::max(sqrt(sum),double(n));
	for (its=0;its<MAXITS;its++) {

		fjac=Vfdjac(x,fvec);

		for (i=0;i<n;i++) {
			sum=0.0;
			for (j=0;j<n;j++)
			{
				sum += fjac(j,i)*fvec(j);
			}
			g(i)=sum;
		}

		for (i=0;i<n;i++) 
			{
				xold(i)=x(i);
			}
		fold=f;
		for (i=0;i<n;i++)
			{
				p(i) = -fvec.coeff(i);
			}

		LUdcmp alu(fjac);

		alu.solve(p,p);

		lnsrch(xold,fold,g,p,x,f,stpmax,check,Vfmin);

		test=0.0;
		for (i=0;i<n;i++)
			if (abs(fvec.coeff(i)) > test) test=abs(fvec.coeff(i));
		if (test < TOLF) {
			check=false;
			return;
		}
		if (check) {
			test=0.0;
			den=std::max(f,0.5*n);
			for (i=0;i<n;i++) {
				temp=abs(g.coeff(i))*std::max(abs(x.coeff(i)),1.0)/den;
				if (temp > test) test=temp;
			}
			check=(test < TOLMIN);
			return;
		}
		test=0.0;
		
		for (i=0;i<n;i++) {
			temp=(abs(x.coeff(i)-xold.coeff(i)))/std::max(abs(x.coeff(i)),1.0);
			
			if (temp > test) 
				{
					test=temp;
				}
		}
		if (test < TOLX)
			return;
	}
	throw("MAXITS exceeded in newt");
}



#endif