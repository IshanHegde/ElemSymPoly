#ifndef JMLE_H

#define JMLE_H
#include <random>
#include <string>
#include <vector>
#include <future>



	struct LUdcmp
	{
		int n;
		std::vector<std::vector<double>> lu;
		std::vector<int> indx;
		double d;
		LUdcmp(std::vector<std::vector<double>> &a);
		void solve(std::vector<double> &b, std::vector<double> &x);
		void solve(std::vector<std::vector<double>> &b, std::vector<std::vector<double>> &x);
		void inverse(std::vector<std::vector<double>> &ainv);
		double det();
		void mprove(std::vector<double> &b, std::vector<double> &x);
		std::vector<std::vector<double>> &aref;
	};
	LUdcmp::LUdcmp(std::vector<std::vector<double> > &a) : n(a.size()), lu(a), aref(a){

		
		indx.reserve(n);
		const double TINY=1.0e-40;
		int i,imax,j,k;
		double big,temp;
		std::vector<double> vv(n);
		d=1.0;

		for (i=0;i<n;i++) {
			big=0.0;
			for (j=0;j<n;j++)
				if ((temp=abs(lu[i][j])) > big) big=temp;
			if (big == 0.0) throw("Singular matrix in LUdcmp");
			vv[i]=1.0/big;
		}
		for (k=0;k<n;k++) {
			big=0.0;
			imax=k;
			for (i=k;i<n;i++) {
				temp=vv[i]*abs(lu[i][k]);
				if (temp > big) {
					big=temp;
					imax=i;
				}
			}
			if (k != imax) {
				for (j=0;j<n;j++) {
					temp=lu[imax][j];
					lu[imax][j]=lu[k][j];
					lu[k][j]=temp;
				}
				d = -d;
				vv[imax]=vv[k];
			}
			indx[k]=imax;
			if (lu[k][k] == 0.0) lu[k][k]=TINY;
			for (i=k+1;i<n;i++) {
				temp=lu[i][k] /= lu[k][k];
				for (j=k+1;j<n;j++)
					lu[i][j] -= temp*lu[k][j];
			}
		}
	}
	void LUdcmp::solve(std::vector<double> &b, std::vector<double> &x)
	{
		int i,ii=0,ip,j;
		double sum;
		
		if (b.size() != n || x.size() != n)
			throw("LUdcmp::solve bad sizes");
		for (i=0;i<n;i++) 
			{
				
				x[i] = b[i];
			}
		for (i=0;i<n;i++) 
		{
			
			ip=indx[i];
			sum=x[ip];
			x[ip]=x[i];
			if (ii != 0)
				for (j=ii-1;j<i;j++) sum -= lu[i][j]*x[j];
			else if (sum != 0.0)
				ii=i+1;
			x[i]=sum;
		}
		for (i=n-1;i>=0;i--) {
			sum=x[i];
			for (j=i+1;j<n;j++) sum -= lu[i][j]*x[j];
			x[i]=sum/lu[i][i];
		}
		
	}
	void LUdcmp::solve(std::vector<std::vector<double> > &b, std::vector<std::vector<double> > &x)
	{
		int i,j,m=b[0].size();
		if (b.size() != n || x.size() != n || b[0].size() != x[0].size())
			throw("LUdcmp::solve bad sizes");
		std::vector<double> xx(n);
		for (j=0;j<m;j++) {
			for (i=0;i<n;i++) xx[i] = b[i][j];
			solve(xx,xx);
			for (i=0;i<n;i++) x[i][j] = xx[i];
		}
	}
	void LUdcmp::inverse(std::vector<std::vector<double> > &ainv)
	{
		int i,j;
		ainv.resize(n,std::vector<double>(n));
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) ainv[i][j] = 0.;
			ainv[i][i] = 1.;
		}
		solve(ainv,ainv);
	}
	double LUdcmp::det()
	{
		double dd = d;
		for (int i=0;i<n;i++) dd *= lu[i][i];
		return dd;
	}
	void LUdcmp::mprove(std::vector<double> &b, std::vector<double> &x)
	{
		int i,j;
		std::vector<double> r(n);
		for (i=0;i<n;i++) {
			double sdp = -b[i];
			for (j=0;j<n;j++)
				sdp += (double)aref[i][j] * (double)x[j];
			r[i]=sdp;
		}
		solve(r,r);
		for (i=0;i<n;i++) x[i] -= r[i];
	}

	template <class T>
	struct Vfdjac {
		const double EPS;
		T &func;
		Vfdjac(T &funcc) : EPS(1.0e-8),func(funcc) {}
		std::vector<std::vector<double>> operator() (const std::vector<double> &x, const std::vector<double> &fvec) {
			int n=x.size();
			std::vector<std::vector<double>> df(n,std::vector<double>(n));
			std::vector<double> xh=x;
			for (int j=0;j<n;j++) {
				double temp=xh[j];
				double h=EPS*abs(temp);
				if (h == 0.0) h=EPS;
				xh[j]=temp+h;
				h=xh[j]-temp;
				std::vector<double> f=func(xh);
				xh[j]=temp;

				std::vector<double> inner_store(n);
				
				for (int i=0;i<n;i++)
					inner_store[i]=(f[i]-fvec[i])/h;
				
				df[j]=inner_store;
			}
			return df;
		}
	};
	template <class T>
	struct Vfmin {
		std::vector<double> fvec;
		T &func;
		int n;
		Vfmin(T &funcc) : func(funcc){}
		double operator() (const std::vector<double> &x) {
			n=x.size();
			double sum=0;
			fvec=func(x);
			for (int i=0;i<n;i++) sum += fvec[i]*fvec[i];
			return 0.5*sum;
		}
	};

	template <class T>
	void lnsrch(const std::vector<double> &xold, const double fold, const std::vector<double> &g, std::vector<double> &p,
	std::vector<double> &x, double &f, const double stpmax, bool &check, T &func) {
		const double ALF=1.0e-4, TOLX=std::numeric_limits<double>::epsilon();
		double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
		double rhs1,rhs2,slope=0.0,sum=0.0,temp,test,tmplam;
		int i,n=xold.size();
		check=false;
		for (i=0;i<n;i++) sum += p[i]*p[i];
		sum=sqrt(sum);
		if (sum > stpmax)
			for (i=0;i<n;i++)
				p[i] *= stpmax/sum;
		for (i=0;i<n;i++)
			slope += g[i]*p[i];
		if (slope >= 0.0) throw("Roundoff problem in lnsrch.");
		test=0.0;
		for (i=0;i<n;i++) {
			temp=abs(p[i])/std::max(abs(xold[i]),1.0);
			if (temp > test) test=temp;
		}
		alamin=TOLX/test;
		alam=1.0;
		for (;;) {
			for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
			f=func(x);
			if (alam < alamin) {
				for (i=0;i<n;i++) x[i]=xold[i];
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
	void newt(std::vector<double> &x, bool &check, T &vecfunc) {
		const int MAXITS=200;
		const double TOLF=1.0e-8,TOLMIN=1.0e-12,STPMX=100.0;
		const double TOLX=std::numeric_limits<double>::epsilon();
		int i,j,its,n=x.size();
		double den,f,fold,stpmax,sum,temp,test;
		std::vector<double> g(n),p(n),xold(n);
		std::vector<std::vector<double>> fjac(n,std::vector<double>(n));
		Vfmin<T> Vfmin(vecfunc);
		Vfdjac<T> Vfdjac(vecfunc);
		std::vector<double> &fvec=Vfmin.fvec;
		f=Vfmin(x);
		test=0.0;
		for (i=0;i<n;i++)
			if (abs(fvec[i]) > test) test=abs(fvec[i]);
		if (test < 0.01*TOLF) {
			check=false;
			return;
		}
		sum=0.0;
		for (i=0;i<n;i++) sum += x[i]*x[i];
		stpmax=STPMX*std::max(sqrt(sum),double(n));
		for (its=0;its<MAXITS;its++) {
			fjac=Vfdjac(x,fvec);
			for (i=0;i<n;i++) {
				sum=0.0;
				for (j=0;j<n;j++) sum += fjac[j][i]*fvec[j];
				g[i]=sum;
			}
			for (i=0;i<n;i++) xold[i]=x[i];
			fold=f;
			for (i=0;i<n;i++) p[i] = -fvec[i];
			LUdcmp alu(fjac);
			alu.solve(p,p);
			lnsrch(xold,fold,g,p,x,f,stpmax,check,Vfmin);
			test=0.0;
			for (i=0;i<n;i++)
				if (abs(fvec[i]) > test) test=abs(fvec[i]);
			if (test < TOLF) {
				check=false;
				return;
			}
			if (check) {
				test=0.0;
				den=std::max(f,0.5*n);
				for (i=0;i<n;i++) {
					temp=abs(g[i])*std::max(abs(x[i]),1.0)/den;
					if (temp > test) test=temp;
				}
				check=(test < TOLMIN);
				return;
			}
			test=0.0;
			for (i=0;i<n;i++) {
				temp=(abs(x[i]-xold[i]))/std::max(abs(x[i]),1.0);
				if (temp > test) test=temp;
			}
			if (test < TOLX)
				return;
		}
		throw("MAXITS exceeded in newt");
	}



#endif