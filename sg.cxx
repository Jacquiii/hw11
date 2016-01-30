#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void step(cmplx* psi1,cmplx*psi0, const double alpha, const double lambda, const double omega, const double dx, const double dt,
          const int Nx, const double xmin);

void init( cmplx* psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);
//-----------------------------------
int main(){

	const int Nx = 300 ;
	const double xmin = -40;
        const double xmax = 40 ;
	const double dx =(xmax-xmin)/(Nx-1) ;

        const double Tend = 10*M_PI;
	const double dt =  0.1*dx ;
        double t = 0;
	
        const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
        const double omega = 0.2;
        const double alpha = sqrt(omega); 


  stringstream strm;

	cmplx* psi1 = new cmplx[Nx];
	cmplx* psi0 = new cmplx[Nx];
	cmplx* h = new cmplx[Nx];

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
                
		step(psi1,psi0,alpha,lambda,omega, dx,dt,Nx,xmin);
		
		h=psi0;
		psi0=psi1;
		psi1=h;
		
                t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;
   
   delete[] psi0;
   delete[] psi1;
	return 0;
}
//-----------------------------------
void step(cmplx* psi1, cmplx* psi0, const double alpha, const double lambda, const double omega, const double dx, const double dt,
          const int Nx, const double xmin)
{

  cmplx* d  = new cmplx[Nx];
  cmplx* au = new cmplx[Nx];
  cmplx* ad = new cmplx[Nx];
  
  cmplx* dconj  = new cmplx[Nx];
  cmplx* auconj = new cmplx[Nx];
  cmplx* adconj = new cmplx[Nx];
  
  
  for(int i=0;i<Nx;i++){
    double x = xmin + i*dx;
    d[i] = 1.0 + cmplx(0.0, dt/(2*dx*dx) + dt/4*pow(omega,2.)*pow(x,2.));
    dconj[i] = conj(d[i]);
  
   au[i] = cmplx(0.0, -dt/(4*dx*dx));
   auconj[i] = conj(au[i]);
   
   ad[i] = cmplx(0.0, -dt/(4*dx*dx));
   adconj[i] = conj(ad[i]);
  }
  
  for(int i=0;i<Nx-1;i++){
    cmplx c = ad[i+1]/d[i];
    
    d[i+1]=d[i+1]- c *au[i];
    ad[i+1]=0;
    
    adconj[i+1]=adconj[i+1]-c*dconj[i];
    dconj[i+1]=dconj[i+1]-c*auconj[i];
    

  }
  //Backward substitution

  psi1[Nx-1]=(adconj[Nx-1]*psi0[Nx-2] + dconj[Nx-1]*psi0[Nx-1])/d[Nx-1]; //letzte Zeile
  
  for(int i=1; i<Nx-1; i++){  
   psi1[Nx-1-i]=(psi0[Nx-2-i] * adconj[Nx-1-i]
                +psi0[Nx-1-i] * dconj [Nx-1-i]
                +psi0[Nx-i]   * auconj[Nx-1-i]
                -psi1[Nx-i]   * au[Nx-1-i])/(d[Nx-1-i]);
   }
 
  psi1[0] = (auconj[0]*psi0[1] + dconj[0]*psi0[0] - au[0]*psi1[1]) / d[0]; //erste Zeile
  
  delete[] d;
  delete[] au;
  delete[] ad;
  
  delete[] dconj;
  delete[] auconj;
  delete[] adconj;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
