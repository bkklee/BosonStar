#include <iostream>
#include <cmath>
using namespace std;

void iniDerBoson(double *der, double x, double *y, double mu){
	der[0] = x*x*y[1];
	der[1] = y[2];
	der[2] = y[3];
	der[3] = 2.0*mu*y[1]*y[2] + 2.0*y[0]*y[1]/x/x + 2.0/y[1]*y[2]*y[3] - 1.0/y[1]/y[1]*y[2]*y[2]*y[2] - 2.0/x*y[3] + 2.0/x/y[1]*y[2]*y[2] + 2.0/x/x*y[2];
	
	return;
}

void getRhoBoson1F(double scattering, double mBoson, double mu, double totalLength){

	const double PI = 3.1415926535897932384626433832795;
	const double hBar = 1.1977151493389159e-76;

	cout << "Solving for a pure boson star" << endl;
	cout << "Big Lambda: " << 2*scattering/mBoson << endl;
	
	const int noOfEqIni = 4;

	double psi = 4.0*scattering/mBoson;
	
	//Step Number
	int stepNum = 100000;
	double dxrk = 30.0/(stepNum-1.0);
	
	//Init
	double *yZero = new double[noOfEqIni];
	double *yOne = new double[noOfEqIni];
	double *yTwo = new double[noOfEqIni];
	double *yThree = new double[noOfEqIni];
	double *der = new double[noOfEqIni];

	double **y = new double *[noOfEqIni];
	for(int i=0;i<noOfEqIni;i++){
		y[i] = new double[stepNum];
	}
	
	double *shell = new double[stepNum];
	double *tempRho = new double [stepNum];
	double *r = new double [stepNum];
	
	int farPoint = -1;
	double minValue = 1.01;
	
	/////////////////////
	
	double start = -15.0;
	double end = -0.0001;
	
	while(abs(start-end) > 1e-15){
		
		double A[7];
		double dA = (end-start)/6.0;
		for(int i=0;i<7;i++){
			A[i] = dA*i + start;
		}
		
		double maxFarPoint = -1.0;
		int maxFarPointAA = -1;
		
		for(int aa=0;aa<7;aa++){
			
			y[0][0] = 0.0;
			y[1][0] = 1.0;
			y[2][0] = 0.0;
			y[3][0] = A[aa];
			
			for(int j=0;j<stepNum;j++){
				
				//k1
				for(int i=0;i<noOfEqIni;i++){
					yZero[i] = y[i][j];
				}
				
				double x = j*dxrk;
				
				iniDerBoson(der, x, yZero, mu);
				
				if(j==0){
					der[3] = 0.0;
				}
				
				double k1[4];
				
				for(int i=0;i<noOfEqIni;i++){
					k1[i] = dxrk*der[i];
					yOne[i] = y[i][j]+k1[i]/2.0;
				}
				
				//k2
				
				x = j*dxrk+(1.0/2.0)*dxrk;
				
				iniDerBoson(der, x, yOne, mu);
				
				double k2[4];
				
				for(int i=0;i<noOfEqIni;i++){
					k2[i] = dxrk*der[i];
					yTwo[i] = y[i][j]+k2[i]/2.0;
				}
				
				//k3
				
				x = j*dxrk+(1.0/2.0)*dxrk;
				
				iniDerBoson(der, x, yTwo, mu);
				
				double k3[4];
				for(int i=0;i<noOfEqIni;i++){
					k3[i] = dxrk*der[i];
					yThree[i] = y[i][j]+k3[i];
				}
				
				//k4
				
				x = j*dxrk+dxrk;
				
				iniDerBoson(der, x, yThree, mu);
				
				double k4[4];
				for(int i=0;i<noOfEqIni;i++){
					k4[i] = dxrk*der[i];
					y[i][j+1] = y[i][j] + k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0;
				}
				
			}
			
			farPoint = -1;
			minValue = 1.01;
			for(int qq=0;qq<stepNum;qq++){
				if(y[1][qq] < minValue){
					farPoint = qq;
					minValue = y[1][qq];
				}else{
					break;
				}
			}
			
			if(farPoint > maxFarPoint){
				maxFarPoint = farPoint;
				maxFarPointAA = aa;
			}
			
		}
		
		start = A[max(maxFarPointAA-1, 0)];
		end = A[min(maxFarPointAA+1, 6)];
	}
	
	//[Post processing]
	
	for(int qq=0;qq<stepNum;qq++){
		shell[qq] = y[1][qq]*qq*dxrk*qq*dxrk;
	}
	
	double volume = shell[0] + shell[farPoint-100];
	for(int i=1;i<farPoint-100;i++){
		if(i%3 == 0){
			volume += 2.0*shell[i];
		}else{
			volume += 3.0*shell[i];
		}
	}
	volume = 3.0*dxrk*volume/8.0;
	
	double n0 = pow(1.0/volume, 4.0);
	double chi = mu/sqrt(n0);
	double mStar = sqrt(chi*hBar*hBar/4.0/scattering/mBoson);
	double b = hBar*hBar/2.0/mStar/mBoson/mBoson;
	double rho0 = mStar*n0/4.0/PI/b/b/b;
	
	double scale = b/pow(n0, 1.0/4.0);
	
	cout << endl;
	cout << "Results of the fourth order ODE solver:" << endl;
	cout << "farPoint (0-1, check if too close to boundary values): " << (double)farPoint/(double)stepNum << endl;
	cout << "A2 (0.001-15.0, check if too close to boundary values): " << -start << endl;
	cout << "Furthest solved star: " << scale*30*(double)farPoint/(double)stepNum << endl;
	cout << "Simulation box size: " << totalLength << endl;
	
	return;
}

/////////////////////////////////////

int main(){

	getRhoBoson1F(1e-73, 2e-77, 0.318, 1000.0);

	return 0;
}