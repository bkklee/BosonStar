#include <iostream>

extern "C"
{
	void findGravRho(double *rho, double *phi, int lengthStep, double dx){
		double temp[lengthStep];
		double phip[lengthStep];
		double phiNew[lengthStep];
		double error[lengthStep];
		const double PI = 3.1415926535897932384626433832795;
		int relaxMax = 1000000;
		double tolerance = 1e-10;

		temp[0] = (4.0/3.0)*PI*((0.0+0.5)*dx)*((0.0+0.5)*dx)*((0.0+0.5)*dx)*rho[0];
		temp[1] = 9.0*temp[0] + (4.0/3.0)*PI*dx*((1.0+0.5)*dx)*((1.0+0.5)*dx)*rho[1];
		
		for(int j=2;j<lengthStep;j++){
			temp[j] = temp[j-2] + 4.0 * PI * (dx/3.0) *
						 (rho[j-2] * (((j-2.0)+0.5)*dx) * (((j-2.0)+0.5)*dx) +
					  	 4.0*rho[j-1]*(((j-1)+0.5)*dx)*(((j-1)+0.5)*dx) +
					  	 rho[j] * ((j+0.5)*dx) * ((j+0.5)*dx));
		}
		
		for(int j=0;j<lengthStep;j++){
			phip[j] = temp[j] / (((j+0.5)*dx)*((j+0.5)*dx));
		}

		//boundary1D(phip, true);
		
		phi[lengthStep-1] = 0.0;
		
		for(int j=lengthStep-2;j>=0;j--){
			phi[j] = phi[j+1] - (phip[j] + phip[j+1]) * dx/2.0;
		}

		//boundary1D(phi, false);
		
		//potentialRelax
		for(int n=0;n<relaxMax;n++){

			phiNew[0] = 0.5 * (phi[1] + phi[0]) + 0.5 * dx/((0.0+0.5)*dx) * (phi[1] - phi[0]) - 
							   2.0*PI*dx*dx*rho[0];
			for(int j=1;j<lengthStep-2;j++){
				phiNew[j] = 0.5 * (phi[j+1] + phi[j-1]) + 0.5 * dx/((j+0.5)*dx) * (phi[j+1] - phi[j-1]) - 
							   2.0*PI*dx*dx*rho[j];
			}
			phiNew[lengthStep-1] = phi[lengthStep-1];
			
			//boundary1D(phiNew, false);
			
			for(int j=0;j<lengthStep-1;j++){
				error[j] = (phiNew[j] - phi[j])/phi[j];
			}
			error[lengthStep-1] = 0.0;
			//boundary1D(error, false);
			
			for(int j=0;j<lengthStep;j++){
				phi[j] = phiNew[j];
			}
			
			phip[0] = (-phi[2] + 8.0*phi[1] - 8.0*phi[0] + phi[1])/(1.2e1*dx);
			phip[1] = (-phi[3] + 8.0*phi[2] - 8.0*phi[0] + phi[0])/(1.2e1*dx);
			phip[lengthStep-2] = (-phi[lengthStep-1] + 8.0*phi[lengthStep-1] - 8.0*phi[lengthStep-3] + phi[lengthStep-4])/(1.2e1*dx);
			phip[lengthStep-1] = (-phi[lengthStep-1] + 8.0*phi[lengthStep-1] - 8.0*phi[lengthStep-2] + phi[lengthStep-3])/(1.2e1*dx);

			for(int j=0;j<lengthStep;j++){
				phip[j] = (-phi[j+2] + 8.0*phi[j+1] - 8.0*phi[j-1] + phi[j-2])/(1.2e1*dx);
			}
			//boundary1D(phip, true);
			
			bool needRelax = false;
			
			for(int j=0;j<lengthStep;j++){
				if(abs(error[j]) > tolerance){
					needRelax = true;
				}
			}
			
			if(needRelax){
				continue;
			}else{
				break;
			}
			
		}
		
		return;

	}

}