#ifndef _GPU_H_
#define _GPU_H_

namespace GPU_interface_routines
{	
	// void setupTDsolve(int N, double* &d_ld, double* &d_d, double* &d_ud, double* &d_x, int device);
	
	// void destroyTDsolve(double* &d_ld, double* &d_d, double* &d_ud, double* &d_x);
	
	void TDsolve( int calculations_per_loop, int n_systems,
                            double *ld, double *dd, double *ud,      
                                  double *fin, int device);
							// double* &d_ld, double* &d_d, double* &d_ud,	double *&d_x);
                                  
};

class FokkerPlanckOnGPU 
{
    public:
        FokkerPlanckOnGPU();
        FokkerPlanckOnGPU(int calculations_per_loop, int n_systems, int device);
        ~FokkerPlanckOnGPU();
        void SolveTridiagonal(double *ld, double *dd, double *ud, double *fin);

    private:
    	int calc_per_loop, n_sys;
    	double *d_ld;
		double *d_d;
		double *d_ud;
		double *d_x;
};

#endif // #ifndef _GPU_H_