#ifndef _GPU_H_
#define _GPU_H_

namespace GPU_interface_routines
{	
    void FreeMatrixSystemOnHost(double *ld, double *dd, double *ud, double *fin);

	void AllocateMatrixSystemOnHost(int totalsize, 
                double *ld, double *dd, double *ud, double *fin);

	void TDsolve( int calculations_per_loop, int n_systems,
                            double *ld, double *dd, double *ud,      
                                  double *fin, int device);
							// double* &d_ld, double* &d_d, double* &d_ud,	double *&d_x);
                                  
};

class FokkerPlanckOnGPU 
{
    public:
        FokkerPlanckOnGPU();
        void initialize(int calculations_per_loop, int n_systems, int device);
        void destroy(int device);
        void SolveTridiagonal(double *ld, double *dd, double *ud, double *fin, int device);

    private:
    	int calc_per_loop, n_sys;
    	double *d_ld;
		double *d_d;
		double *d_ud;
		double *d_x;
};

#endif // #ifndef _GPU_H_