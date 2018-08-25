/*! \brief Numerical Methods - Declarations
* \author  PICKSC
 * \date   2018
 * \file   nmethods.cpp
 * 
 * This cpp file contains the definitions for the functions
 * required for the numerical methods. 
 */  


//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <algorithm>
    #include <cstdlib>
    #include <float.h>

    #include <math.h>
    #include <map>

//  My libraries
    #include "lib-array.h"
    #include "input.h"


//-------------------------------------------------------------------
     bool Gauss_Seidel(Array2D<double>& A, 
                       valarray<complex<double> >& b,
                       valarray<complex<double> >& xk) {
//-------------------------------------------------------------------
        double     tol(1.0e-1);   //< Tolerance for absolute error
        size_t        MAXiter(5);    //< Maximum iteration allowed


//      The Matrices all have the right dimensions
//      -------------------------------------------------------------
        if ( ( A.dim1() != A.dim2()  ) || 
             ( A.dim1() != b.size()  ) ||
             ( A.dim1() != xk.size() )    )  {
            cout << "Error: The Matrices don't have the right dimensions!" << endl;
            exit(1);
        }
//      Check if the matrix A is diagonally dominant
//      -------------------------------------------------------------
        for (size_t i(0); i < A.dim1(); ++i){
            double rowi(0.0);
            for (size_t j(0); j < A.dim2(); ++j){
                rowi += A(i,j);
            }
            if (!(rowi < 2.0*A(i,i))) return false;
        }
//      -------------------------------------------------------------


//      Calculate and invert the diagonal elements once and for all
//      -------------------------------------------------------------      
        valarray<double> invDIAG(A.dim1());
        for (size_t i(0); i < invDIAG.size(); ++i) invDIAG[i] = 1.0 / A(i,i);

        valarray<complex<double> > xold(xk);
        size_t iteration(0);         // used to count iterations
        size_t conv(0);              // used to test convergence

//      Start the iteration loop
//      -------------------------------------------------------------      
        while ( (iteration++ < MAXiter) && (conv < b.size()) ) {

            xold = xk;
            for (size_t i(0); i < A.dim1(); ++i){
                complex<double> sigma(0.0);    // Temporary sum
                for (size_t j(0); j < i; ++j){
                    sigma += A(i,j)*xk[j];   
                }
                for (size_t j(i+1); j < A.dim2(); ++j){
                    sigma += A(i,j)*xk[j];   
                }
                xk[i] = invDIAG[i] * (b[i] - sigma);
            }

            // Calculate Dx = x_old - x_new
            xold -= xk;

//          If the relative error < prescribed tolerance everywhere the method has converged
//          |Dx(i)| < t*|x(i)| + eps 
            conv = 0;
            while ( ( conv < b.size() ) && 
                    ( abs(xold[conv]) < (tol*abs(xk[conv] + 20.0*DBL_MIN)) ) ){ 
                ++conv;
            } 

            //----> Output for testing
            //--------------------------------
            //cout << "iteration = " << iteration << "    ";
            //for (size_t i(0); i < b.size(); ++i){
            //    cout << xk[i] << "        ";
            //}
            //cout << "\n";
            //--------------------------------

        }

        //----> Output for testing
        //--------------------------------
        // cout << "Iterations = " << iteration-1  <<"\n";
        //for (size_t i(0); i < b.size(); ++i) {
        //    cout << "Error |Dx| = " << abs(xold[i]) 
        //         << ",    " 
        //         << "Tolerance * |x| = " << tol*abs(xk[i]) <<"\n";
        //}
        //--------------------------------

        return true;
    }
//-------------------------------------------------------------------
//-------------------------------------------------------------------
     bool Gauss_Seidel(Array2D<double>& A, 
                       valarray<double >& b,
                       valarray<double >& xk) {
//-------------------------------------------------------------------
        double     tol(1.0e-1);   //< Tolerance for absolute error
        size_t        MAXiter(5);    //< Maximum iteration allowed


//      The Matrices all have the right dimensions
//      -------------------------------------------------------------
        if ( ( A.dim1() != A.dim2()  ) || 
             ( A.dim1() != b.size()  ) ||
             ( A.dim1() != xk.size() )    )  {
            cout << "Error: The Matrices don't have the right dimensions!" << endl;
            exit(1);
        }
//      Check if the matrix A is diagonally dominant
//      -------------------------------------------------------------
        for (size_t i(0); i < A.dim1(); ++i){
            double rowi(0.0);
            for (size_t j(0); j < A.dim2(); ++j){
                rowi += A(i,j);
            }
            if (!(rowi < 2.0*A(i,i))) return false;
        }
//      -------------------------------------------------------------


//      Calculate and invert the diagonal elements once and for all
//      -------------------------------------------------------------      
        valarray<double> invDIAG(A.dim1());
        for (size_t i(0); i < invDIAG.size(); ++i) invDIAG[i] = 1.0 / A(i,i);

        valarray<double > xold(xk);
        size_t iteration(0);         // used to count iterations
        size_t conv(0);              // used to test convergence

//      Start the iteration loop
//      -------------------------------------------------------------      
        while ( (iteration++ < MAXiter) && (conv < b.size()) ) {

            xold = xk;
            for (size_t i(0); i < A.dim1(); ++i){
                double sigma(0.0);    // Temporary sum
                for (size_t j(0); j < i; ++j){
                    sigma += A(i,j)*xk[j];   
                }
                for (size_t j(i+1); j < A.dim2(); ++j){
                    sigma += A(i,j)*xk[j];   
                }
                xk[i] = invDIAG[i] * (b[i] - sigma);
            }

            // Calculate Dx = x_old - x_new
            xold -= xk;

//          If the relative error < prescribed tolerance everywhere the method has converged
//          |Dx(i)| < t*|x(i)| + eps 
            conv = 0;
            while ( ( conv < b.size() ) && 
                    ( abs(xold[conv]) < (tol*abs(xk[conv] + 20.0*DBL_MIN)) ) ){ 
                ++conv;
            } 

//             //----> Output for testing
//             //--------------------------------
//             //cout << "iteration = " << iteration << "    ";
//             //for (size_t i(0); i < b.size(); ++i){
//             //    cout << xk[i] << "        ";
//             //}
//             //cout << "\n";
//             //--------------------------------

        }

//         //----> Output for testing
//         //--------------------------------
//         // cout << "Iterations = " << iteration-1  <<"\n";
//         //for (size_t i(0); i < b.size(); ++i) {
//         //    cout << "Error |Dx| = " << abs(xold[i]) 
//         //         << ",    " 
//         //         << "Tolerance * |x| = " << tol*abs(xk[i]) <<"\n";
//         //}
//         //--------------------------------

        return true;
    }
//-------------------------------------------------------------------
//*******************************************************************
//-------------------------------------------------------------------
     void TridiagonalSolve (valarray<double>& a, 
                            valarray<double>& b, 
                                  valarray<double>& c,      
                                  valarray<complex<double> >&  d,
                                  valarray<complex<double> >& x) {
//-------------------------------------------------------------------
//   Fills solution into x. Warning: will modify c and d! 
//-------------------------------------------------------------------
        size_t n(x.size());
	   // Modify the coefficients. 
    	c[0] /= b[0];                            // Division by zero risk. 
    	d[0] /= b[0];                            // Division by zero would imply a singular matrix. 
        // 
    	for (size_t i(1); i < n; ++i){
            double id(1.0/(b[i]-c[i-1]*a[i]));   // Division by zero risk. 
    	    c[i] *= id;	                         // Last value calculated is redundant.
    	    d[i] -= d[i-1] * a[i];
    	    d[i] *= id;                          // d[i] = (d[i] - d[i-1] * a[i]) * id 
    	}
     
    	// Now back substitute. 
    	x[n-1] = d[n-1];
        for (int i(2); i < n+1; ++i)
        {
            x[n-i]  = d[n-i];
            x[n-i] -= c[n-i] * x[n-i+1];               // x[i] = d[i] - c[i] * x[i + 1];
        }
    }
//-------------------------------------------------------------------
//*******************************************************************
//-------------------------------------------------------------------
     void TridiagonalSolve ( size_t calculations_per_loop,
                            valarray<double>& a, 
                            valarray<double>& b, 
                                  valarray<double>& c,      
                                  valarray<complex<double> >&  d,
                                  valarray<complex<double> >& x) {
//-------------------------------------------------------------------
//   Fills solution into x. Warning: will modify c and d! 
//-------------------------------------------------------------------
    size_t total_loops(x.size()/calculations_per_loop);
    
    #pragma omp parallel for num_threads(Input::List().ompthreads)
    for (size_t iLoop = 0; iLoop < total_loops; ++iLoop)
    {
        size_t offset(iLoop*calculations_per_loop);
        // Modify the coefficients. 
        c[offset] /= b[offset];                            // Division by zero risk. 
        d[offset] /= b[offset];                            // Division by zero would imply a singular matrix. 
                                                               // 
        for (size_t i(offset+1); i < offset+calculations_per_loop; ++i)
        {
            double id(1.0/(b[i]-c[i-1]*a[i]));   // Division by zero risk. 
            c[i] *= id;                          // Last value calculated is redundant.
            d[i] -= d[i-1] * a[i];
            d[i] *= id;                          // d[i] = (d[i] - d[i-1] * a[i]) * id 
        }
     
        // Now back substitute. 
        x[offset+calculations_per_loop-1] = d[offset+calculations_per_loop-1];

        #pragma novector
        for (int i(offset+calculations_per_loop-2); i > offset-1; --i)
        {
            x[i]  = d[i];
            x[i] -= c[i] * x[i+1];               // x[i] = d[i] - c[i] * x[i + 1];
        }
    }
}
//-------------------------------------------------------------------
//*******************************************************************
//-------------------------------------------------------------------
     void TridiagonalSolve ( size_t calculations_per_loop,
                            valarray<double>& a, 
                            valarray<double>& b, 
                                  valarray<double>& c,      
                                  valarray<double>& d,
                                  valarray<double>& x) {
//-------------------------------------------------------------------
//   Fills solution into x. Warning: will modify c and d! 
//-------------------------------------------------------------------
    size_t total_loops(x.size()/calculations_per_loop);
    
    #pragma omp parallel for num_threads(Input::List().ompthreads)
    for (size_t iLoop = 0; iLoop < total_loops; ++iLoop)
    {
        size_t offset(iLoop*calculations_per_loop);
        // Modify the coefficients. 
        c[offset] /= b[offset];                            // Division by zero risk. 
        d[offset] /= b[offset];                            // Division by zero would imply a singular matrix. 
                                                               // 
        for (size_t i(offset+1); i < offset+calculations_per_loop; ++i)
        {
            double id(1.0/(b[i]-c[i-1]*a[i]));   // Division by zero risk. 
            c[i] *= id;                          // Last value calculated is redundant.
            d[i] -= d[i-1] * a[i];
            d[i] *= id;                          // d[i] = (d[i] - d[i-1] * a[i]) * id 
        }
     
        // Now back substitute. 
        x[offset+calculations_per_loop-1] = d[offset+calculations_per_loop-1];

        #pragma novector
        for (int i(offset+calculations_per_loop-2); i > offset-1; --i)
        {
            x[i]  = d[i];
            x[i] -= c[i] * x[i+1];               // x[i] = d[i] - c[i] * x[i + 1];
        }
    }
}
//*******************************************************************
//-------------------------------------------------------------------
void TridiagonalSolve ( valarray<double>& a,
                        valarray<double>& b,
                       valarray<double>& c,
                       valarray<double>& d,
                       valarray<double>& x) {
//-------------------------------------------------------------------
//   Fills solution into x. Warning: will modify c and d!
//-------------------------------------------------------------------
    // size_t n(d.size());
    // // // Modify the coefficients.
    // c[0] /= b[0];                            // Division by zero risk.
    // d[0] /= b[0];                            // Division by zero would imply a singular matrix.
    // for (size_t i(1); i < n; ++i){
    //     double id(1.0/(b[i]-c[i-1]*a[i]));   // Division by zero risk.
    //     c[i] *= id;	                         // Last value calculated is redundant.
    //     d[i] -= d[i-1] * a[i];
    //     d[i] *= id;                          // d[i] = (d[i] - d[i-1] * a[i]) * id
    // }

    // // Now back substitute.
    // x[n-1] = d[n-1];
    // for (int i(2); i < n+1; ++i)
    // {
    //     x[n-i]  = d[n-i];
    //     x[n-i] -= c[n-i] * x[n-i+1];               // x[i] = d[i] - c[i] * x[i + 1];
    // }

    int j, n(a.size());
    double bet(b[0]);
    valarray <double> gam(0.,n);

    if (b[0] == 0.) throw("Error 1 in TridiagonalSolve");

    x[0] = d[0]/b[0];

    for (j=1;j<n;++j)
    {
        gam[j] = c[j-1]/bet;
        bet = b[j]-a[j]*gam[j];
        if (bet == 0.) throw("Error 2 in TridiagonalSolve");   
        x[j] = (d[j]-a[j]*x[j-1])/bet;
    }
    #pragma novector
    for (j=n-2;j>=0;j--)
    {
        x[j] -= gam[j+1]*d[j+1];
    }
}
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//*******************************************************************
//-------------------------------------------------------------------
bool Thomas_Tridiagonal(Array2D<double>& A,
                        valarray<double> & d,
                        valarray<double> & xk) {
//-------------------------------------------------------------------
//   Fills solution into xk. The other matrices are not modified
//   The function returns "false" if the matrix A is not diagonally
//   dominant
//-------------------------------------------------------------------

//      The Matrices all have the right dimensions
//      -------------------------------------------------------------
    if ( ( A.dim1() != A.dim2()  ) ||
         ( A.dim1() != d.size()  ) ||
         ( A.dim1() != xk.size() )    )  {
        cout << "Error: The Matrices don't have the right dimensions!" << endl;
        exit(1);
    }
//      -------------------------------------------------------------

    valarray<double> a(d.size()), b(d.size()), c(d.size());

    for (size_t i(1); i < A.dim1(); ++i){
        a[i] = A(i,i-1);
        // std::cout << "\n a[" << i+1 << "] = " << a[i+1];
    }
    for (size_t i(0); i < A.dim1(); ++i){
        b[i] = A(i,i);
        // std::cout << "\n b[" << i << "] = " << b[i];
    }
    for (size_t i(0); i < A.dim1()-1; ++i){
        c[i] = A(i,i+1);
        // std::cout << "\n c[" << i << "] = " << c[i];
    }

    // for (size_t i(0); i < A.dim1(); ++i)
    // {
    //     std::cout << "\n d[" << i << "] = " << d[i];
    // }

    TridiagonalSolve(a,b,c,d,xk);


    // for (size_t i(0); i < A.dim1(); ++i)
    // {
    //     std::cout << "\n xk[" << i << "] = " << xk[i];
    // }

    // exit(1);



    return true;
}

//-------------------------------------------------------------------
bool Thomas_Tridiagonal(Array2D<double>& A,
                        valarray<complex<double> >& d,
                        valarray<complex<double> >& xk) {
//-------------------------------------------------------------------
//   Fills solution into xk. The other matrices are not modified
//   The function returns "false" if the matrix A is not diagonally
//   dominant
//-------------------------------------------------------------------

//      The Matrices all have the right dimensions
//      -------------------------------------------------------------
    if ( ( A.dim1() != A.dim2()  ) ||
         ( A.dim1() != d.size()  ) ||
         ( A.dim1() != xk.size() )    )  {
        cout << "Error: The Matrices don't have the right dimensions!" << endl;
        exit(1);
    }
//      -------------------------------------------------------------

    valarray<double> a(d.size()), b(d.size()), c(d.size());

    for (size_t i(0); i < A.dim1()-1; ++i){
        a[i+1] = A(i+1,i);
    }
    for (size_t i(0); i < A.dim1(); ++i){
        b[i] = A(i,i);
    }
    for (size_t i(0); i < A.dim1()-1; ++i){
        c[i] = A(i,i+1);
        // std::cout << "\n c[" << i << "] = " << c[i];
    }

//        valarray< double > dcopy(d);
    TridiagonalSolve(a,b,c,d,xk);
    // xk = d;
    return true;
}
//*******************************************************************
//-------------------------------------------------------------------
    complex <double> Det33(/*const valarray<double>& D, */
                          Array2D<complex <double> >& A) {           // Determinant for a 3*3 system
//-------------------------------------------------------------------
        return A(0,0) * ( A(1,1)*A(2,2) - A(2,1)*A(1,2) ) -
               A(1,0) * ( A(0,1)*A(2,2) - A(2,1)*A(0,2) ) +
               A(2,0) * ( A(0,1)*A(1,2) - A(1,1)*A(0,2) );
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    complex <double> Detx33(valarray<complex <double> >& D, 
                           Array2D<complex <double> >& A) {         // Determinant x for a 3*3 system
//-------------------------------------------------------------------
        return D[0] * ( A(1,1)*A(2,2) - A(2,1)*A(1,2) ) -
               D[1] * ( A(0,1)*A(2,2) - A(2,1)*A(0,2) ) +
               D[2] * ( A(0,1)*A(1,2) - A(1,1)*A(0,2) );
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    complex <double> Dety33(valarray<complex <double> >& D, 
                           Array2D<complex <double> >& A) {         // Determinant y for a 3*3 system
//-------------------------------------------------------------------
        return A(0,0) * ( D[1]*A(2,2) - D[2]*A(1,2) ) -
               A(1,0) * ( D[0]*A(2,2) - D[2]*A(0,2) ) +
               A(2,0) * ( D[0]*A(1,2) - D[1]*A(0,2) );
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    complex <double> Detz33(valarray<complex <double> >& D,
                           Array2D<complex <double> >& A) {         // Determinant z for a 3*3 system
//-------------------------------------------------------------------
        return A(0,0) * ( A(1,1)*D[2] - A(2,1)*D[1] ) -
               A(1,0) * ( A(0,1)*D[2] - A(2,1)*D[0] ) +
               A(2,0) * ( A(0,1)*D[1] - A(1,1)*D[0] );
    }
//-------------------------------------------------------------------


/**
 * @brief      Convert data structure to float structure
 *
 * @param[in]  vDouble  The v double
 *
 * @return     float structure
 */

vector<double> vdouble_real(const vector<complex<double> >& vDouble) {
    vector<double> vd;
    for (size_t i(0); i < vDouble.size(); ++i) {
        vd.push_back(static_cast<double>(vDouble[i].real()));
    }
    return vd;
}

vector<double> vdouble_imag(const vector<complex<double> >& vDouble) {
    vector<double> vd;
    for (size_t i(0); i < vDouble.size(); ++i) {
        vd.push_back(static_cast<double>(vDouble[i].imag()));
    }
    return vd;
}

vector<double> valtovec(const valarray<double>& vDouble) {
    vector<double> vf;
    for (size_t i(0); i < vDouble.size(); ++i) {
        vf.push_back(static_cast<double>(vDouble[i]));
    }
    return vf;
}

