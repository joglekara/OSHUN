/*! \brief Functors for various methods
 * \author PICKSC
 * \file   functors.cpp
 *
 * Includes spatial advection, electric field advection, and electric field update routines
 *
 */
//--------------------------------------------------------------
//  Standard libraries

#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <algorithm>
#include <cstdlib>
#include <omp.h>
#include <math.h>
#include <map>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"

//  Declerations
#include "state.h"
#include "setup.h"
#include "input.h"
#include "fluid.h"
#include "vlasov.h"
#include "functors.h"

//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods with explicit
//  E-field solver

//--------------------------------------------------------------
//  Constructor
VlasovFunctor1D_explicitE::VlasovFunctor1D_explicitE(vector<size_t> Nl,vector<size_t> Nm,vector<valarray<double> > dp,
                                                     double xmin, double xmax, size_t Nx):WD(xmin, xmax, Nx, 0., 1., 1)
{
//--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        SA.push_back( Spatial_Advection(Nl[s], Nm[s], dp[s], xmin, xmax, Nx, 0., 1., 1) );

        EF.push_back( Electric_Field(Nl[s], Nm[s], dp[s]) );        

        JX.push_back( Current() );

        BF.push_back( Magnetic_Field(Nl[s], Nm[s], dp[s]) );
    }

    AM.push_back( Ampere(xmin, xmax, Nx, 0., 1., 1) );

    FA.push_back( Faraday(xmin, xmax, Nx, 0., 1., 1) );

    
}
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Collect all of the terms
void VlasovFunctor1D_explicitE::operator()(const State1D& Yin, State1D& Yslope){
//--------------------------------------------------------------
    bool debug(0);

    Yslope = static_cast<complex<double> > (0.0);

    for (size_t s(0); s < Yin.Species(); ++s) 
    {

        if (Yin.DF(s).m0() == 0) {

            // GA[s].es1d(Yin.DF(s),Yslope.EMF().Ex());
            
            // if (debug) 
            // {
            //     std::cout << "\n\n f at start:";
            //     for (size_t ip(0); ip < Yin.SH(0,0,0).nump(); ++ip){
            //         std::cout << "\nf(" << ip << ") = " << Yin.SH(0,1,0)(ip,4);
            //     }

            //     std::cout << "\n\n E at start:";
            //     for (size_t ix(0); ix < Yin.SH(0,0,0).numx(); ++ix){
            //         std::cout << "\nEx(" << ix << ") = " << Yin.EMF().Ex()(ix);
            //     }            
            // }
            

            // if (Input::List().flm_acc)
            // {
                // EF[s].gpu1d(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));
                // SA[s].gpu1d(Yin.DF(s),Yslope.DF(s));
            // }
            if (Input::List().vgradf)
                SA[s].es1d(Yin.DF(s),Yslope.DF(s));
            // { 
            
            if (Input::List().Edfdv)
                EF[s].es1d(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));
                
            // }
            
            // if (debug) 
            // {
            //     std::cout << "\n\nf after E:";
            //     for (size_t ip(0); ip < Yin.SH(0,0,0).nump(); ++ip){
            //         std::cout << "\nf(" << ip << ") = " << Yslope.SH(0,1,0)(ip,4);
            //     }
            // }
            if (Input::List().dEdt)
                JX[s].es1d(Yin.DF(s),Yslope.EMF().Ex());



            // if (debug) 
            // {
            //     std::cout << "\n\n after J:";
            //     for (size_t ix(0); ix < Yin.SH(0,0,0).numx(); ++ix){
            //         std::cout << "\nEx(" << ix << ") = " << Yslope.EMF().Ex()(ix);
            //     }            
            // }
            
            

            // if (debug) 
            // {
            //     std::cout << "\n\n after SA:";
            //     for (size_t ip(0); ip < Yin.SH(0,0,0).nump(); ++ip){
            //         std::cout << "\nf(" << ip << ") = " << Yslope.SH(0,1,0)(ip,4);
            //     }
            // }
            // 
            // collide.advance(Yin,Yslope,time,dt);
        }
        else
        {
            if (Yin.DF(s).l0() == 1) 
            {
                SA[s].f1only(Yin.DF(s),Yslope.DF(s));

                EF[s].f1only(Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

                BF[s].f1only(Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

                JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

            }

            else 
            {
                SA[s](Yin.DF(s),Yslope.DF(s));

                EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

                BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

                JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());
            }
        }
    }
    
    AM[0](Yin.EMF(),Yslope.EMF());
    FA[0](Yin.EMF(),Yslope.EMF());
            
}

void VlasovFunctor1D_explicitE::operator()(const State1D& Yin, State1D& Yslope, size_t direction){}
void VlasovFunctor1D_explicitE::operator()(const State1D& Yin, State1D& Yslope, double time, double dt){
    bool debug(0);

    Yslope = static_cast<complex<double> > (0.0);

    

    for (size_t s(0); s < Yin.Species(); ++s) 
    {

        EMF1D EMF_ext(Yin.DF(s)(0,0).numx()); EMF_ext = static_cast<complex<double> > (0.0);
        if (Input::List().trav_wave) WD.applytravelingwave(EMF_ext,time + dt*0.5);

        if (Yin.DF(s).m0() == 0) 
        {
            // std::cout << "\n10\n";
            // EF[s].es1d(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));
            
            if (Input::List().dEdt)
                JX[s].es1d(Yin.DF(s),Yslope.EMF().Ex());
            // SA[s].es1d(Yin.DF(s),Yslope.DF(s));

            size_t l0(Yin.DF(s).l0());
            
            #pragma omp parallel num_threads(Input::List().ompthreads)
            {   
                size_t this_thread  = omp_get_thread_num();
                // std::cout << "\n hi i'm " << this_thread << "\n";

                size_t f_start_thread(EF[s].get_f_start(this_thread));
                size_t f_end_thread(EF[s].get_f_end(this_thread));

                
                // std::cout << "\n els[ " << this_thread << "] = " << f_start_thread << "....\n";
                // std::cout << "\n ele[ " << this_thread << "] = " << f_end_thread << "....\n";
                //  Initialize work variables
                SHarmonic1D fd1(SA[s].get_vr().size(),Yin.DF(s)(0,0).numx()),fd2(SA[s].get_vr().size(),Yin.DF(s)(0,0).numx());
                
                valarray<complex<double> > vtemp(SA[s].get_vr());
                vtemp /= Yin.DF(s).mass();


                valarray<complex<double> > Ex(Yin.FLD(0).array());
                
                // valarray<complex<double> > Ex_ext(EMF_ext(0).array());
                Ex = Ex + EMF_ext(0).array();
                Ex *= Yin.DF(s).q();

                if (!Input::List().Edfdv)
                    Ex = 0.;

                //  -------------------------------------------------------- //
                //   First thread takes the boundary conditions (l = 0)
                //   Last thread takes the boundary condition (l = l0)
                //  Rest proceed to chunks
                //  -------------------------------------------------------- //
                if (this_thread == 0)
                {
                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    //      m = 0, l = 0
                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    EF[s].MakeG00(Yin.DF(s)(0,0),fd1);
                    Ex *= EF[s].getA1(0,0);  Yslope.DF(s)(1,0) += fd1.mxaxis(Ex);


                    fd1 = Yin.DF(s)(0,0);                           fd1.Dx(Input::List().dbydx_order); //fd1 = complex<double>(0.);
                    vtemp *= SA[s].getA1(0,0);                      Yslope.DF(s)(1,0) += fd1.mpaxis(vtemp);
                    vtemp /= SA[s].getA1(0,0);

                    f_start_thread = 1;
                }

                if (this_thread == Input::List().ompthreads - 1)    
                {
                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    //      m = 0,  l = l0
                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    EF[s].MakeGH(Yin.DF(s)(l0,0),fd1,fd2,l0);
                    Ex *= EF[s].getA2(l0,0);  Yslope.DF(s)(l0-1,0) += fd2.mxaxis(Ex);
                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Ex /= EF[s].getA2(l0,0);             // Reset Ex
                    
                        
                    fd1 = Yin.DF(s)(l0,0);                          fd1.Dx(Input::List().dbydx_order); //fd1 = complex<double>(0.);
                    vtemp *= SA[s].getA2(l0,0);                     Yslope.DF(s)(l0-1,0) += fd1.mpaxis(vtemp);
                    vtemp /= SA[s].getA2(l0,0);

                    f_end_thread -= 1;
                }

                //  -------------------------------------------------------- //
                //  Do the chunks
                //  Initialize vtemp so that it starts correctly
                //  -------------------------------------------------------- //
                vtemp *= SA[s].getA1(f_start_thread-1,0);
                Ex *= EF[s].getA1(f_start_thread-1,0);


                for (size_t l = f_start_thread; l < f_end_thread; ++l)
                {
                    EF[s].MakeGH(Yin.DF(s)(l,0),fd1,fd2,l);

                    Ex *= EF[s].getA2(l,0) / EF[s].getA1(l-1,0);    Yslope.DF(s)(l-1,0) += fd2.mxaxis(Ex);
                    Ex *= EF[s].getA1(l,0) / EF[s].getA2(l,0);      Yslope.DF(s)(l+1,0) += fd1.mxaxis(Ex);

                    fd1 = Yin.DF(s)(l,0);  //std::cout << "\n \n before dx, l = " << l << " \n";          
                    fd1.Dx(Input::List().dbydx_order);  //fd1 = complex<double>(0.); //std::cout << " \n after dx\n";

                    vtemp *= SA[s].getA2(l,0)/SA[s].getA1(l-1,0);    fd2 = fd1;  Yslope.DF(s)(l-1,0) += fd1.mpaxis(vtemp);
                    vtemp *= SA[s].getA1(l,0)/SA[s].getA2(l  ,0);                Yslope.DF(s)(l+1,0) += fd2.mpaxis(vtemp);
                }
            }

            //  -------------------------------------------------------- //
            //  Do the boundaries between the chunks
            //  -------------------------------------------------------- //
            #pragma omp parallel for num_threads(Input::List().ompthreads-1)
            for (size_t threadboundaries = 0; threadboundaries < Input::List().ompthreads - 1; ++threadboundaries)
            {
                SHarmonic1D fd1(Yin.DF(s)(0,0).nump(),Yin.DF(s)(0,0).numx()),fd2(Yin.DF(s)(0,0).nump(),Yin.DF(s)(0,0).numx());
                valarray<complex<double> > vtemp(SA[s].get_vr());
                vtemp /= Yin.DF(s).mass();    

                valarray<complex<double> > Ex(Yin.FLD(0).array());

                Ex = Ex + EMF_ext(0).array();

                Ex *= Yin.DF(s).q();

                if (!Input::List().Edfdv)
                    Ex = 0.;

            //  Initialize Ex so that it its ready for loop iteration l
                Ex *= EF[s].getA1(EF[s].get_f_end(threadboundaries)-1,0);    
                
                vtemp *= SA[s].getA1(EF[s].get_f_end(threadboundaries)-1,0);

                for (size_t l = EF[s].get_f_end(threadboundaries); l < EF[s].get_f_start(threadboundaries+1); ++l)
                {   
                    EF[s].MakeGH(Yin.DF(s)(l,0),fd1,fd2,l);

                    Ex *= EF[s].getA2(l,0) / EF[s].getA1(l-1,0);     Yslope.DF(s)(l-1,0) += fd2.mxaxis(Ex);
                    Ex *= EF[s].getA1(l,0) / EF[s].getA2(l,0);       Yslope.DF(s)(l+1,0) += fd1.mxaxis(Ex); 

                    fd1 = Yin.DF(s)(l,0);  //std::cout << "\n \n before dx, l = " << l << " \n";          
                    fd1.Dx(Input::List().dbydx_order); //fd1 = complex<double>(0.);//std::cout << " \n after dx\n";

                    vtemp *= SA[s].getA2(l,0)/SA[s].getA1(l-1,0);    fd2 = fd1;  Yslope.DF(s)(l-1,0) += fd1.mpaxis(vtemp);
                    vtemp *= SA[s].getA1(l,0)/SA[s].getA2(l  ,0);                Yslope.DF(s)(l+1,0) += fd2.mpaxis(vtemp);
                }
            }         

        }
        else
        {
            if (Yin.DF(s).l0() == 1) 
            {
                SA[s].f1only(Yin.DF(s),Yslope.DF(s));

                EF[s].f1only(Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

                BF[s].f1only(Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

                JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

            }

            else 
            {
                SA[s](Yin.DF(s),Yslope.DF(s));

                EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

                BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

                JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());
            }

            AM[0](Yin.EMF(),Yslope.EMF());
            FA[0](Yin.EMF(),Yslope.EMF());
        }
        
        if (Input::List().filterdistribution)  Yslope.DF(s).Filterp();

    }

    // if (Input::List().trav_wave) WD.applytravelingwave(Yslope.EMF(),time);
    
}

// //**************************************************************
// //--------------------------------------------------------------
// //  Functor to be used in the Runge-Kutta methods with explicit
// //  E-field solver

// //--------------------------------------------------------------
// //  Constructor
// VlasovFunctor1D_spatialAdvection::VlasovFunctor1D_spatialAdvection(vector<size_t> Nl,vector<size_t> Nm,
//                                                     // vector<double> pmax, vector<size_t> Np,
//                                                     vector<valarray<double> > dp,
//                                                      double xmin, double xmax, size_t Nx) {
// //--------------------------------------------------------------

//     for (size_t s(0); s < Nl.size(); ++s){

//         SA.push_back( Spatial_Advection(Nl[s], Nm[s], dp[s], xmin, xmax, Nx, 0., 1., 1) );

//     }
// }
// //--------------------------------------------------------------


// //--------------------------------------------------------------
// //  Collect all of the terms
// void VlasovFunctor1D_spatialAdvection::operator()(const State1D& Yin, State1D& Yslope){
// //--------------------------------------------------------------
//     bool debug(0);

//     Yslope = 0.0;

//     for (size_t s(0); s < Yin.Species(); ++s) {

//         if (Yin.DF(s).l0() == 1) {

//             SA[s].f1only(Yin.DF(s),Yslope.DF(s));

//         }
//         else if (Yin.DF(s).m0() == 0) {
            
//             SA[s].es1d(Yin.DF(s),Yslope.DF(s));

//         }

//         else {

//             SA[s](Yin.DF(s),Yslope.DF(s));            
//         }

//     }

// }

// void VlasovFunctor1D_spatialAdvection::operator()(const State1D& Yin, State1D& Yslope, size_t direction){}

// //--------------------------------------------------------------
// //  Constructor
// VlasovFunctor1D_fieldUpdate::VlasovFunctor1D_fieldUpdate(vector<size_t> Nl,vector<size_t> Nm,
//                                                     // vector<double> pmax, vector<size_t> Np,
//                                                     vector<valarray<double> > dp,
//                                                      double xmin, double xmax, size_t Nx) {
// //--------------------------------------------------------------

//     for (size_t s(0); s < Nl.size(); ++s)
//     {

//         JX.push_back( Current() );

//         AM.push_back( Ampere(xmin, xmax, Nx, 0., 1., 1) );

//         FA.push_back( Faraday(xmin, xmax, Nx, 0., 1., 1) );

//     }
// }
// //--------------------------------------------------------------


// //--------------------------------------------------------------
// //  Collect all of the terms
// void VlasovFunctor1D_fieldUpdate::operator()(const State1D& Yin, State1D& Yslope){
// //--------------------------------------------------------------
//     bool debug(0);

//     Yslope = 0.0;

//     for (size_t s(0); s < Yin.Species(); ++s) {

//         if (Yin.DF(s).l0() == 1) {

//             JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

//             AM[s](Yin.EMF(),Yslope.EMF());

//             FA[s](Yin.EMF(),Yslope.EMF());

//         }
//         else if (Yin.DF(s).m0() == 0) {

//             JX[s].es1d(Yin.DF(s),Yslope.EMF().Ex());
//         }

//         else {

//             JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

//             AM[s](Yin.EMF(),Yslope.EMF());

//             FA[s](Yin.EMF(),Yslope.EMF());            
//         }

//     }

// }

// void VlasovFunctor1D_fieldUpdate::operator()(const State1D& Yin, State1D& Yslope, size_t direction){}

// //  Constructor
// VlasovFunctor1D_momentumAdvection::VlasovFunctor1D_momentumAdvection(vector<size_t> Nl,vector<size_t> Nm,
//                                                     // vector<double> pmax, vector<size_t> Np,
//                                                     vector<valarray<double> > dp,
//                                                      double xmin, double xmax, size_t Nx) {
// //--------------------------------------------------------------

//     for (size_t s(0); s < Nl.size(); ++s){

//         EF.push_back( Electric_Field(Nl[s], Nm[s], dp[s]) );        

//         BF.push_back( Magnetic_Field(Nl[s], Nm[s], dp[s]) );

//     }
// }
// //--------------------------------------------------------------


// //--------------------------------------------------------------
// //  Collect all of the terms
// void VlasovFunctor1D_momentumAdvection::operator()(const State1D& Yin, State1D& Yslope){
// //--------------------------------------------------------------
//     bool debug(0);

//     Yslope = 0.0;

//     for (size_t s(0); s < Yin.Species(); ++s) {

//         if (Yin.DF(s).l0() == 1) {

//             EF[s].f1only(Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

//             BF[s].f1only(Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

//         }
//         else if (Yin.DF(s).m0() == 0) {
            
//             EF[s].es1d(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));

//         }

//         else {

//             EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

//             BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));
                     
//         }

//     }

// }

// void VlasovFunctor1D_momentumAdvection::operator()(const State1D& Yin, State1D& Yslope, size_t direction){}


//--------------------------------------------------------------
//  Constructor
VlasovFunctor2D_explicitE::VlasovFunctor2D_explicitE(vector<size_t> Nl,vector<size_t> Nm,
                                                    // vector<double> pmax, vector<size_t> Np,
                                                    vector<valarray<double> > dp,
                                                     double xmin, double xmax, size_t Nx,
                                                     double ymin, double ymax, size_t Ny):
                                                     WD(xmin, xmax, Nx, ymin, ymax, Ny) {
//--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        SA.push_back( Spatial_Advection(Nl[s], Nm[s], dp[s], xmin, xmax, Nx, ymin, ymax, Ny) );

        EF.push_back( Electric_Field(Nl[s], Nm[s], dp[s]) );        

        JX.push_back( Current() );

        BF.push_back( Magnetic_Field(Nl[s], Nm[s], dp[s]) );    

    }
    AM.push_back( Ampere(xmin, xmax, Nx, ymin, ymax, Ny) );

    FA.push_back( Faraday(xmin, xmax, Nx, ymin, ymax, Ny) );
}
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Collect all of the terms
void VlasovFunctor2D_explicitE::operator()(const State2D& Yin, State2D& Yslope){
//--------------------------------------------------------------
    bool debug(0);

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) 
    {

        if (Yin.DF(s).l0() == 1) 
        {

            SA[s].f1only(Yin.DF(s),Yslope.DF(s));

            EF[s].f1only(Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

            BF[s].f1only(Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

            JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

        }
        
        else 
        {
            // EMF2D EMF_ext(Yin.DF(s)(0,0).numx(),Yin.DF(s)(0,0).numy()); EMF_ext = static_cast<complex<double> > (0.0);
            // if (Input::List().trav_wave) WD.applytravelingwave(EMF_ext,time + dt*0.5);

            SA[s](Yin.DF(s),Yslope.DF(s));

            EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

            BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

            JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());
            
        }

    }

    AM[0](Yin.EMF(),Yslope.EMF());

    FA[0](Yin.EMF(),Yslope.EMF());

}

void VlasovFunctor2D_explicitE::operator()(const State2D& Yin, State2D& Yslope, size_t direction){}
void VlasovFunctor2D_explicitE::operator()(const State2D& Yin, State2D& Yslope, double time, double dt){
    bool debug(0);

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) 
    {

        if (Yin.DF(s).l0() == 1) 
        {

            SA[s].f1only(Yin.DF(s),Yslope.DF(s));

            EF[s].f1only(Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

            BF[s].f1only(Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

            JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

        }
        
        else 
        {
            // EMF2D EMF_ext(Yin.DF(s)(0,0).numx(),Yin.DF(s)(0,0).numy()); EMF_ext = static_cast<complex<double> > (0.0);
            if (Input::List().trav_wave) WD.applytravelingwave(Yin.EMF(),time + dt*0.5);

            SA[s](Yin.DF(s),Yslope.DF(s));

            EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

            BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

            JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());
            
        }

    }

    AM[0](Yin.EMF(),Yslope.EMF());

    FA[0](Yin.EMF(),Yslope.EMF());

}
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods with implicit
//  E-field solver

//--------------------------------------------------------------
//  Constructor
VlasovFunctor1D_implicitE_p1::VlasovFunctor1D_implicitE_p1(vector<size_t> Nl,vector<size_t> Nm,
                                                            // vector<double> pmax, vector<size_t> Np,
                                                    vector<valarray<double> > dp,
                                                           double xmin, double xmax, size_t Nx) {
// //--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        SA.push_back( Spatial_Advection(Nl[s], Nm[s], dp[s], xmin, xmax, Nx, 0., 1., 1) );

        // HA.push_back( Hydro_Advection_1D(Nl[s], Nm[s], dp[s], xmin, xmax, Nx) );

        BF.push_back( Magnetic_Field(Nl[s], Nm[s], dp[s]) );

    }

}

//--------------------------------------------------------------
//
//
void VlasovFunctor1D_implicitE_p1::operator()(const State1D& Yin, State1D& Yslope){
//--------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {

        if (Yin.DF(s).l0() == 1) 
        {
            SA[s].f1only(Yin.DF(s),Yslope.DF(s));
            BF[s].f1only(Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));
            // HA[s](Yin.DF(s),Yin.HYDRO(),Yslope.DF(s));
        }
        else {
            SA[s](Yin.DF(s),Yslope.DF(s));
            BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));
            // HA[s](Yin.DF(s),Yin.HYDRO(),Yslope.DF(s));
        }

        // if (Input::List().filterdistribution)  Yslope.DF(s) = Yslope.DF(s).Filterp();


    }

}
void VlasovFunctor1D_implicitE_p1::operator()(const State1D& Yin, State1D& Yslope, size_t direction){}
void VlasovFunctor1D_implicitE_p1::operator()(const State1D& Yin, State1D& Yslope, double time, double dt){}

//--------------------------------------------------------------
//  Constructor
VlasovFunctor2D_implicitE_p1::VlasovFunctor2D_implicitE_p1(vector<size_t> Nl,vector<size_t> Nm,
                                            
                                                    vector<valarray<double> > dp,
                                                    double xmin, double xmax, size_t Nx,
                                                    double ymin, double ymax, size_t Ny) {
// //--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        SA.push_back( Spatial_Advection(Nl[s], Nm[s], dp[s], xmin, xmax, Nx, ymin, ymax, Ny) );

        // HA.push_back( Hydro_Advection_1D(Nl[s], Nm[s], dp[s], xmin, xmax, Nx) );

        BF.push_back( Magnetic_Field(Nl[s], Nm[s], dp[s]) );

    }

}

//--------------------------------------------------------------
//
//
void VlasovFunctor2D_implicitE_p1::operator()(const State2D& Yin, State2D& Yslope){
//--------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {

        if (Yin.DF(s).l0() == 1) 
        {
            SA[s].f1only(Yin.DF(s),Yslope.DF(s));
            BF[s].f1only(Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));
            // HA[s](Yin.DF(s),Yin.HYDRO(),Yslope.DF(s));
        }
        else {
            SA[s](Yin.DF(s),Yslope.DF(s));
            BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));
            // HA[s](Yin.DF(s),Yin.HYDRO(),Yslope.DF(s));
        }

        // if (Input::List().filterdistribution)  Yslope.DF(s) = Yslope.DF(s).Filterp();


    }

}
void VlasovFunctor2D_implicitE_p1::operator()(const State2D& Yin, State2D& Yslope, size_t direction){}
void VlasovFunctor2D_implicitE_p1::operator()(const State2D& Yin, State2D& Yslope, double time, double dt){}


//--------------------------------------------------------------
//  Constructor
VlasovFunctor1D_implicitE_p2::VlasovFunctor1D_implicitE_p2(vector<size_t> Nl,vector<size_t> Nm,
                                                                vector<valarray<double> > dp,
                                                           double xmin, double xmax, size_t Nx) {
// //--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        FA.push_back( Faraday(xmin, xmax, Nx, 0., 1., 1) );
        EF.push_back( Electric_Field(Nl[s], Nm[s], dp[s]) );

    }

}
//------------------------------------------------------------------------------------------------------
void VlasovFunctor1D_implicitE_p2::operator()(const State1D& Yin, State1D& Yslope){
// //---------------------------------------------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {
        FA[s](Yin.EMF(),Yslope.EMF());

        if (Yin.DF(s).l0() == 1) EF[s].f1only(Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));
        else                     EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

        // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();

    }

}

//--------------------------------------------------------------------------------------------------
void VlasovFunctor1D_implicitE_p2::operator()(const State1D& Yin, State1D& Yslope, size_t direction){
//--------------------------------------------------------------------------------------------------

    Yslope = 0.0;

    if (direction == 1)
    {
        for (size_t s(0); s < Yin.Species(); ++s) {
            if (Yin.DF(s).l0() == 1) EF[s].Implicit_Ex_f1only(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));
            else                     EF[s].Implicit_Ex(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));
            // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }
    else if (direction == 2)
    {
        for (size_t s(0); s < Yin.Species(); ++s) {
            if (Yin.DF(s).l0() == 1) EF[s].Implicit_Ey_f1only(Yin.DF(s),Yin.EMF().Ey(),Yslope.DF(s));
            else                     EF[s].Implicit_Ey(Yin.DF(s),Yin.EMF().Ey(),Yslope.DF(s));
            // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }
    else
    {
        // complex<double> ii(0.0,1.0);
        for (size_t s(0); s < Yin.Species(); ++s) {
            if (Yin.DF(s).l0() == 1) EF[s].Implicit_Ez_f1only(Yin.DF(s),Yin.EMF().Ez(),Yslope.DF(s));
            else                     EF[s].Implicit_Ez(Yin.DF(s),Yin.EMF().Ez(),Yslope.DF(s));
            // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }

}
void VlasovFunctor1D_implicitE_p2::operator()(const State1D& Yin, State1D& Yslope, double time, double dt){}
//--------------------------------------------------------------
//  Constructor
VlasovFunctor2D_implicitE_p2::VlasovFunctor2D_implicitE_p2(vector<size_t> Nl,vector<size_t> Nm,
                                                                vector<valarray<double> > dp,
                                                           double xmin, double xmax, size_t Nx,
                                                           double ymin, double ymax, size_t Ny) {
// //--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        FA.push_back( Faraday(xmin, xmax, Nx, ymin, ymax, Ny) );
        EF.push_back( Electric_Field(Nl[s], Nm[s], dp[s]) );

    }

}
//------------------------------------------------------------------------------------------------------
void VlasovFunctor2D_implicitE_p2::operator()(const State2D& Yin, State2D& Yslope){
// //---------------------------------------------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {
        FA[s](Yin.EMF(),Yslope.EMF());

        if (Yin.DF(s).l0() == 1) EF[s].f1only(Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));
        else                     EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

        // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();

    }

}

//--------------------------------------------------------------------------------------------------
void VlasovFunctor2D_implicitE_p2::operator()(const State2D& Yin, State2D& Yslope, size_t direction){
//--------------------------------------------------------------------------------------------------

    Yslope = 0.0;

    if (direction == 1)
    {
        for (size_t s(0); s < Yin.Species(); ++s) {
            if (Yin.DF(s).l0() == 1) EF[s].Implicit_Ex_f1only(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));
            else                     EF[s].Implicit_Ex(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));
            // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }
    else if (direction == 2)
    {
        for (size_t s(0); s < Yin.Species(); ++s) {
            if (Yin.DF(s).l0() == 1) EF[s].Implicit_Ey_f1only(Yin.DF(s),Yin.EMF().Ey(),Yslope.DF(s));
            else                     EF[s].Implicit_Ey(Yin.DF(s),Yin.EMF().Ey(),Yslope.DF(s));
            // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }
    else
    {
        // complex<double> ii(0.0,1.0);
        for (size_t s(0); s < Yin.Species(); ++s) {
            if (Yin.DF(s).l0() == 1) EF[s].Implicit_Ez_f1only(Yin.DF(s),Yin.EMF().Ez(),Yslope.DF(s));
            else                     EF[s].Implicit_Ez(Yin.DF(s),Yin.EMF().Ez(),Yslope.DF(s));
            // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }

}
void VlasovFunctor2D_implicitE_p2::operator()(const State2D& Yin, State2D& Yslope, double time, double dt){}


