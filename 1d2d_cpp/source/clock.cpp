/*! \brief Clock manager
 * \author PICKSC
 * \file   clock.cpp
 *
 * Includes spatial advection, electric field advection, and electric field update routines
 *
 */
//--------------------------------------------------------------
//  Standard libraries
#include <mpi.h>
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <algorithm>
#include <cstdlib>
#include <mpi.h>

#include <math.h>
#include <map>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"

//  Declerations
#include "state.h"
#include "input.h"
#include "formulary.h"
#include "collisions.h"
#include "vlasov.h"
#include "functors.h"
#include "parallel.h"
#include "stepper.h"
#include "clock.h"

//**************************************************************

//**************************************************************
//--------------------------------------------------------------
Clock::Clock(double starttime, double __dt, double abs_tol, double rel_tol, size_t _maxfails,
                    State1D& Y): 
    current_time(starttime), dt_next(0.5*__dt), _dt(0.5*__dt),
    atol(abs_tol), rtol(rel_tol), 
    acceptability(0.), err_val(0.), 
    failed_steps(0), max_failures(_maxfails), _success(false),
    Nbc(Input::List().BoundaryCells), world_rank(0), world_size(1),
    Solver(Y)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        acceptabilitylist = new double[world_size];



    }
//--------------------------------------------------------------
Clock:: ~Clock(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    delete[] acceptabilitylist;
}
//--------------------------------------------------------------
Clock& Clock::operator++() 
{

    // std::cout << "\n time = "  << current_time;
    // end_of_loop_output();
    // end_of_loop_dist_output();
    // end_of_loop_dist_output();
    end_of_loop_time_updates();    
    
    return *this;
}
//--------------------------------------------------------------
//  Collect all of the terms
void Clock::end_of_loop_time_updates()
{
    dt_next = min(dt_next,1.01*_dt);
    dt_next = min(dt_next,Input::List().dt);
    failed_steps = 0;
    current_time += _dt;
    _dt = dt_next; 
    _success = 0;
}
//--------------------------------------------------------------
//  Collect all of the terms
void Clock::do_step(State1D& Ystar, State1D& Y_new, State1D& Y_old,
                        VlasovFunctor1D_explicitE& vF, collisions_1D& cF, Parallel_Environment_1D& PE)
{
    if (Input::List().adaptive_dt)
    {    
        while (_success == 0)
        {
            Solver.take_step(Ystar, Y_new, current_time, _dt, vF, cF, PE);
            update_dt(Y_old, Ystar, Y_new);
        }
        cF.advance(Y_new,current_time,_dt);
        PE.Neighbor_Communications(Y_new);
    }
    
    else Solver.take_step(Ystar, Y_new, current_time, _dt, vF, cF, PE);
}
//--------------------------------------------------------------
//  Collect all of the terms
void Clock::update_dt(State1D& Y_old, const State1D& Ystar, State1D& Y_new){
//--------------------------------------------------------------

    // err_val =  check_temperature(Ystar,Y);  

    /// Ex is checked for convergence
    err_val =  check_flds(Ystar,Y_new,acceptability);
    // err_val =  check_last_harmonic(Ystar,Y_new,acceptability);

    /// Error is determined
    acceptability = err_val/(atol + rtol*acceptability);
    if (acceptability > 1) _success = 0;
    else _success = 1;

    // std::cout << "\n acc = " << acceptability;
    std::cout << "\n dt = " << _dt;

    /// Error is shared so that global timestep can be determined on rank 0
    MPI_Gather(&acceptability, 1, MPI_DOUBLE, acceptabilitylist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (world_rank == 0)
    {   

        /// Determine whether to proceed
        for (size_t iprocess(0); iprocess < world_size; ++iprocess)
        {
            _success = int(((!(acceptabilitylist[iprocess] > 1)) && _success));
            acceptability = max(acceptabilitylist[iprocess],acceptability);
        }

        /// Update timestep and if failed, add an iteration
        if (_success == 1)
        {
            dt_next = 0.9*_dt/pow(acceptability,0.25);   
        }
        else
        {
            dt_next = 0.9*_dt/pow(acceptability,1./3.);
            
            ++failed_steps;
            if (failed_steps > max_failures) 
            {
                fprintf(stderr, "Solution failed to converge within %d steps \n", max_failures);
                MPI_Finalize();
                exit(1);
            }
        }
    }

    /// Share success and new timestep
    MPI_Bcast(&_success, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt_next, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /// If failed, restore old state and updated time step.
    /// Success time-step is updated at the end of outer loop.
    if (_success == 0)
    {
        _dt = dt_next;
        Y_new = Y_old;
    }
    else
    {
        Y_old = Y_new;
    }
}
//--------------------------------------------------------------
//  Collect all of the terms
void Clock::update_dt(const State2D& Y_old, const State2D& Ystar, State2D& Y_new){
//--------------------------------------------------------------

    // err_val =  check_temperature(Ystar,Y);  

    /// Ex is checked for convergence
    err_val =  check_flds(Ystar,Y_new,acceptability);
    // err_val =  check_last_harmonic(Ystar,Y_new,acceptability);

    /// Error is determined
    acceptability = err_val/(atol + rtol*acceptability);
    if (acceptability > 1) _success = 0;
    else _success = 1;

    // std::cout << "\n acc = " << acceptability;
    std::cout << "\n dt = " << _dt;

    /// Error is shared so that global timestep can be determined on rank 0
    MPI_Gather(&acceptability, 1, MPI_DOUBLE, acceptabilitylist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (world_rank == 0)
    {   

        /// Determine whether to proceed
        for (size_t iprocess(0); iprocess < world_size; ++iprocess)
        {
            _success = int(((!(acceptabilitylist[iprocess] > 1)) && _success));
            acceptability = max(acceptabilitylist[iprocess],acceptability);
        }

        /// Update timestep and if failed, add an iteration
        if (_success == 1)
        {
            dt_next = 0.9*_dt/pow(acceptability,0.25);   
        }
        else
        {
            dt_next = 0.9*_dt/pow(acceptability,1./3.);
            
            ++failed_steps;
            if (failed_steps > max_failures) 
            {
                fprintf(stderr, "Solution failed to converge within %d steps \n", max_failures);
                MPI_Finalize();
                exit(1);
            }
        }
    }

    /// Share success and new timestep
    MPI_Bcast(&_success, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt_next, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /// If failed, restore old state and updated time step.
    /// Success time-step is updated at the end of outer loop.
    if (_success == 0)
    {
        _dt = dt_next;
        // Y_new = Y_old;
    }
    else
    {
        // Y_old = Y_new;
    }
}

//--------------------------------------------------------------
//  Collect all of the terms
double Clock::check_temperature(const State1D& Ystar, const State1D& Y, double& maxval){
//--------------------------------------------------------------
    maxval = 0.;
    valarray<double> tmp_err(Y.FLD(0).numx());


    // Add up temperature over whole grid
    for (size_t s(0); s < Ystar.Species(); ++s)
    {
        tmp_err = Ystar.DF(s).getpressure() - Y.DF(s).getpressure();
        maxval = (Ystar.DF(s).getpressure()).max();
        
        tmp_err *= tmp_err;
        return sqrt(tmp_err.sum());
    }
    
}

//--------------------------------------------------------------
//  Collect all of the terms
double Clock::check_last_harmonic(const State1D& Ystar, const State1D& Y, double& maxval){
//--------------------------------------------------------------
    maxval = 0.;
    
    valarray<double> fave(Y.SH(0,0,0).nump());
    valarray<double> fstarave(Y.SH(0,0,0).nump());


    size_t nl(Ystar.DF(0).dim());
    --nl;

    
    for (size_t ix(0); ix < Y.FLD(0).numx(); ++ix)
    {
        for (size_t ip(0); ip < Y.SH(0,0,0).nump(); ++ip)
        {
            fave[ip] += Y.SH(0,nl,0)(ix,ip).real();
            fstarave[ip] += Ystar.SH(0,nl,0)(ix,ip).real();
        }
    }

    fave *= fave;
    fstarave *= fstarave;

    fstarave -= fave;
    return sqrt(fstarave.sum());
}

//--------------------------------------------------------------
//  Collect all of the terms
double Clock::check_flds(const State1D& Ystar, const State1D& Y, double& maxval){
//--------------------------------------------------------------
    maxval = 0.;
    valarray<double> tmp_err(Y.FLD(0).numx());

    for (size_t i(0); i < 1; ++i)
    {
        for (size_t ix(Nbc); ix < tmp_err.size() - Nbc; ++ix)
        {
            tmp_err[ix] = Ystar.FLD(i)(ix).real()-Y.FLD(i)(ix).real();
            maxval = max(maxval,abs(Y.FLD(i)(ix).real()));
            maxval = max(maxval,abs(Ystar.FLD(i)(ix).real()));
        }

        tmp_err *= tmp_err;
        return sqrt(tmp_err.sum());
    }
}
//--------------------------------------------------------------
//  Collect all of the terms
double Clock::check_flds(const State2D& Ystar, const State2D& Y, double& maxval){
//--------------------------------------------------------------
    maxval = 0.;
    valarray<double> tmp_err(Y.FLD(0).numx()*Y.FLD(0).numy());
    size_t ixy;
    for (size_t i(0); i < 1; ++i)
    {
        ixy = 0;
        for (size_t ix(Nbc); ix < Y.FLD(i).numx() - Nbc; ++ix)
        {
            for (size_t iy(Nbc); iy < Y.FLD(i).numy() - Nbc; ++iy)
            {
                tmp_err[ixy] = Ystar.FLD(i)(ix,iy).real()-Y.FLD(i)(ix,iy).real();
                maxval = max(maxval,abs(Y.FLD(i)(ix,iy).real()));
                maxval = max(maxval,abs(Ystar.FLD(i)(ix,iy).real()));
                ++ixy;
            }
        }

        tmp_err = abs(tmp_err);

        return tmp_err.sum();
    }
}


//--------------------------------------------------------------
//  Collect all of the terms
double Clock::check_js(const State1D& Ystar, const State1D& Y){
//--------------------------------------------------------------
    
    // for (size_t s(0); s < Ystar.Species(); ++s)
    // {
    //     for (size_t i(0); i < 3; ++i)
    //     {
    //         tmp_err = Ystar.DF(s).getcurrent(i) - Y.DF(s).getcurrent(i);
    //         tmp_err *= tmp_err;

    //         // if (tmp_err.sum() < total_tol) return true;        
    //         return tmp_err.sum();
    //     }
    // }

}
