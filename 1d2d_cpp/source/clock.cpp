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

#include <math.h>
#include <map>
#include <iomanip>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"

#include "external/exprtk.hpp"
#include "external/highfive/H5DataSet.hpp"

//  Declarations
#include "state.h"
#include "input.h"
#include "formulary.h"
#include "collisions.h"
#include "vlasov.h"
#include "setup.h"
#include "functors.h"
#include "parallel.h"
#include "export.h"
#include "stepper.h"
#include "clock.h"

//**************************************************************

//**************************************************************
//--------------------------------------------------------------
Clock::Clock(double starttime, double __dt, double abs_tol, double rel_tol, size_t _maxfails,
                    State1D& Y): 
    current_time(starttime), dt_next(__dt), _dt(__dt),
    atol(abs_tol), rtol(rel_tol), 
    acceptability(0.), err_val(0.), 
    failed_steps(0), max_failures(_maxfails), _success(0),
    Nbc(Input::List().BoundaryCells), world_rank(0), world_size(1),
    Solver(Y)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        acceptabilitylist = new double[world_size];

        /// /// /// /// /// /// /// /// /// //
        // int tout_start;
        if (Input::List().isthisarestart) 
        {
            tout_start = Input::List().restart_time;
        }
        else
        {
            tout_start = 0;
        }
        t_out = tout_start+1;
        
        dt_out = Input::List().t_stop / (Input::List().n_outsteps);
        dt_dist_out = Input::List().t_stop / (Input::List().n_distoutsteps);
        dt_big_dist_out = Input::List().t_stop / (Input::List().n_bigdistoutsteps);
        dt_restart = Input::List().t_stop / (Input::List().n_restarts);
        
        next_out = t_out*dt_out;
        next_dist_out = (tout_start*dt_out)+dt_dist_out;
        next_big_dist_out = (tout_start*dt_out)+dt_big_dist_out;

        next_restart = (tout_start*dt_out)+dt_restart;

        start_time = tout_start*dt_out;

        timings_at_current_timestep.push_back(0.);      // Vlasov
        timings_at_current_timestep.push_back(0.);      // Fokker-Planck
        timings_at_current_timestep.push_back(0.);      // Output Routines
        timings_at_current_timestep.push_back(0.);      // Big Output Routines

        timing_indices.push_back(0.);
        timing_indices.push_back(1.);
        timing_indices.push_back(2.);
        timing_indices.push_back(3.);
    }
//--------------------------------------------------------------
Clock::Clock(double starttime, double __dt, double abs_tol, double rel_tol, size_t _maxfails,
                    State2D& Y): 
    current_time(starttime), dt_next(__dt), _dt(__dt),
    atol(abs_tol), rtol(rel_tol), 
    acceptability(0.), err_val(0.), 
    failed_steps(0), max_failures(_maxfails), _success(0),
    Nbc(Input::List().BoundaryCells), world_rank(0), world_size(1),
    Solver(Y)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        acceptabilitylist = new double[world_size];

        /// /// /// /// /// /// /// /// /// //
        // int tout_start;
        if (Input::List().isthisarestart) 
        {
            tout_start = Input::List().restart_time;
        }
        else
        {
            tout_start = 0;
        }
        t_out = tout_start+1;
        
        dt_out = Input::List().t_stop / (Input::List().n_outsteps);
        dt_dist_out = Input::List().t_stop / (Input::List().n_distoutsteps);
        dt_big_dist_out = Input::List().t_stop / (Input::List().n_bigdistoutsteps);
        dt_restart = Input::List().t_stop / (Input::List().n_restarts);
        
        next_out = t_out*dt_out;
        next_dist_out = (tout_start*dt_out)+dt_dist_out;
        next_big_dist_out = (tout_start*dt_out)+dt_big_dist_out;

        next_restart = (tout_start*dt_out)+dt_restart;

        start_time = tout_start*dt_out;

        timings_at_current_timestep.push_back(0.);      // Vlasov
        timings_at_current_timestep.push_back(0.);      // Fokker-Planck
        timings_at_current_timestep.push_back(0.);      // Output Routines
        timings_at_current_timestep.push_back(0.);      // Big Output Routines

        timing_indices.push_back(0.);
        timing_indices.push_back(1.);
        timing_indices.push_back(2.);
        timing_indices.push_back(3.);
    }
//--------------------------------------------------------------
Clock:: ~Clock(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    delete[] acceptabilitylist;
}
//--------------------------------------------------------------
Clock& Clock::advance(State1D& Y_current, Grid_Info& grid, 
    Output_Data::Output_Preprocessor &output, Export_Files::Restart_Facility &Re,
    Parallel_Environment_1D& PE) 
{
    // std::cout << "\n dt = " << _dt;

    end_of_loop_time_updates();

    timings_at_current_timestep[3] -= MPI_Wtime(); 
    if (current_time >= next_dist_out)
    {    
        if (!(PE.RANK())) cout << " \n Dist Output #" << t_out << "\n";
        output.distdump(Y_current, grid, t_out, current_time, _dt, PE);
        next_dist_out += dt_dist_out;
    }
    
    if (current_time >= next_big_dist_out)
    {
        if (!(PE.RANK())) cout << " \n Big Dist Output #" << t_out << "\n";
        output.bigdistdump(Y_current, grid, t_out, current_time, _dt, PE);
        next_big_dist_out += dt_big_dist_out;
    }
    
    if (current_time >= next_restart)
    {
        if (!(PE.RANK())) cout << " \n Restart Output #" << t_out << "\n";        
        Re.Write(PE.RANK(), t_out, Y_current, current_time);
        next_restart += dt_restart;
    }
    timings_at_current_timestep[3] += MPI_Wtime(); 


    
    if (Input::List().o_Exhist) Ex_history1D.push_back(Y_current.FLD(0).array());
    if (Input::List().o_Eyhist) Ey_history1D.push_back(Y_current.FLD(1).array());
    if (Input::List().o_Ezhist) Ez_history1D.push_back(Y_current.FLD(2).array());
    if (Input::List().o_Bxhist) Bx_history1D.push_back(Y_current.FLD(3).array());
    if (Input::List().o_Byhist) By_history1D.push_back(Y_current.FLD(4).array());
    if (Input::List().o_Bzhist) Bz_history1D.push_back(Y_current.FLD(5).array());
    if (Input::List().o_fhat0hist)  fhat0_history1D.push_back(output.px_radial_hat0(Y_current,grid));
    // if (Input::List().o_fhat1)  fhat0_history1D.push_back(Y_current.DF(0).getfhat1());

    timing_history.push_back(timings_at_current_timestep); std::fill(timings_at_current_timestep.begin(),timings_at_current_timestep.end(),0.);
    time_history.push_back(current_time);

    if (current_time >= next_out)
    {
        timings_at_current_timestep[2] -= MPI_Wtime(); 
        if (!(PE.RANK()))
        {
            cout << "\n dt = " << _dt;
            cout << " , Output #" << t_out;                        
        }
            
        if (Input::List().o_fhat0hist) 
        {
            output.histdump(fhat0_history1D, time_history, grid,  t_out, current_time, _dt, PE, "fhat0hist");
            fhat0_history1D.clear();
        }
        if (Input::List().o_Exhist) 
        {
            output.histdump(Ex_history1D, time_history, grid,  t_out, current_time, _dt, PE, "Exhist");
            Ex_history1D.clear();
        }
        if (Input::List().o_Eyhist) 
        {
            output.histdump(Ey_history1D, time_history, grid,  t_out, current_time, _dt, PE, "Eyhist");
            Ey_history1D.clear();
        }
        if (Input::List().o_Ezhist) 
        {
            output.histdump(Ez_history1D, time_history, grid,  t_out, current_time, _dt, PE, "Ezhist");
            Ez_history1D.clear();
        }
        if (Input::List().o_Bxhist) 
        {
            output.histdump(Bx_history1D, time_history, grid,  t_out, current_time, _dt, PE, "Bxhist");
            Bx_history1D.clear();
        }
        if (Input::List().o_Byhist) 
        {
            output.histdump(By_history1D, time_history, grid,  t_out, current_time, _dt, PE, "Byhist");
            By_history1D.clear();
        }
        if (Input::List().o_Bzhist) 
        {
            output.histdump(Bz_history1D, time_history, grid,  t_out, current_time, _dt, PE, "Bzhist");
            Bz_history1D.clear();
        }
        
        output(Y_current, grid, t_out, current_time, _dt, PE);

        timings_at_current_timestep[2] += MPI_Wtime(); 
        output.histdump(timing_history, time_history, timing_indices,  t_out, current_time, _dt, PE, "Timings");

        timing_history.clear();
        Y_current.checknan();

        next_out += dt_out;
        ++t_out;

        time_history.clear();

    }
    
    return *this;
}
//--------------------------------------------------------------
//--------------------------------------------------------------
Clock& Clock::advance(State2D& Y_current, Grid_Info& grid, 
    Output_Data::Output_Preprocessor &output, Export_Files::Restart_Facility &Re,
    Parallel_Environment_2D& PE) 
{
    end_of_loop_time_updates();

    timings_at_current_timestep[2] -= MPI_Wtime(); 
    if (current_time >= next_dist_out)
    {    
        if (!(PE.RANK())) cout << " \n Dist Output #" << t_out << "\n";
        output.distdump(Y_current, grid, t_out, current_time, _dt, PE);
        next_dist_out += dt_dist_out;
    }
    
    if (current_time >= next_big_dist_out)
    {
        if (!(PE.RANK())) cout << " \n Big Dist Output #" << t_out << "\n";
        output.bigdistdump(Y_current, grid, t_out, current_time, _dt, PE);
        next_big_dist_out += dt_big_dist_out;
    }
    
    if (current_time >= next_restart)
    {
        if (!(PE.RANK())) cout << " \n Restart Output #" << t_out << "\n";        
        Re.Write(PE.RANK(), t_out, Y_current, current_time);
        next_restart += dt_restart;
    }
    timings_at_current_timestep[2] += MPI_Wtime(); 



    if (Input::List().o_Exhist) Ex_history2D.push_back(Y_current.FLD(0).array());
    if (Input::List().o_Eyhist) Ey_history2D.push_back(Y_current.FLD(1).array());
    if (Input::List().o_Ezhist) Ez_history2D.push_back(Y_current.FLD(2).array());
    if (Input::List().o_Bxhist) Bx_history2D.push_back(Y_current.FLD(3).array());
    if (Input::List().o_Byhist) By_history2D.push_back(Y_current.FLD(4).array());
    if (Input::List().o_Bzhist) Bz_history2D.push_back(Y_current.FLD(5).array());

    timing_history.push_back(timings_at_current_timestep); std::fill(timings_at_current_timestep.begin(),timings_at_current_timestep.end(),0.);
    time_history.push_back(current_time);

    if (current_time >= next_out)
    {
        if (!(PE.RANK()))
        {
            cout << "\n dt = " << _dt;
            cout << " , Output #" << t_out;                        
        }
        
        if (Input::List().o_Exhist) 
        {
            output.histdump(Ex_history2D, time_history, grid,  t_out, current_time, _dt, PE, "Exhist");
            Ex_history2D.clear();
        }
        if (Input::List().o_Eyhist) 
        {
            output.histdump(Ey_history2D, time_history, grid,  t_out, current_time, _dt, PE, "Eyhist");
            Ey_history2D.clear();
        }
        if (Input::List().o_Ezhist) 
        {
            output.histdump(Ez_history2D, time_history, grid,  t_out, current_time, _dt, PE, "Ezhist");
            Ez_history2D.clear();
        }
        if (Input::List().o_Bxhist) 
        {
            output.histdump(Bx_history2D, time_history, grid,  t_out, current_time, _dt, PE, "Bxhist");
            Bx_history2D.clear();
        }
        if (Input::List().o_Byhist) 
        {
            output.histdump(By_history2D, time_history, grid,  t_out, current_time, _dt, PE, "Byhist");
            By_history2D.clear();
        }
        if (Input::List().o_Bzhist) 
        {
            output.histdump(Bz_history2D, time_history, grid,  t_out, current_time, _dt, PE, "Bzhist");
            Bz_history2D.clear();
        }
        
        output.histdump(timing_history, time_history, timing_indices,  t_out, current_time, _dt, PE, "Timings");
        timing_history.clear();

        output(Y_current, grid, t_out, current_time, _dt, PE);
        Y_current.checknan();

        next_out += dt_out;
        ++t_out;

        time_history.clear();

    }
    
    return *this;
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
    // dt_next = max(0.01,dt_next);
    failed_steps = 0;
    current_time += _dt;

    // std::cout << "\ndt_next = " << _dt;

    _dt = dt_next; 
    _success = 0;
}
//--------------------------------------------------------------
void Clock::do_step(State1D& Ystar, State1D& Y_new, State1D& Y_old,
                        VlasovFunctor1D_explicitE& vF, collisions_1D& cF, Parallel_Environment_1D& PE)
{
    if (Input::List().adaptive_dt)
    {    
        while (_success == 0)
        {
            Solver.take_step(Ystar, Y_new, current_time, _dt, vF, cF, PE);
            if (current_time > Input::List().adaptive_tmin) 
                update_dt(Y_old, Ystar, Y_new);
            else 
            {
                Y_old = Y_new;
                _success = 1;
            }                
        }
        if (Input::List().collisions)  
        {
            cF.advance(Y_new,current_time,_dt);
            PE.Neighbor_Communications(Y_new);
        }
    }
    else 
    {   
        
                                    timings_at_current_timestep[0] -= MPI_Wtime(); 
        Solver.take_step(Ystar, Y_new, current_time, _dt, vF, cF, PE);
                                    timings_at_current_timestep[0] += MPI_Wtime(); 
        
                                    timings_at_current_timestep[1] -= MPI_Wtime(); 
        if (Input::List().collisions)   
        {
            cF.advance(Y_new,current_time,_dt);
            PE.Neighbor_Communications(Y_new);
        }
                                    timings_at_current_timestep[1] += MPI_Wtime(); 

        //                             timings_at_current_timestep[0] -= MPI_Wtime(); 
        // Solver.take_step(Ystar, Y_new, current_time, 0.5*_dt, vF, cF, PE);
        //                             timings_at_current_timestep[0] += MPI_Wtime(); 
    }

}
//-------------------------------------------------------------------------------------------------------------------
// //-------------------------------------------------------------------------------------------------------------------
void Clock::do_step(State2D& Ystar, State2D& Y_new, State2D& Y_old,
                        VlasovFunctor2D_explicitE& vF, collisions_2D& cF, Parallel_Environment_2D& PE)
{
    if (Input::List().adaptive_dt)
    {    
        while (_success == 0)
        {
            Solver.take_step(Ystar, Y_new, current_time, _dt, vF, cF, PE);
            if (current_time > Input::List().adaptive_tmin) 
                update_dt(Y_old, Ystar, Y_new);
            else 
            {
                Y_old = Y_new;
                _success = 1;
            }                
        }
        if (Input::List().collisions)  
        {
            cF.advance(Y_new,current_time,_dt);
            PE.Neighbor_Communications(Y_new);
        }
    }
    else 
    {   
        
                                    timings_at_current_timestep[0] -= MPI_Wtime(); 
        Solver.take_step(Ystar, Y_new, current_time, _dt, vF, cF, PE);
        PE.Neighbor_Communications(Y_new);
                                    timings_at_current_timestep[0] += MPI_Wtime(); 
        
                                    timings_at_current_timestep[1] -= MPI_Wtime(); 
        if (Input::List().collisions)   
        {
            cF.advance(Y_new,current_time,_dt);
            PE.Neighbor_Communications(Y_new);
        }
                                    timings_at_current_timestep[1] += MPI_Wtime(); 

        //                             timings_at_current_timestep[0] -= MPI_Wtime(); 
        // Solver.take_step(Ystar, Y_new, current_time, 0.5*_dt, vF, cF, PE);
        //                             timings_at_current_timestep[0] += MPI_Wtime(); 
    }

}
// //-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
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

    std::cout << "\n acc = " << acceptability;
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
            dt_next = 0.9*_dt/pow(acceptability,0.2);   
        }
        else
        {
            dt_next = 0.9*_dt/pow(acceptability,0.25);
            
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
    // std::cout << "\n dt = " << _dt;

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

    maxval = max(fave.max(),fstarave.max());
    maxval = sqrt(maxval);

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
// //  Collect all of the terms
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
