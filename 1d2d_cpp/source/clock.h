/*!\brief  Functors for various time-integration methodds - Declarations
 * \author PICKSC
 * \file   clock.h
 *
 * Includes functors for fully explicit, implicit B, implicit E
 *
 */
#ifndef OSHUN_CLOCK_H
#define OSHUN_CLOCK_H

//**************************************************************
//  This Clock controls an iteration loop, given an initial
//  time tout_start*dt_out it evaluates the number of time-
//  steps it takes to get to (tout_start+1)*dt_out and can
//  be incremented;

class Clock {
public:
//      Constructor
    Clock(double starttime, double __dt, double abs_tol, double rel_tol, size_t _maxfails, State1D& Y);
    ~Clock();

    void end_of_loop_time_updates();
    // void end_of_loop_output(Output_Data::Output_Preprocessor output);
    // void end_of_loop_dist_output(Output_Data::Output_Preprocessor output);
    // void end_of_loop_bigdist_output(Output_Data::Output_Preprocessor output);

    void do_step(State1D& Ystar, State1D& Y_new, State1D& Y_old,
                        VlasovFunctor1D_explicitE& vF, collisions_1D& cF, Parallel_Environment_1D& PE);

    double check_temperature(const State1D& Ystar, const State1D& Y, double& maxval);
    double check_flds(const State1D& Ystar, const State1D& Y, double& maxval);
    double check_flds(const State2D& Ystar, const State2D& Y, double& maxval);
    double check_last_harmonic(const State1D& Ystar, const State1D& Y, double& maxval);
    double check_js(const State1D& Ystar, const State1D& Y);

    void update_dt(State1D& Y_old, const State1D& Ystar, State1D& Y_new);
    void update_dt(const State2D& Y_old, const State2D& Ystar, State2D& Y_new);

    Clock& operator++();
    Clock& advance(State1D& Y_current, Grid_Info& grid, 
    Output_Data::Output_Preprocessor &output, Export_Files::Restart_Facility &Re,
    Parallel_Environment_1D& PE);

    double dt() {return _dt;}
    double nextdt() {return dt_next;}
    double time() {return current_time;}
    int success() {return _success;}

private:

    double current_time, dt_next, _dt;
    double atol, rtol, acceptability, err_val;
    
    size_t failed_steps, max_failures;

    int _success;

    size_t Nbc;
    int world_rank, world_size;

    // ARK32 Solver;
    // ARK43 Solver;

    // RKCK45 Solver;
    // RKDP85 Solver;
    RK4C Solver;

    int tout_start;
    size_t t_out;
    double dt_out, dt_dist_out, dt_big_dist_out, dt_restart;
    double next_out, next_dist_out, next_big_dist_out;
    double next_restart;
    double start_time;

    vector<valarray<complex<double> > > Ex_history, Ey_history, Ez_history;
    vector<valarray<complex<double> > > Bx_history, By_history, Bz_history;
    vector<vector<double> > timing_history;

    vector<double>                      time_history;
    vector<double>                      timings_at_current_timestep;
    vector<double>                      timing_indices;

    double* acceptabilitylist;

};
//--------------------------------------------------------------

//
//**************************************************************
#endif //OSHUN_CLOCK_H

