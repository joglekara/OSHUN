/*!\brief  Particle Tracker - Declarations
 * \author PICKSC
 * \file   particletracker.h
 *
 * Includes declarations for spatial advection, electric field advection, bfield and current
 *
 */

#ifndef DECL_PARTICLETRACKER_H
#define DECL_PARTICLETRACKER_H

/** \addtogroup vfp1d
 *  @{
 */
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Spatial advection
class Particle_Pusher {
//--------------------------------------------------------------
public:
//      Constructors/Destructors
    Particle_Pusher(std::vector<double> xpos, 
                    std::vector<double> px,
                    std::vector<double> py,
                    std::vector<double> pz,
                    valarray<double> _xaxis,
                    Particle1D& Pin);
//          Advance
    void push(State1D& Y, double dt);
    
private:
    
    size_t numpar;

    vector<double> xaxis;

    double chargetomass;
};
//--------------------------------------------------------------
#endif
