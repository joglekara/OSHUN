/*!\brief  Particle Pusher Routines - Definitions
* \author  PICKSC
 * \date   March, 2017
 * \file   particletracker.cpp
 *
 * In here are the routines for the particle tracker object.
 *
 */

//  Standard libraries
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <algorithm>
#include <cstdlib>
#include <cfloat>

#include <math.h>
#include <map>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"
#include "external/spline.h"

// Declarations
#include "input.h"
#include "state.h"
#include "particletracker.h"


////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////
//
Particle_Pusher::Particle_Pusher(std::vector<double> xpos, 
	std::vector<double> px, std::vector<double> py, std::vector<double> pz,
	valarray<double> _xaxis,
	Particle1D& Pin): numpar(py.size()),
	chargetomass(Pin.charge()/Pin.mass())
{
	
	for (int ix(0); ix < _xaxis.size(); ++ix)	
	{
		xaxis.push_back(_xaxis[ix]);
	}


	for (int ip(0); ip < numpar; ++ip){
		
		Pin.x(ip) = xpos[ip];
		if (Pin.x(ip) < xaxis[xaxis.size()-1] && Pin.x(ip) > xaxis[0])
		{
			// std::cout << "\nParticle " << ip << " is in cell with xmin = " << xmin;
			Pin.ishere(ip) = 1;
		}
		else Pin.ishere(ip) = 0;

		Pin.px(ip) = px[ip];
		Pin.py(ip) = py[ip];
		Pin.pz(ip) = pz[ip];


		// std::cout << "min = " << xmin << ", and Pin.ishere(ip)" << Pin.ishere(ip) << "\n";
	}
}
/**
 * Pusher
 */
void Particle_Pusher::push(State1D& Y, double dt){


	tk::spline splEx;
	vector<double> exvec(xaxis.size());

	for (size_t ix(0); ix < xaxis.size(); ++ix)
	{
		exvec[ix] = Y.EMF().Ex()(ix).real();
	}

	splEx.set_points(xaxis,exvec);

	#pragma omp parallel for num_threads(Input::List().ompthreads)
	for (int ip = 0; ip < numpar; ++ip)
	{
		double tempEx(0.);
		// std::cout << "min = " << xmin << ", position = " << Y.particles().x(ip) << ", Pin.ishere = " << Y.particles().ishere(ip) << "\n";
		if (Y.particles().ishere(ip))
		{	

			Y.particles().x(ip) = Y.particles().x(ip) + 0.5*dt*Y.particles().px(ip);
			
			tempEx = splEx(Y.particles().x(ip));
			Y.particles().px(ip) = Y.particles().px(ip) + dt*chargetomass*tempEx;

			Y.particles().x(ip) = Y.particles().x(ip) + 0.5*dt*Y.particles().px(ip);
			
			if (Y.particles().x(ip) >= xaxis[xaxis.size()-1] )
			{	
				// std::cout << "\n Leaving right! \n";
				Y.particles().ishere(ip) = 0;
				Y.particles().goingright(ip) = 1;
			}
			else if (Y.particles().x(ip) < xaxis[0])
			{
				// std::cout << "\n Leaving left! \n";
				Y.particles().ishere(ip) = 0;
				Y.particles().goingright(ip) = -1;
			}
			// if (ip == 2) 
			// std::cout << "xnew = " << Y.particles().x(ip) << "\n";
		}
	}
}
//
