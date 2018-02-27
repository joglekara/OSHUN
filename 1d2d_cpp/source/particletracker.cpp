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
void Particle_Pusher::push(State1D& Y, double dt)
{
	tk::spline splEx, splEy, splEz;
	tk::spline splBx, splBy, splBz;
	vector<double> exvec(xaxis.size()), eyvec(xaxis.size()), ezvec(xaxis.size());
	vector<double> bxvec(xaxis.size()), byvec(xaxis.size()), bzvec(xaxis.size());

	for (size_t ix(0); ix < xaxis.size(); ++ix)
	{
		exvec[ix] = Y.EMF().Ex()(ix).real();
		eyvec[ix] = Y.EMF().Ey()(ix).real();
		ezvec[ix] = Y.EMF().Ez()(ix).real();
		bxvec[ix] = Y.EMF().Bx()(ix).real();
		byvec[ix] = Y.EMF().By()(ix).real();
		bzvec[ix] = Y.EMF().Bz()(ix).real();
	}

	/// Create Splines
	splEx.set_points(xaxis,exvec); splEy.set_points(xaxis,eyvec); splEz.set_points(xaxis,ezvec);
	splBx.set_points(xaxis,bxvec); splBy.set_points(xaxis,byvec); splBz.set_points(xaxis,bzvec);

	#pragma omp parallel for num_threads(Input::List().ompthreads)
	for (int ip = 0; ip < numpar; ++ip)
	{
		double tempEx(0.), tempEy(0.), tempEz(0.);
		double tempBx(0.), tempBy(0.), tempBz(0.);
		// std::cout << "min = " << xmin << ", position = " << Y.particles().x(ip) << ", Pin.ishere = " << Y.particles().ishere(ip) << "\n";
		if (Y.particles().ishere(ip))
		{	

			Y.particles().x(ip) = Y.particles().x(ip) + 0.5*dt*Y.particles().px(ip);
			
			tempEx = splEx(Y.particles().x(ip)); tempEy = splEy(Y.particles().x(ip)); tempEz = splEz(Y.particles().x(ip));
			tempBx = splBx(Y.particles().x(ip)); tempBy = splBy(Y.particles().x(ip)); tempBz = splBz(Y.particles().x(ip));

			double Lorentz_x = tempEx + Y.particles().py(ip) * tempBz - Y.particles().pz(ip) * tempBy;
			double Lorentz_y = tempEy + Y.particles().pz(ip) * tempBx - Y.particles().px(ip) * tempBz;
			double Lorentz_z = tempEz + Y.particles().px(ip) * tempBy - Y.particles().py(ip) * tempBx;

			Y.particles().px(ip) = Y.particles().px(ip) + dt*chargetomass*Lorentz_x;
			Y.particles().py(ip) = Y.particles().py(ip) + dt*chargetomass*Lorentz_y;
			Y.particles().pz(ip) = Y.particles().pz(ip) + dt*chargetomass*Lorentz_z;			

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
