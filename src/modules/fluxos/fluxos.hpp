// Copyright 1992: Cornel Beffa and Roland Faeh
// Copyright 2013: Kashif Shaad and Diogo Costa
// Copyright 2019, Diogo Costa

// This program, FLUXOS, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"

#include<iostream>
#include<math.h>
#include<armadillo>
#include<string>
#include<memory>


/**
 * Fluxos
 */
class fluxos : public module_base
{
REGISTER_MODULE_HPP(fluxos)
public:
    fluxos(config_file cfg);

    ~fluxos();

    void run(mesh& domain);
    void init(mesh& domain);


    std::unique_ptr<arma::Mat<double>> forcing;
    std::unique_ptr<arma::Mat<double>> mapping;
    class data : public face_info {
    public:
       arma::uvec rastercells;
    };
    struct declavar
    {
        int nx,ny;
        int mx,my,dxy,arbase,
                ntim;                                       // maximum time step (seconds)
        float gacc = 9.80665,                           // gravitational acceleration
                cfl,                                        // Courant condition
                cvdef, nuem;                                // for turbulence calc; nuem is molecular viscosity
        std::unique_ptr<arma::Mat<double>> z,zb,h,
                u,v,                                        // velocities
                p,q,                                        // discharge at cell center: u*h [m3/s/m]
                pj,qj,                                      // discharges at cell face: u*h [m3/s/m]
                us,                                         // shear stress velocity
                dh,dp ,dq,                                  // changes in h[ix][iy], p[ix][iy] and q[ix][iy]
        //sbmx,sbmy,                                  // for calc of weight of water (bed slope term) (solver_wet)
                ks, //cfri                                  // Friction (Chezy model is not being used for now)
                fe_1,fe_2,fe_3,fn_1,fn_2,fn_3;
        std::unique_ptr<arma::Mat<float>> ldry;
        double hdry,                                    //minimum water depth
                dtfl,tim;                                   // timestep for flow computation (s)
    } ds;

    size_t print_step;



    // helper fluxos functions


    void initiation(declavar& ds);
    void solver_dry(declavar& ds, unsigned int ix, unsigned int iy);
    void solver_wet(declavar& ds, unsigned int ix, unsigned int iy);
    void flow_solver(declavar& ds);
    void write_results(declavar& ds, int print_tag);
};
