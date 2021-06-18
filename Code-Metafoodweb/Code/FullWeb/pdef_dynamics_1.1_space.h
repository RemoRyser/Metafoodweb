/*
    Code for the article 'Landscape heterogeneity buffers biodiversity of meta-food-webs under global change through rescue and drainage effects'

    Copyright (C) 2016 Christian Guill & Florian D. Schneider
    Copyright (C) 2020 Remo Ryser
 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    last edit: 01.03.2020 by RR



 */

/*
// ************** Local foodweb code extended into space ***************************

Remo Ryser, Johanna Häussler, Markus Stark

*/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>      // for input/output from or to files
//#include <stdlib.h>
#include <iostream>		// for input/output on terminal
#include <sstream>
#include <vector>
#include <cstdlib>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>				// random number generator
#include <gsl/gsl_randist.h>			// random number distributions
#include <gsl/gsl_blas.h>				// linear algebra routines
#include <gsl/gsl_linalg.h>				// linear algebra
#include <gsl/gsl_sort_vector.h>		// vector operations
#include <gsl/gsl_odeiv.h>              // ODE solver

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */

double *web_calc(gsl_rng *r);
static void dynamics(double B[], double Bdot[], void *params, double t);
static void pdef_structure(gsl_rng *r,gsl_matrix *Ap, gsl_vector *mass, gsl_matrix *Up, gsl_vector *DmaxPlant);
static void output(gsl_matrix *SW, gsl_vector *mass, gsl_vector *con, double iniB[], double meanB[], double meanB_tot[], double CV[], double CV_tot[], double params[], int paramslength, double Biomass_Temp[], double RGG_params[], double net_growth_Temp[]);
static void getInputSpp();
static void getInputRicker();
static void getInputMass(gsl_vector *mass, gsl_matrix *Up,gsl_vector *DmaxPlant);
static void getInputLandscape();
static void getInputCoords(gsl_vector *coords_in_v);
static void getInputSW(gsl_matrix *SW_input);
static void getInputNutrient(gsl_vector *NS_input);

double get_random_parameter(double mean, double sigma, double low_cutoff, double high_cutoff ,gsl_rng *r);
static void mean_body_mass(double meanB[], gsl_vector *mass, double mean_body_masses[]);

static void set_parameters(gsl_rng *r, gsl_matrix *Ap, gsl_vector *mass, double params[], gsl_matrix *Up);
static void set_spatial_parameters(gsl_rng *r, gsl_matrix *SW, gsl_vector *Loc, gsl_vector *NS, double params[]);
static void initialise_biomass_densities(double B[], double iniB[], double params[], gsl_rng *r);
static void solve_ode(double B[], double meanB[],  double meanB_tot[], double CV[],  double CV_tot[], double params[], gsl_matrix *CovM, double Biomass_Temp[], double net_growth_Temp[]);
static void Extinction(double B[], int Num);
static void Extinct_Species(double B[], int Num, int Pat);
static int Cvode_rates(realtype t, N_Vector y, N_Vector ydot, void *params);

static void net_rates(gsl_vector *TBvec, gsl_vector *Ivec, gsl_vector *Dvec, gsl_vector *Xvec, gsl_vector *Emvec, gsl_vector *N_in, gsl_vector *N_out, gsl_vector *N_up, void *params, double t);

static void Prepare_timeseries_file(FILE *timeseries);
static void Write_timeseries_to_file(double B[], double t, FILE *timeseries);
static void Prepare_netrates_file(FILE *netrates);
static void Write_netrates_to_file(int p, gsl_vector *netvec, gsl_vector *Emvec_p, FILE *netrates);
static void Prepare_disprates_file(FILE *disprates);
static void Write_disprates_to_file(int i, int j, gsl_vector *Imvec_i, FILE *disprates);
static void Write_params_to_file(double params[], int paramslength);
static void Write_landscape_to_file(double params[]);
static void Write_SW_matrix_to_file(gsl_matrix *SW);


static void Generate_RGG_structure(gsl_rng *ran, gsl_vector *mass, gsl_vector *Loc, gsl_matrix *SW, gsl_matrix *SW_input, gsl_vector *con, double params[], double RGG_params[], gsl_vector *DmaxPlant);

static void calc_patch_distances(gsl_matrix *P_Dist, gsl_vector *Loc, double RGG_params[]);
static void calc_spatial_connectance(gsl_matrix *SW, gsl_vector *con, double RGG_params[]);


static void show_matrix(gsl_matrix *A, int Num, int N2);
static void show_vector(gsl_vector *A, int Num);
double calc_sd(gsl_vector *con, int START, int P, int div);
double calc_mean(gsl_vector *con, int START, int P, int div);

