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

    last edit: 01.07.2021 by RR



 */

/*
// ************** Local foodweb code extended into space ***************************

Remo Ryser, Johanna Häussler, Markus Stark

*/



#include "pdef_dynamics_1.1_space.h"

//  basic simulation variables:

int S;                        // total number of species
int S_c;                      // number consumers (animals, predators)
int S_b;                      // number basal species (plants)
int N=1;                      // number of nutrients

double tend=100000;//5e5;           // time at which simulation runs stop
double teval=10000;//2e4;           // length of time interval from which to determine mean(B) and CV(B)
double Delta_t=0.1;//0.01;         // step size in ODE integration (internal step size can be less or larger)

const double eps_rel= 1e-10;  // relative error tolerance for ODE solver
const double eps_abs = 1e-10; // absolute error tolerance for ODE solver
const double EXTINCT = 1e-20;  // extinction threshold for biomass densities
const double GLOBAL_EXTINCT = EXTINCT;//1e-5;  // global extinction threshold for species biomass densities
//const double Max_Disp = 1e8;
int EXT_FLAG;


//  parameters that determine the network topology
double zeta_c=6;              // log_10 of the mean of the consumer (animal) body masses
double sigma_c=3;             // width of the distribution of consumer (animal) body masses
double cutoff_c=1e5;          // half relative cutoff for distribution of consumer body masses
double zeta_b=5;              // log_10 of the mean of the basal (plant) body masses
double sigma_b=3;             // width of the distribution of basal (plant) body masses
double cutoff_b=1e5;          // half relative cutoff for distribution of basal body masses

double m_p_min = 0;           // minimla and maximal log10 body masses in case of uniform distributions
double m_a_min = 2;
double m_p_max = 6;
double m_a_max = 12;

double cutoff=0.02;           // cutoff of the Ricker curve for setting a link between predator and prey
double R_opt=100;             // optimal predator-prey body-mass ratio
double g;                   // width of the Ricker curve (higher g-value -> narrower curve) read in from txt file

double f_herbiv = 0.0;//0.50;       // fraction of species that are strict herbivores
double f_pred = 0.0;          // fraction of species that are strict predators

//  parameters of the functional response
double a_0=15;                // scaling factor for the attack rate
double a_c= 0.42;                   // exponent for predator body-mass scaling of the attack rate (value is drawn from distribution)
//double mean_a_c = 0.19;       // (mean for the above distribution)
//double sigma_a_c = 0.04;      // (standard deviation for the above distribution)
double a_p= 0.19;                   // exponent for prey body-mass scaling of the attack rate (value is drawn from distribution)
//double mean_a_p = 0.19;       // (mean for the above distribution)
//double sigma_a_p = 0.03;      // (standard deviation for the above distribution)
double a_plant = 3500;          // constant attack rate coefficient for plants
//double f;                       // environmental factor to move the intercept of the attack rate

double h_0=0.4;               // scaling factor for the handling time
double h_c= -0.48;                   // exponent for predator body-mass scaling of the handling time (value is drawn from distribution)
//double mean_h_c = -0.48;      // (mean for the above distribution)
//double sigma_h_c = 0.03;      // (standard deviation for the above distribution)
double h_p= -0.66;                   // exponent for prey body-mass scaling of the handling time (value is drawn from distribution)
//double mean_h_p = -0.66;      // (mean for the above distribution)
//double sigma_h_p = 0.02;      // (standard deviation for the above distribution)

double hill=1.1;                  // Hill coefficient (value is drawn from distribution)
//double mu_hill=1.5;           // (mean for the above distribution)
//double sigma_hill=0.2;        // (standard deviation for the above distribution)
double C_interference=0;        // interference competition (value is drawn from distribution)
//double mu_Cint=0.8;           // (mean for the above distribution)
//double sigma_Cint=0.2;        // (standard deviation for the above distribution)

double x_resp_a = 0.141;      // intercept of animal respiration rates
double x_resp_p = 0.138;      // intercept of producer respiration rates
double e_p=0.545;              // assimilation efficiency for plant resources
double e_a=0.906;              // assimilation efficiency for animal resources

double B_b=10;                // intercept initial basal biomass densities
double B_c=10;                // intercept initial consumer biomass densities

//  parameters of the nutrient model
double C1=1;                  // content of first nutrient in plants
double C2=0.5;                // content of second nutrient in plants
//double D_N;              // nutrient turnover rate
double K_min=0.1;             // minimal nutrient uptake half saturation density
double K_rel=0.1;             // width of interval for nutrient uptake half saturation densities
//double mu_S=10;               // mean nutrient supply concentration
//double sigma_S=2;             // standard deviation of nutrient supply concentration

//  default parameters for loops and flags
int GF=0;                     // global flag, for testing purposes
double cv_flag;               // indicates potentially pathological coefficients of variation
int TIMESERIES = 0;           // if true, for the last network a file with the time series will be created
int WRITENETRATESTOFILE = 0;   // if true, for the last network a file with the net growth rate will be created
int WRITEPARAMSTOFILE = 0;    // if true, the parameter array is written to a file
int WRITELANDSCAPETOFILE = 1; // 1: generated RGG is written to a file
int VAR_COEFF = 1;            // if true, the coefficient of variation will be calculated; if false, the simulations run a bit faster
int UNIFORM_MASSES = 1;       // 0: body masses drawn from log-normal distribution; 1: log-body masses drawn from uniform distribution
int CANNIBALISM = 0;          // 0: all cannibalism links are explicitly removed // do we need this here ???
int VARIABILITY = 0;		  // 0: no calculation of alpha-, beta- and gamma-variability 1: Calculation of variability introduced by Wang&Loreau 2014
int INPUT_SPP = 1;                // 0: body masses randomly generated; 1: body masses imported
int INPUT_LANDSCAPE = 1;          // 0: landscape randomly generated; 1: landscape imported
int INPUT_SW =0;                    // 0: Success weighted dispersal generated; 1: Sucess weighted dispersal imported


// basic simulation variables for spatial dynamics:
int Z; 						  // number of patches (=number of subpopulations)
int D;	 					  // dimension of the biotic system (D=S*Z)
int DN; 					  // dimension of the abiotic system (D=N*Z)


const double exp_suc = 0.33;                // exponent for scaling dispersal success half saturation density ?? extinction_boundary


// parameters for plant migration dynamics
double emigr_a_S_b=0;                         // max. emigration rate
double emigr_b_S_b=0;                         // shape of curve for emigration
//double mean_emigr_b_S_b=0;                 // (mean for the above distribution)
//double sigma_emigr_b_S_b=10;                 // (standard deviation for the above distribution)
//double emigr_a_S_c;                         // max. emigration rate
double emigr_b_S_c=10;                         // shape of curve for emigration
//double mean_emigr_b_S_c=15;                 // (mean for the above distribution)
//double sigma_emigr_b_S_c=5;                 // (standard deviation for the above distribution)
// parameters for animal migration dynamics


const double MIN_MIGRATION = 1e-10;         // minimum migration threshold

// disperal parameters for RGG
double D_0;                     // minimal dispersal distance -- defined in Generate_RGG_structure
double eps = 0.05;              // scaling factor of body mass for maximal dispersal distance d_max_i=D_0*m_i^eps
double theta = 1;               // shape of curve for dispersal success

const double mass_crit = 0;     // body mass threshold at which scaling for extinction boundary changes
char* RGG_Id = getenv("LANDSCAPE");
char* Web_Id = getenv("WEB");
char* landscape_s = getenv("LANDSCAPESIZE"); 		// read in landscape size from bash
char* directory = getenv("OUTPUTDIR"); 		// diretory for output files
char* inputdir = getenv("INPUTDIR");        // diretory for input files
char* name= getenv("NAME"); 			// name completion for output files
char* bash_seed= getenv("SEED");            //read in seed from bash
//char* directory = "/home/leco35/Uni/GIT/Data/";
int seed = atoi(bash_seed);              // different replicates need different random numbers!
char* bash_D_N = getenv("NUTRATE");
double NutRate = atof(bash_D_N);
double D_N = NutRate;
char* bash_Supply = getenv("SUPPLY");
double NutSupply=atof(bash_Supply);
char* f_bash= getenv("F_BASH");
double f = atof(f_bash);            // environmental factor to move the intercept of the attack rate
char* emigr_bash=getenv("EMIGR");
double emigr_a_S_c=atof(emigr_bash);
char* bash_loss=getenv("LOSS");
int EMIGRATIONLOSS = atoi(bash_loss);
char* Disp=getenv("DISP");
char* NUT=getenv("NUT");
char* bash_kill=getenv("KILL");
int KILLSWITCH = atoi(bash_kill);           // 0: initialising all species on all patches; 1: kill some species on patches during initialisation

//int constseed = 125;
//int RGG_Id= atoi(RGG_Id);
//int Web_Id= atoi(Web_Id);
//char directory[256] = "/Users/Remo/Patchdistance/Simulation/";  // directory for the output files
//char filename[256]="pdef_";          // name of output file
//char name[256] = "run03";               // name completion for output files
//int  seed = 1450;
//char inputdir[256]      = "/Users/Remo/Patchdistance/Simulation"; // directory for the input files
//FILE *netrates;
//FILE *disprates;



//  ********** Main function **********
int main(int argc, char* argv[])
{
    int i;


    gsl_rng_default_seed = seed;                  // seed the rng
    gsl_rng *r=gsl_rng_alloc(gsl_rng_default);

    if(INPUT_SPP)
        getInputSpp();
    else
    {

        S_c = 6;                       // number of consumer (animal) species
        S_b = 2;                       // number of basal (plant) species // set basal species to S-S_c
        S = S_c + S_b;                  // total number of species
    }
    
    if(INPUT_LANDSCAPE)
        getInputLandscape();
    
    else
    {
    Z = 50 ; 	    // number of patches
    }
    

    D = S*Z;								// dimension of the biotic system
    DN = N*Z;								// dimension of the abiotic system

    web_calc(r);					// initialise model

    gsl_rng_free(r);

    return(0);
}

//  ********** generate a network, simulate population dynamics, and evaluate the final network **********
double *web_calc(gsl_rng *r)
{

    //gsl_rng_default_seed = constseed;                  // seed the rng with a constant
    //gsl_rng *f=gsl_rng_alloc(gsl_rng_default);

    gsl_matrix *Ap=gsl_matrix_calloc(S,S);                     // adjacency matrix
    gsl_vector *mass=gsl_vector_calloc(S);                     // mean body masses of the species
    gsl_vector *con=gsl_vector_calloc(S);                      // spatial connectance of each species
    gsl_matrix *Up=gsl_matrix_calloc(S_b,N);                    // Nutrient uptake half saturation density
    gsl_vector *DmaxPlant=gsl_vector_calloc(S_b);               //Maximal dispersal distance for Plants
    
    gsl_matrix *SW        = gsl_matrix_calloc(S*Z,Z);		     	// success-weighted dispersal matrix
    gsl_matrix *SW_input        = gsl_matrix_calloc(S*Z,Z);                 // success-weighted dispersal matrix
    gsl_matrix *CovM 		= gsl_matrix_calloc((S+1)*Z, Z);		// Covariance-Matrix, Submarix for each species + one for every species + one for all
    gsl_vector *Loc       = gsl_vector_calloc(2*Z);		     		// patch location (XY)
    gsl_vector *NS        =gsl_vector_calloc(DN);                   // patch-specific nutrient supply
   

    double *B=(double *) calloc(D+DN,sizeof(double));          		// biomass densities
    double *iniB=(double *) calloc(D+DN,sizeof(double));       		// initial biomass denisties
    double *meanB=(double *) calloc(D+DN,sizeof(double));      		// mean biomass densities
    double *meanB_tot=(double *) calloc(S+N,sizeof(double));   		// mean total biomass densities
    double *VarCoeffB=(double *) calloc(D+DN,sizeof(double));  		// coefficients of variation
    double *VarCoeffB_tot=(double *) calloc(S+N,sizeof(double)); 	// coefficients of variation
    double *net_growth_Temp=(double *) calloc(3*D,sizeof(double));		// values for the net growth rate
    double *Biomass_Temp= (double *) calloc(3*(D+DN),sizeof(double));   // biomass densities of three time steps which are separated 10.000 time steps
    
    double *RGG_params=(double *) calloc(10,sizeof(double));          // RGG_parameters (spatial connectance, mean distance, etc.)

    int paramslength = 2*S*S+7*S+S_b*N+N+S_b+2+S*Z*Z+2*Z+7+D+DN+DN;
    double *params=(double *) calloc(paramslength,sizeof(double));  // parameters for population dynamics

    
    printf("1 Start simulation: S_b = %d; S_c = %d; N = %d; Z = %d\n",S_b,S_c,N,Z);

    pdef_structure(r,Ap,mass,Up,DmaxPlant);                                   // generate random network structure
    printf("2 Foodweb generated\n");

    set_parameters(r,Ap,mass,params,Up);                           // set random parameters for the simulation run
    printf("3 Parameters set\n");


    Generate_RGG_structure(r,mass,Loc,SW, SW_input,con,params,RGG_params,DmaxPlant);
    printf("4 RGG generated\n");
    

    set_spatial_parameters(r,SW,Loc,NS, params);                // set random spatial parameters for the simulation run
    printf("5 Spatial parameters set\n");
    

    initialise_biomass_densities(B, iniB, params, r);            // set random initial values of biomass densities
    printf("6 Biomass initialised\n");

    solve_ode(B,meanB,meanB_tot,VarCoeffB,VarCoeffB_tot,params,CovM,Biomass_Temp,net_growth_Temp);   // apply solver for integration of the ODE

    printf("7 ODE solved\n");

    output(SW,mass,con,iniB,meanB,meanB_tot,VarCoeffB,VarCoeffB_tot,params,paramslength,Biomass_Temp,RGG_params,net_growth_Temp);
    printf("8 Output generated\n");

    double Surv=0;
    double con_count=0;
    double bas_count=0;
    for(int i=0; i<S; i++)
    {
        double B_Species = 0;

        for (int m=0; m<Z; m++)
        {
            B_Species = B_Species + B[i+S*m];
        }
        if(B_Species > GLOBAL_EXTINCT)
            Surv++;
        if(i >= S_b && B_Species > GLOBAL_EXTINCT)
            con_count++;
        if(i < S_b && B_Species > GLOBAL_EXTINCT)
            bas_count++;
    }

    printf("global persistence = %.7g, number of basal species = %.7g, number of consumer species = %.7g\n",Surv/S,bas_count,con_count);

//  *********** Free memories ***********
    gsl_matrix_free(Ap);
    gsl_matrix_free(Up);
    gsl_vector_free(DmaxPlant);
    gsl_matrix_free(SW);
    gsl_matrix_free(SW_input);
    gsl_matrix_free(CovM);
    gsl_vector_free(mass);
    gsl_vector_free(con);
    gsl_vector_free(Loc);
    gsl_vector_free(NS);
    free(params);
    free(RGG_params);
    free(B);
    free(iniB);
    free(meanB);
    free(meanB_tot);
    free(net_growth_Temp);
    free(VarCoeffB);
    free(VarCoeffB_tot);
    free(Biomass_Temp);

    return 0;
}

//  ********** generate the body-mass based network structure **********
static void pdef_structure(gsl_rng *r, gsl_matrix *Ap, gsl_vector *mass, gsl_matrix *Up, gsl_vector *DmaxPlant)
{

    int i,j,flag=0;
    double Wkeit,sigma_i,zeta_act;
    double temp1,temp2,m_opt,a_max,R;
    double m_crit;

    while(flag==0)
    {
        flag=1;

        gsl_matrix_set_zero(Ap);

        if(INPUT_SPP)
        {
            getInputMass(mass, Up, DmaxPlant);
            getInputRicker();
        }
        else
        {
            gsl_vector_set_all(mass,0);

//  ********* determine body masses *********
            if(UNIFORM_MASSES == 0)
            {
                zeta_act=zeta_c*log(10);                                          // calculate log(mean) from log10(mean)
                for(i=0; i<S_c; i++)
                {
                    temp1=gsl_ran_lognormal(r,zeta_act,sigma_c);
                    if(temp1>exp(zeta_act)/cutoff_c&&temp1<exp(zeta_act)*cutoff_c)  // check whether mass is within certain boundaries
                        gsl_vector_set(mass,S_b+i,temp1);
                    else
                        i--;
                }
                gsl_sort_vector(mass);

                gsl_vector_view mass_b_vec=gsl_vector_subvector(mass,0,S_b);
                gsl_vector *mass_b=&mass_b_vec.vector;

                zeta_act=zeta_b*log(10);                                          // calculate log(mean) from log10(mean)
                for(i=0; i<S_b; i++)
                {
                    temp1=gsl_ran_lognormal(r,zeta_act,sigma_b);
                    if(temp1>exp(zeta_act)/cutoff_b&&temp1<exp(zeta_act)*cutoff_b)  // check whether mass is within certain boundaries
                        gsl_vector_set(mass_b,i,temp1);
                    else
                        i--;
                }
                gsl_sort_vector(mass_b);
            }
            else
            {
                for(i = 0; i< S_c; i++)
                    gsl_vector_set(mass,S_b+i,pow(10.,gsl_ran_flat(r,m_a_min,m_a_max)));
                gsl_sort_vector(mass);

                gsl_vector_view mass_b_vec=gsl_vector_subvector(mass,0,S_b);
                gsl_vector *mass_b=&mass_b_vec.vector;

                for(i = 0; i< S_b; i++)
                    gsl_vector_set(mass_b,i,pow(10.,gsl_ran_flat(r,m_p_min,m_p_max)));
                gsl_sort_vector(mass_b);
            }
        }

//  ********** fill the adjacency matrix with Ricker attack rates *************
        for(i=0; i<S_c; i++)
        {
            temp1=gsl_vector_get(mass,S_b+i);
            a_max=pow(gsl_vector_get(mass,S_b+i)/R_opt,0.25);

            for(j=0; j<S; j++)
            {
                temp2=gsl_vector_get(mass,j);
                R = temp1/temp2;
                if(pow((R / R_opt) * exp(1 - (R / R_opt)), g) >= cutoff)
                    gsl_matrix_set(Ap,S_b+i,j,1);
            }

            if(gsl_rng_uniform(r) < (f_herbiv + f_pred))                      // combined probability to mess with this species' links
            {
                if(gsl_rng_uniform(r) < f_herbiv/(f_herbiv+f_pred))             // either make it a strict herbivore...
                {
                    if(UNIFORM_MASSES == 0)
                        m_crit = pow(10,zeta_b) * cutoff_b * R_opt;
                    else
                        m_crit = pow(10,m_p_max) * R_opt;

                    if(gsl_vector_get(mass,S_b+i) < m_crit)
                    {
                        for(j = S_b; j < S; j++)
                            gsl_matrix_set(Ap,S_b+i,j,0);                             // and remove all links from non-plant resources
                    }
                }
                else                                                            // ... or make it a strict carnivore
                {
                    for(j = 0; j < S_b; j++)
                        gsl_matrix_set(Ap,S_b+i,j,0);                               // and remove all links from plant resources
                }
            }

            if(!CANNIBALISM)
                gsl_matrix_set(Ap,S_b+i,S_b+i,0);                               // remove canibalistic links if necessary

            gsl_vector_view tempp=gsl_matrix_row(Ap,S_b+i);                   // reject networks with consumers or predators without prey
            flag=flag*(1-gsl_vector_isnull(&tempp.vector));
        }

        for(i=0; i<S_b; i++)
        {
            gsl_vector_view tempp = gsl_matrix_column(Ap,i);                  // reject networks with uncontrolled basal species
            flag = flag*(1-gsl_vector_isnull(&tempp.vector));
        }
    }

    return;
}

//  ********** write all parameters required for the dynamics to the array 'params' **********
static void set_parameters(gsl_rng *r, gsl_matrix *Ap, gsl_vector *mass, double params[], gsl_matrix *Up)
{
    int i,j,k;
    double temp1,temp2,temp3,temp4,temp5,R;

    gsl_matrix *A=gsl_matrix_calloc(S,S);                                 // attack rates
    gsl_matrix *H=gsl_matrix_calloc(S,S);                                 // handling times

    gsl_matrix_memcpy(A,Ap);

    //hill = get_random_parameter(mu_hill,sigma_hill,1,2,r);
    //C_interference = get_random_parameter(mu_Cint,sigma_Cint,mu_Cint - 3*sigma_Cint, mu_Cint + 3*sigma_Cint,r);

    //a_c = get_random_parameter(mean_a_c, sigma_a_c, mean_a_c - 3*sigma_a_c, mean_a_c + 3*sigma_a_c,r);
    //a_p = get_random_parameter(mean_a_p, sigma_a_p, mean_a_p - 3*sigma_a_p, mean_a_p + 3*sigma_a_p,r);

    //h_c = get_random_parameter(mean_h_c, sigma_h_c, mean_h_c - 3*sigma_h_c, mean_h_c + 3*sigma_h_c,r);
    //h_p = get_random_parameter(mean_h_p, sigma_h_p, mean_h_p - 3*sigma_h_p, mean_h_p + 3*sigma_h_p,r);

    for(i=0; i<S; i++)
    {
        if(i >= S_b)                                                        // the following lines are only for non-basal species
        {
            gsl_vector_view tempp=gsl_matrix_row(A,i);
            temp1=gsl_blas_dasum(&tempp.vector);                              // temp1 stores the number of prey species of predator i
            gsl_vector_scale(&tempp.vector,1/temp1);                          // reduce attack rates for generalists

            

            for(j=0; j<S; j++)
            {
                if(j < S_b)
                {
                    temp1=pow(gsl_vector_get(mass,i),a_p);
                    temp2 = a_plant;
                }
                else
                {
                    temp1 = pow(gsl_vector_get(mass,i),a_c);
                    temp2 = a_0 * pow(gsl_vector_get(mass,j),a_c);
                }
                
                R = gsl_vector_get(mass,i) / gsl_vector_get(mass,j);            // predator-prey body-mass ratio for Ricker curve
                temp3 = f * temp2 * temp1 * pow((R / R_opt) * exp(1 - (R / R_opt)),g);
                gsl_matrix_set(A,i,j,temp3*gsl_matrix_get(A,i,j));          // attack rates

                temp4 = h_0 * pow(gsl_vector_get(mass,i),h_c) * pow(gsl_vector_get(mass,j),h_p);
                gsl_matrix_set(H,i,j,temp4);                                    // handling times
            }
        }

        for(j=0; j<S; j++)
        {
            *(params+i*S+j)=gsl_matrix_get(A,i,j);
            *(params+S*S+i*S+j)=gsl_matrix_get(H,i,j);
        }

        if(i < S_b)
        {
            *(params+2*S*S+i)=x_resp_p*pow(gsl_vector_get(mass,i),-0.25);     // plant respiration rates
            *(params+2*S*S+2*S+i)=1;                                          // identifyer for basal species
            *(params+2*S*S+5*S+i) = e_p;                                      // assimilation efficiency for plant resources
        }
        else
        {
            *(params+2*S*S+i)=x_resp_a*pow(gsl_vector_get(mass,i),-0.305);     // animal respiration rates
            *(params+2*S*S+5*S+i) = e_a;                                      // assimilation efficiency for animal resources
        }

        *(params+2*S*S+S+i)=C_interference;                                 // interference competition
        *(params+2*S*S+3*S+i)=hill;                                         // Hill coefficient
        *(params+2*S*S+4*S+i)=gsl_vector_get(mass,i);                       // body masses
    }

    for(i=0; i<S_b; i++)                                                  // the following lines are only for basal species
    {
        if(INPUT_SPP)
        {
        for(j=0; j<N; j++)
            *(params+2*S*S+6*S+i*N+j)=gsl_matrix_get(Up,i,j);         // nutrient uptake half saturation densities
        }
        else
        {
        for(j=0; j<N; j++)
        *(params+2*S*S+6*S+i*N+j)=K_min+K_rel*gsl_rng_uniform(r);         // nutrient uptake half saturation densities
        }
    }

    for(i=0; i<N; i++)
    {
        temp1=NutSupply;
        
        *(params+2*S*S+6*S+S_b*N+i)=temp1;                                // nutrient supply concentrations
    }


    for(i=0; i<S_b; i++)
        *(params+2*S*S+6*S+S_b*N+N+i)=pow(gsl_vector_get(mass,i),-0.25);    // max. nutrient uptake rates

    *(params+2*S*S+6*S+S_b*N+N+S_b) = C1;
    *(params+2*S*S+6*S+S_b*N+N+S_b+1) = C2;


    for(i=0; i<S; i++) 						// include differences between S_c and S_b ??? -- already included by mass_crit (which is now 0)
    {
        if(gsl_vector_get(mass,i) < mass_crit)
            *(params+2*S*S+6*S+S_b*N+N+S_b+2+i) = EXTINCT*pow(gsl_vector_get(mass,i),exp_suc);      // extinction boundary for small animals
        else
            *(params+2*S*S+6*S+S_b*N+N+S_b+2+i) = EXTINCT*gsl_vector_get(mass,i);                   // extinction boundary for large animals

    }

    gsl_matrix_free(A);
    gsl_matrix_free(H);

}

static void set_spatial_parameters(gsl_rng *r, gsl_matrix *SW, gsl_vector *Loc, gsl_vector *NS,  double params[])
{

    
    int i=0;
    for(int j=0; j<S*Z; j++)
    {
        for(int k=0; k<Z; k++)
        {
            *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+i) = gsl_matrix_get(SW, j, k);     //success-weighted dispersal matrix
            i++;
        }
    }



    for(int j=0; j<Z; j++)
    {
        *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+0+2*j) = gsl_vector_get(Loc, 0+j*2);  // x-coordinates
        *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+1+2*j) = gsl_vector_get(Loc, 1+j*2);  // y-coordinates
    }


    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z) = emigr_a_S_c;
    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+1) = emigr_b_S_c;
    
    
    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+2) = emigr_a_S_b;
    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+3) = emigr_b_S_b;
    
    // parameters for the RGG
    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+4) = eps;
    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+5) = D_0;
    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+6) = theta;
    
    
    // parameters for net growth rate
    for(int j=0; j<D; j++)
    {
        *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+7+j) = 0;
    }
   
    getInputNutrient(NS);
    for(int j=0; j<DN; j++)
    {
        *(params+2*S*S+7*S+S_b*N+N+S_b+2+S*Z*Z+2*Z+7+D+DN+j) = gsl_vector_get(NS,j);
    }
    
    
    return;
}


// ********** get a random number from a Gaussian distribution with specified parameters and within specified limits **********
double get_random_parameter(double mean, double sigma, double low_cutoff, double high_cutoff ,gsl_rng *r)
{
    double par = low_cutoff - 1;                                          // initialise parameter with a value outside the desired range

    while(par < low_cutoff || par > high_cutoff)
        par = gsl_ran_gaussian(r,sigma) + mean;

    return par;
}

//  ********** set the initial biomass densities and save them in an extra array **********
static void initialise_biomass_densities(double B[], double iniB[], double params[], gsl_rng *r)
{

    int i,j;
    double temp1;
    double *pparams=(double *)params;          // pointer to first element of params...



    for(i=0; i<S; i++)
    {
       
        for(j=0; j<Z; j++)
        {
           
                if(i<S_b)
                    B[j*S+i]=B_b*gsl_rng_uniform_pos(r);                    // basal species
                if(i>=S_b)
                    B[j*S+i]=B_c*gsl_rng_uniform_pos(r);                    // consumer and predator species

           

            //printf("%.9g\n",B[j*S+i]);
        }
    }

    for(i=0; i<N; i++)
    {

        for(j=0; j<Z; j++)
        {

                temp1=(double)*(params+2*S*S+5*S+2*S_b*N+i);
                B[D+j*N+i] = 0.5*temp1+0.5*temp1*gsl_rng_uniform_pos(r);    //nutrients
            }
            

        }
    



    if(KILLSWITCH)
    {

        for(j=0; j<Z; j++)
        {

            int no_species_killed = gsl_rng_uniform_int(r,(S/2+1));         // draw number of species to be killed on each patch [0,S/2]
            double temp3[no_species_killed],temp2[S];                       // initialise arrays for gsl_ran_choose

            for (i = 0; i < S; i++)
                temp2[i] = (double) i;                                      // fill array with species numbers

            gsl_ran_choose(r, temp3, no_species_killed, temp2, S, sizeof (double));    // randomly select species

//            printf("no species to be killed on patch %d: %d\n",j,no_species_killed);
//            for (i = 0; i < no_species_killed; i++)
//                printf("which species: %g\n",temp3[i]);

            for(i=0; i<S; i++)
            {
                for(int m=0; m<no_species_killed; m++)
                {
                    if(i==temp3[m])
                        B[j*S+i]=0;//0;//EXTINCT;        // initialise biomass much smaller than EXTINCT (initialising with 0 does not work -- solver stops at t=0 !!);
                }
//                    printf("%.9g\n",B[j*S+i]);
            }
//        printf("\n");
        }
    }

    for(i=0; i<D+DN; i++)
        iniB[i]=B[i];                                                       // saving initial biomass values

    return;
}


//***************************************************************************************
//  efficient environment for solving the ODE, including reducing the dimension if possible
//***************************************************************************************
static void solve_ode(double B[], double meanB[], double meanB_tot[], double CV[], double CV_tot[], double params[], gsl_matrix *CovM, double Biomass_Temp[], double net_growth_Temp[])
{
    int i,j,m,flag,ERR_FLAG,tempBm;
    double t = 0;                                                       		// starting time
    double *pparams=(double *)params;          									// pointer to first element of params...

    double *meanB2=(double *) calloc(D+DN,sizeof(double));             			// mean biomass densities (for comparison)
    gsl_vector_view Mvec2 = gsl_vector_view_array(meanB2,D+DN);       			// mean biomass denisties (including resources), for comparison

    gsl_vector_view Bvec = gsl_vector_view_array(B,D+N*Z);              		// biomass densities (including resources)
    gsl_vector_view Mvec = gsl_vector_view_array(meanB,D+N*Z);          		// mean biomass denisties (including resources)
    gsl_vector_view MTvec = gsl_vector_view_array(meanB_tot,S);         		// mean total biomass densities
    gsl_vector_view netvec = gsl_vector_view_array(net_growth_Temp,3*D);        // net growth rates temp 1&2 and their Standard deviation from 2*D onwards
    gsl_vector_view MTNvec = gsl_vector_view_array_with_stride(meanB_tot,S,N);   // mean total nutrient densities
    gsl_vector_view Bvec_temp = gsl_vector_view_array(Biomass_Temp,3*(D+N*Z));   // biomass densities (including resources)

    gsl_vector *Mvec_all = gsl_vector_calloc(Z);								// aggregate mean Biomasses of all species per patch
    gsl_vector *Nvec_all = gsl_vector_calloc(Z);								// aggregate mean nutrient level per patch
    //gsl_vector *StD_temp = gsl_vector_calloc(D);								// termporal values for standard deviations

    gsl_vector_view net_growthrate = gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+7,D); // net growth rate
    gsl_vector *netgrowthrate = &net_growthrate.vector;

    gsl_vector_view net_temp1 = gsl_vector_subvector(&netvec.vector,0,D);       // net growth rates temporally summed
    gsl_vector_view net_temp2 = gsl_vector_subvector(&netvec.vector,D,D);       // net growth rates temporally summed
    gsl_vector_view StD_all = gsl_vector_subvector(&netvec.vector,D,D);       	// Standard deviation from the net growth rate

    FILE *timeseries;
    if(TIMESERIES)
        Prepare_timeseries_file(timeseries);

    double tout = 0.1;//0.01;                                        // initial step size

    N_Vector y = N_VNew_Serial(D+DN);
    N_VSetArrayPointer(B, y);

    void *cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    CVodeSetUserData(cvode_mem,params);
    CVodeInit(cvode_mem, Cvode_rates, t, y);
    CVodeSStolerances(cvode_mem, eps_rel, eps_abs);
    CVDense(cvode_mem, D+DN);

    
    // ********* integrate ODE to reach an attractor; optional: output of time series ********
    ERR_FLAG = 0;
    
    while(t < tend)
    {
        
        if(TIMESERIES)
            flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        else
            flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

        if(flag != CV_SUCCESS)
        {
            ERR_FLAG = 1;
            break;
        }

        if(tout < Delta_t)
            tout *= 10;
        else
            tout += Delta_t;

        // Set biomass to zero under a certain threshold
        EXT_FLAG = 0;
        for (i=0; i<D+DN; i++)
        {
                if (B[i] < EXTINCT && B[i] != 0)                                     // determine whether a species just fell below the extinction threshold
                    EXT_FLAG = 1;
            //if (B[i] > 1e100)
            //EXT_FLAG = 1;
        }
        if (EXT_FLAG)                                          			// if an extinction occured:
        {
            Extinct_Species(B,S,Z);									// set the biomass of a whole species to 0...
            Extinction(B,D+DN);                                     	// set the biomass of the species to 0...
            tout = t+0.1;                                              // set the next output time to a reasonable value...
            CVodeReInit(cvode_mem, t, y);             		// and re-initialise the solver
            CVDense(cvode_mem, D+DN);
            
           
        }


       if(TIMESERIES && t >= tend - teval)
        //if(TIMESERIES)
            Write_timeseries_to_file(B,t,timeseries);
        else                                                            // calculate the coefficient of variation only if the step size is constant
        {

            if(fabs(t - (tend - 2*teval)) < 0.1*Delta_t)
                // 	copy the current biomasses to a vector for evaluation (tend-20.000)
            {
                gsl_vector_view Bvec_Temp1 = gsl_vector_subvector(&Bvec_temp.vector,0,D+N*Z);
                gsl_vector_memcpy(&Bvec_Temp1.vector,&Bvec.vector);

            }

            if(t >= tend - (1+VAR_COEFF)*teval && t < tend - teval)
            {
                gsl_vector_add(&Mvec.vector,&Bvec.vector);
                gsl_vector_add(&net_temp1.vector,netgrowthrate);
            }

            if(fabs(t - (tend - teval)) < 0.1*Delta_t)
            {
//                show_vector(&net_temp1.vector,D);

                gsl_vector_scale(&Mvec.vector,Delta_t/teval);					// scale the biomasses to one time step
                gsl_vector_scale(&net_temp1.vector,Delta_t/teval);				// scale the net growth rate to one time step

                gsl_vector_view M1vec = gsl_vector_subvector(&Mvec.vector,0,S); // mean biomass densities in first patch
                gsl_vector_view M2vec = gsl_vector_subvector(&Mvec.vector,S,S); // mean biomass densities in second patch

                gsl_vector_memcpy(&MTvec.vector,&M2vec.vector);
                gsl_vector_add(&MTvec.vector,&M1vec.vector);        // total biomasses

                gsl_vector_view M1Nvec = gsl_vector_subvector_with_stride(&Mvec.vector,0,D,N);    // mean nutrient densities in first patch
                gsl_vector_view M2Nvec = gsl_vector_subvector_with_stride(&Mvec.vector,D,N,N);   // mean nutrient densities in second patch

                gsl_vector_memcpy(&MTNvec.vector,&M2Nvec.vector);
                gsl_vector_add(&MTNvec.vector,&M1Nvec.vector);        // total nutrient densities

                // 	copy the current biomasses to a vector for evaluation (tend-10.000)
                gsl_vector_view Bvec_Temp2 = gsl_vector_subvector(&Bvec_temp.vector,D+N*Z,D+N*Z);
                gsl_vector_memcpy(&Bvec_Temp2.vector,&Bvec.vector);

                for(m=0; m<Z; m++)
                {

                    gsl_vector_view Mvec_sp = gsl_vector_subvector(&Mvec.vector,m*(S+N),S); // mean biomass densities of a single patch without resource
                    gsl_vector_view Mvec_nu = gsl_vector_subvector(&Mvec.vector,S+m*(S+N),N); // mean biomass densities of a single patch without resource


                    tempBm = 0;
                    for(i=0; i<S; i++)
                    {
                        tempBm += gsl_vector_get(&Mvec_sp.vector, i);

                    }
                    gsl_vector_set(Mvec_all, m, tempBm);				// total biomass of one patch without resource

                    tempBm = 0;
                    for(i=0; i<N; i++)
                    {
                        tempBm += gsl_vector_get(&Mvec_nu.vector, i);

                    }
                    gsl_vector_set(Nvec_all, m, tempBm);				// total nutrient level of each patch


                }

                gsl_vector_scale(Mvec_all, 1./S); // Mean species biomass of every patch		=> does this make sense?!
                gsl_vector_scale(Nvec_all, 1./N); // Mean nutrient level of every patch

                //show_vector(&Mvec.vector,D+N*Z);
            }

            if(t >= tend - teval)
                //  Calcualte mean biomass second time
            {
                
                gsl_vector_add(&Mvec2.vector,&Bvec.vector);                     // compute mean biomass densities second time
                
                //     copy the current net growth rates to a vector for evaluation (tend-10.000)
                gsl_vector_add(&net_temp2.vector,netgrowthrate);
                
            }
            
            if(fabs(t - tend) < 0.1*Delta_t)
            {
                gsl_vector_scale(&net_temp2.vector,Delta_t/teval);                // scale the net growth rate to one time step
                
                
            }
        }
    }
    
    //     copy the current biomasses to a vector for evaluation (tend-10.000)
    gsl_vector_view Bvec_Temp3 = gsl_vector_subvector(&Bvec_temp.vector,2*(D+N*Z),D+N*Z);
    gsl_vector_memcpy(&Bvec_Temp3.vector,&Bvec.vector);
    
    gsl_vector_scale(&Mvec2.vector,Delta_t/teval);                              // Scale the mean biomass vector to one time step
    
    double norm1 = gsl_blas_dasum(&Mvec.vector);
    double norm2 = gsl_blas_dasum(&Mvec2.vector);
    
    cv_flag = norm1 - norm2;
    printf("%g\n", cv_flag);

    CVodeFree(&cvode_mem);
    free(meanB2);
    gsl_vector_free(Mvec_all);
    gsl_vector_free(Nvec_all);
    

    return;


}



//***************************************************************************************
//  Set biomass densities below the extinction threshold to 0
//***************************************************************************************
static void Extinction(double B[], int Num)
{
    int i;
    for (i=0; i<Num; i++)
    {
        if (B[i] < EXTINCT)
        {
            B[i] = 0;
        }
        //printf("Local_Extinct\n");
        //if (B[i] > 1e100){
        //B[i] = 1e90;
        ////printf("super");
        //}
    }

    return;
}


//***************************************************************************************
//  Set biomass densities of whole species below the extinction threshold to 0
//***************************************************************************************
static void Extinct_Species(double B[], int Num, int Pat)
{
    int i,m;
    for (i=0; i<Num; i++)
    {
        double B_Species = 0;

        for (m=0; m<Pat; m++)
        {
            B_Species = B_Species + B[i+Num*m];
        }
        //printf("TOTBiomassSpecies: %g\n",B_Species);

        if (B_Species < GLOBAL_EXTINCT)
        {
            for (m=0; m<Pat; m++)
            {
                B[i+Num*m] = 0;
                //printf("Global_Extinct");
            }
        }
    }

    return;
}


//  ********** this function is required by the cvode-solver **********
static int Cvode_rates(double t, N_Vector y, N_Vector ydot, void *params)
{
    dynamics(N_VGetArrayPointer(y),N_VGetArrayPointer(ydot),params, t);
    return(0);
}

//  ********** helper function that prints a matrix in the standard output (e.g. a terminal) **********
static void show_matrix(gsl_matrix *M, int rows, int cols)
{

    int i, j;
    for(i=0; i<(rows); i++)
    {
        for(j=0; j<(cols); j++)
        {

            if(j==cols-1)
                printf("%.3g \n",
                       (double)gsl_matrix_get(M, i, j));
            else
                printf("%.3g \t",
                       (double)gsl_matrix_get(M, i, j));

        }
    }
    return;
}

static void show_vector(gsl_vector *A, int Num)
{
    int i,j;
    for(i = 0; i<Num; i++)
    {
        printf("%g \n",gsl_vector_get(A,i));
    }
    return;
}

//  ********** define population dynamics equations **********
static void dynamics(double B[], double Bdot[], void *params, double t)
{
    int i,j;
    double *pparams=(double *)params;                             // pointer to first element of params

    gsl_vector_view TB_vec = gsl_vector_view_array(B,D+DN);      // vector with all biomass densities and nutrient concentration
    gsl_vector *TBvec = &TB_vec.vector;

    gsl_vector_view Bdot_vec=gsl_vector_view_array(Bdot,D);     // vector for time derivatives of biomass densities
    gsl_vector *Bdotvec=&Bdot_vec.vector;

    gsl_vector_view Ndot_vec=gsl_vector_view_array(Bdot+D,DN);   // vector for time derivatives of nutrient concentrations
    gsl_vector *Ndotvec=&Ndot_vec.vector;

    gsl_vector *Ivec=gsl_vector_calloc(D);                      // ingestion: nutrient uptake, consumption, predation
    gsl_vector *Dvec=gsl_vector_calloc(D);                      // death (predation losses)
    gsl_vector *Xvec=gsl_vector_calloc(D);                      // respiration

    gsl_vector *N_in=gsl_vector_calloc(DN);                      // influx of nutrients
    gsl_vector *N_out=gsl_vector_calloc(DN);                     // outflux of nutrients
    gsl_vector *N_up=gsl_vector_calloc(DN);                      // uptake of nutrients by plants

    gsl_vector *Imvec=gsl_vector_calloc(D);                      // immigration
    gsl_vector *Emvec=gsl_vector_calloc(D);                      // emigration

    gsl_matrix_view SW_mat=gsl_matrix_view_array(pparams+2*S*S+6*S+S_b*N+N+S_b+2+S,D,Z);  // dispersal succes weighted matrix
    gsl_matrix *SW =&SW_mat.matrix;
   
    net_rates(TBvec,Ivec,Dvec,Xvec,Emvec,N_in,N_out,N_up,params,t);     // calculate net rates

//  ************* calculate dispersal success *************
    for(j=0; j<Z; j++)
    {
        for(i=0; i<S; i++)
        {
            gsl_vector_view Im_vec_i = gsl_vector_subvector_with_stride(Imvec, i, S, Z);
            gsl_vector *Imvec_i=&Im_vec_i.vector;

            gsl_vector_view Em_vec_i = gsl_vector_subvector_with_stride(Emvec, i, S, Z);
            gsl_vector *Emvec_i=&Em_vec_i.vector;

            gsl_matrix_view SW_i_mat = gsl_matrix_submatrix(SW,i*Z,0,Z,Z);
            gsl_matrix *SW_i = &SW_i_mat.matrix;

            gsl_blas_dgemv(CblasNoTrans,1,SW_i,Emvec_i,0,Imvec_i);
//		printf("s=%d, p=%d, %0.9g\n",i,j,gsl_vector_get(Imvec_i,j));

            //if(WRITENETRATESTOFILE)
            //  Write_disprates_to_file(i,j,Imvec_i,disprates); // does not work yet !!!
        }
    }


// intercept too small migrations
    for(i=0; i<D; i++)
    {
        if(gsl_vector_get(Imvec,i)<MIN_MIGRATION && gsl_vector_get(Imvec,i) != 0)
            gsl_vector_set(Imvec,i,0);
        //if(gsl_vector_get(Imvec,i)>Max_Disp)
        //gsl_vector_set(Imvec,i,Max_Disp);
    }

//  ************* assmble dynamical equations (Bdot) *************
    gsl_vector_memcpy(Ndotvec,N_in);
    
    gsl_vector_sub(Ndotvec,N_out);
    
    gsl_vector_sub(Ndotvec,N_up);
    

    gsl_vector_memcpy(Bdotvec,Ivec);
    
    gsl_vector_sub(Bdotvec,Dvec);
    
    gsl_vector_sub(Bdotvec,Xvec);
    
    gsl_vector_sub(Bdotvec,Emvec);
    
    
    for(i=0; i<D; i++)
    {
        if (B[i] < EXTINCT)
            Bdot[i] = 0;
        
    }
    gsl_vector_add(Bdotvec,Imvec);

    //  ************* control of negative biomass densities ************* ###
    //for(i=0; i<D; i++)
    //{
        //if (B[i] < EXTINCT)
           // Bdot[i] = Bdot[i]- 20*B[i];
        
   // }

    /*
    //  ************* apply extinction boundary *************
       gsl_vector_view EXT_vec=gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N+N+S_b+2,S);    // extinction boundaries (depending on body mass)
       gsl_vector *EXTvec=&EXT_vec.vector;
       for(i=0; i<D; i++){
    		if(B[i] < gsl_vector_get(EXTvec,i%S))
                Bdot[i] = Bdot[i] - 20*B[i];  // warum -20 ???
    	}
    */

//  ********** free memory ************
    gsl_vector_free(Ivec);
    gsl_vector_free(Dvec);
    gsl_vector_free(Xvec);
    gsl_vector_free(Emvec);
    gsl_vector_free(Imvec);
    gsl_vector_free(N_in);
    gsl_vector_free(N_out);
    gsl_vector_free(N_up);

    return;
}




//  ********** calculate biomass flow rates between species and resources **********
static void net_rates(gsl_vector *TBvec, gsl_vector *Ivec, gsl_vector *Dvec, gsl_vector *Xvec, gsl_vector *Emvec,
                      gsl_vector *N_in, gsl_vector *N_out, gsl_vector *N_up, void *params, double t)
{

    int i,p,k;
    double temp;
    double *pparams=(double *)params;                                     // pointer to the first element of params

    gsl_matrix_view A_mat=gsl_matrix_view_array(pparams,S,S);             // a_ij*f_ij
    gsl_matrix *Amat=&A_mat.matrix;

    gsl_matrix_view H_mat=gsl_matrix_view_array(pparams+S*S,S,S);         // H_ij
    gsl_matrix *Hmat=&H_mat.matrix;

    gsl_vector_view R_vec=gsl_vector_view_array(pparams+2*S*S,S);         // respiration rates
    gsl_vector *Rvec=&R_vec.vector;

    gsl_vector_view C_vec=gsl_vector_view_array(pparams+2*S*S+1*S,S);     // interference vector
    gsl_vector *Cvec=&C_vec.vector;

    gsl_vector_view P_vec=gsl_vector_view_array(pparams+2*S*S+2*S,S);     // basal vector (1 for basal species, 0 otherwise)
    gsl_vector *Pvec=&P_vec.vector;

    gsl_vector_view H_vec=gsl_vector_view_array(pparams+2*S*S+3*S,S);     // Hill coefficients
    gsl_vector *Hvec=&H_vec.vector;

    gsl_vector_view M_vec=gsl_vector_view_array(pparams+2*S*S+4*S,S);     // body masses
    gsl_vector *Mvec=&M_vec.vector;

    gsl_vector_view E_vec=gsl_vector_view_array(pparams+2*S*S+5*S,S);     // assimilation efficiencies
    gsl_vector *Evec=&E_vec.vector;

    gsl_matrix_view K_mat=gsl_matrix_view_array(pparams+2*S*S+6*S,S_b,N); // nutrient uptake half saturation densities
    gsl_matrix *Kmat=&K_mat.matrix;                                       // (P rows, N columns)

    gsl_vector_view S_vec=gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N,N);        // nutrient supply concentrations
    gsl_vector *Svec=&S_vec.vector;

    gsl_vector_view U_vec=gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N+N,S_b);    // max. nutrient uptake rates
    gsl_vector *Uvec=&U_vec.vector;

    gsl_vector_view Cc_vec=gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N+N+S_b,N); // nutrient contents in plants
    gsl_vector *Ccvec=&Cc_vec.vector;

    gsl_vector_view net_growthrate = gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+7,D); // net growth rate
    gsl_vector *netgrowthrate = &net_growthrate.vector;

    gsl_vector_view Nut_Syn = gsl_vector_view_array(pparams+2*S*S+7*S+S_b*N+N+S_b+2+S*Z*Z+2*Z+7+D+DN,DN); // Nutreint syn## adjust to mean
    gsl_vector *NutSyn = &Nut_Syn.vector;

    
    gsl_vector *Prey=gsl_vector_calloc(S);                                // prey biomasses

    gsl_matrix *mmat=gsl_matrix_calloc(S,S);
    gsl_vector *rvec=gsl_vector_calloc(S);
    gsl_vector *svec=gsl_vector_calloc(S);
    gsl_vector *nvec=gsl_vector_calloc(N);
    gsl_vector *mvec=gsl_vector_calloc(N);
    gsl_vector *pvec=gsl_vector_calloc(S_b);
    gsl_vector *netvec=gsl_vector_calloc(S);								// net growth rate without migration

    for(p=0; p<Z; p++)
    {

        gsl_vector_view B_vec_p	= gsl_vector_subvector(TBvec, p*S, S); 		// all biomasses (plants + animals) on patch p
        gsl_vector *Bvec_p = &B_vec_p.vector;

        gsl_vector_view PB_vec_p	= gsl_vector_subvector(TBvec, p*S, S_b); // basal biomasses (plants) on patch p
        gsl_vector *PBvec_p = &PB_vec_p.vector;

        gsl_vector_view NB_vec_p	= gsl_vector_subvector(TBvec, D+p*N, N); // nutrient biomasses on patch p
        gsl_vector *NBvec_p = &NB_vec_p.vector;

        gsl_vector_view D_vec_p	= gsl_vector_subvector(Dvec, p*S, S); 		// death rates on patch p
        gsl_vector *Dvec_p = &D_vec_p.vector;

        gsl_vector_view X_vec_p	= gsl_vector_subvector(Xvec, p*S, S); 		// respiration rates on patch p
        gsl_vector *Xvec_p = &X_vec_p.vector;

        gsl_vector_view I_vec_p	= gsl_vector_subvector(Ivec, p*S, S); 		// ingestion: nutrient uptake, consumption, predation on patch p
        gsl_vector *Ivec_p = &I_vec_p.vector;

        gsl_vector_view Em_vec_p = gsl_vector_subvector(Emvec, p*S, S); 	// emigratio rates on patch p
        gsl_vector *Emvec_p = &Em_vec_p.vector;

        gsl_vector_view N_in_p = gsl_vector_subvector(N_in, p*N, N); 	// nutrient influx on patch p
        gsl_vector *Nin_p = &N_in_p.vector;

        gsl_vector_view N_out_p = gsl_vector_subvector(N_out, p*N, N); 	// nutrient outflux on patch p
        gsl_vector *Nout_p = &N_out_p.vector;

        gsl_vector_view N_up_p	= gsl_vector_subvector(N_up, p*N, N); 	// nutrient uptake on patch p
        gsl_vector *Nup_p = &N_up_p.vector;

        gsl_vector_view Nut_Syn_p    = gsl_vector_subvector(NutSyn, p*N, N);     // nutrient syn on patch ## addjust
        gsl_vector *NutSyn_p = &Nut_Syn_p.vector;

        gsl_vector_view net_growthrate_p = gsl_vector_subvector(netgrowthrate,p*S,S); // species-specific subvector for each patch
        gsl_vector *netgrowth_p = &net_growthrate_p.vector;

//  ************* if patch is habitat then feeding etc... *************
//~ for(k = 0; k<S; k++){       // later loop over speices if species-specific differences in patch state
        k = 0;
        

//  ************ in case of Type 3 functional response ***********

            gsl_vector_memcpy(Prey,Bvec_p);

            for(i=0; i<S; i++)
            {
                if(gsl_vector_get(Prey,i)<0)
                    gsl_vector_set(Prey,i,0);                            // paranoia (this doesn't mean that you can safely remove these lines!)

                gsl_vector_set(Prey,i,pow(gsl_vector_get(Prey,i),gsl_vector_get(Hvec,i)));  // filling prey vector (prey biomasses to the power of hill coefficient)
            }

//  ************* calcualate functional responses *************

            gsl_matrix_memcpy(mmat,Hmat);
            gsl_matrix_mul_elements(mmat,Amat);

            gsl_blas_dgemv(CblasNoTrans,1,mmat,Prey,0,rvec);         // handling time term: r_i=sum_j a_ij*h_ij*Prey_j
            gsl_vector_memcpy(svec,Bvec_p);
            gsl_vector_mul(svec,Cvec);                               // predator interference: s_i=C_i*B_i
            gsl_vector_add(rvec,svec);
            gsl_vector_add_constant(rvec,1);
            gsl_vector_mul(rvec,Mvec);                               // rvec: denominator of functional response

            gsl_vector_memcpy(svec,Evec);                            // assimilation efficiency (prey-specific)
            gsl_vector_mul(svec,Prey);
            gsl_blas_dgemv(CblasNoTrans,1,Amat,svec,0,Ivec_p);       // numerator of prey intake term: I_i = sum_j e_j*A_ij*Prey_j

            gsl_vector_div(Ivec_p,rvec);                             // divide by denominator...
            gsl_vector_mul(Ivec_p,Bvec_p);                           // multiply with target species' biomass

            gsl_vector_memcpy(Xvec_p,Rvec);                          // compute the total respiration rate:
            gsl_vector_mul(Xvec_p,Bvec_p);                           // per unit biomass respiration times biomass density

            gsl_vector_memcpy(svec,Bvec_p);                          // predator biomass
            gsl_vector_div(svec,rvec);                               // ... divide by denominator of the functional response
            gsl_blas_dgemv(CblasTrans,1,Amat,svec,0,Dvec_p);         // per capita predation losses:
            gsl_vector_mul(Dvec_p,Prey);                             // multiply with target species' biomass

            //  ************* terms for basal species *************
            gsl_vector_view u_vec=gsl_vector_subvector(rvec,0,S_b);  // basal species growth rates
            gsl_vector *uvec=&u_vec.vector;

            for(i=0; i<S_b; i++)
            {
                gsl_vector_view K_vec=gsl_matrix_row(Kmat,i);
                gsl_vector *Kvec=&K_vec.vector;                        // Kvec is a vector of length N (Kvec_j=Kmat_ij)

                gsl_vector_memcpy(nvec,NBvec_p);
                gsl_vector_add(nvec,Kvec);
                gsl_vector_memcpy(mvec,NBvec_p);
                gsl_vector_div(mvec,nvec);                             // m_j=NB_j/(Kvec_j+NB_j)

                gsl_vector_set(uvec,i,gsl_vector_min(mvec));           // remember: uvec is just a vector_view of rvec!
            }

            gsl_vector_mul(uvec,Uvec);
            gsl_vector_mul(uvec,PBvec_p);                            // uvec_i=PB_i*mas_i^-0.25*Min_j(NB_j/(NB_j+Kmat_ij))

            gsl_vector_mul(rvec,Pvec);                               // growth function only for plant species
            gsl_vector_add(Ivec_p,rvec);


//  ********** put the rates for the nutrient dynamics together ********** // nutrient dynamics on non-habitat patches ???
            gsl_vector_memcpy(Nin_p,Svec);
        

            temp=gsl_vector_get(NutSyn_p,0);          //##Sin
                //printf("test %.9g\n",temp);
                gsl_vector_scale(Nin_p,temp);
        
            gsl_vector_scale(Nin_p,D_N);                              // nutrient influx: D_N * S_i
        
            gsl_vector_memcpy(Nout_p,NBvec_p);
            gsl_vector_scale(Nout_p,D_N);                             // nutrient outflux: D_N * N_i

            gsl_vector_memcpy(Nup_p,Ccvec);
            gsl_vector_scale(Nup_p,gsl_blas_dasum(uvec));             // nutrient uptake: C_i * sum_j (r_j*G_j*P_j)

        



//  ********** calculate net growth rates for basal and consumer species **********
        gsl_vector_memcpy(netvec,Ivec_p);                              // copy biomass gain by ingestion into netvec
        gsl_vector_sub(netvec,Dvec_p);                                 // substract biomass lost by death
        gsl_vector_sub(netvec,Xvec_p);                                 // substract biomass lost by respiration
        
        
        gsl_vector_div(netvec,Bvec_p);                                 // divide by biomass to get net growth rate
        
        for (i=0;i<S;i++)
        {
            if (gsl_vector_get(netvec, i)!= gsl_vector_get(netvec, i)) // Replace potential NaNs with 0
                gsl_vector_set(netvec, i, 0);
        }

//  ************* dispersal rates *************
        for(i=0; i<S_b+1; i++)		    // basal (plant) species -- passive dispersal (random) TEST DIFFERENT FUNCTIONS FOR PLANT DISPERSAL (include different dispersal directions??)
            gsl_vector_set(Emvec_p,i,(emigr_a_S_b/(1+(exp(-emigr_b_S_b*(-gsl_vector_get(Rvec, i) - gsl_vector_get(netvec, i)))))));

        for(i=S_b+1; i<S; i++)		    // consumer (animal) species -- active dispersal (body mass dependent)
            gsl_vector_set(Emvec_p,i,(emigr_a_S_c/(1+(exp(-emigr_b_S_c*(-gsl_vector_get(Rvec, i) - gsl_vector_get(netvec, i)))))));

//  ************* cp net growth rates *************
        gsl_vector_memcpy(netgrowth_p,netvec);


// if(WRITENETRATESTOFILE && gsl_vector_get(Pstatevec_p, k) == 0)
//   Write_netrates_to_file(p,netvec,Emvec_p,netrates);

        gsl_vector_mul(Emvec_p,Bvec_p);
        
    } // end loop patches


//  ********** free memory ************
    gsl_vector_free(Prey);
    gsl_vector_free(svec);
    gsl_vector_free(rvec);
    gsl_vector_free(nvec);
    gsl_vector_free(mvec);
    gsl_vector_free(netvec);
    gsl_vector_free(pvec);
    gsl_matrix_free(mmat);

    return;
}


/* The following contains functions to generate the spatial topology structures.
 * The network structures are returned in the form of weighted success matrices SW (S*Z,Z),
 * i.e., SW is non-zero if and only if patch i is linked with patch j, and the
 * actual value equals the success weighted dispersal rate of this link.
 *
 * Spatial topology generating functions:
 *  - Generate_RGG_structure - random geometric graph (RGG) build

 *
 * Requirements for spatial topology generating functions:
 *  1. pointer to an instance of the random number generator; gsl_rng *r
 *  2. double array of length S for the body masses of the species; gsl_vector *mass
 *  3. pointer to gsl_vector of length 2*Z to be filled  with the location of the patches; gsl_vector *Loc (returned)
 *  4. pointer to gsl_matrix of size S*Z,Z to be filled with success weighted dispersal rates; gsl_vector *SW (returned)
 *
 * The following also contains functions which are used within the above functions.
*/

//***************************************************************************************
//  basic RGG network structure
//***************************************************************************************
static void Generate_RGG_structure(gsl_rng *ran, gsl_vector *mass, gsl_vector *Loc, gsl_matrix *SW,gsl_matrix *SW_input,  gsl_vector *con, double params[], double RGG_params[], gsl_vector *DmaxPlant)
{

    int i,m,n,j,k;
    double X_1,X_2,Y_1,Y_2,DistX,DistY, d_mn, d_spec, sum,temp5, D_max_scale, D_max_max;
    double landscape_size = atof(landscape_s);           // read in landscape size from bash
    //landscape_size = 1;
    //if (landscape_size >= 1)                          // for Landscape sizes >= 1 and random number between 0 and 1; for Landsacpe sizes < 1 add random number between 0 an 0.1
    //{
    //stratran_D_0 = gsl_rng_uniform(ran);
    //}
    //else
    //{
    //temp11 = gsl_rng_uniform(ran);
    //stratran_D_0 = temp11/10;
    //}
    
    D_0 = 1/landscape_size;  // to scale D_0 to landscape size


    D_max_max=0.5;      // maximal distance for the largest species
    D_max_scale=D_max_max/pow(1e12,eps);  // scaling factor acording to D_max_max for D_max for maximal body size of 1e12 for consumers

// printf("%g\t%g\t%g",eps,D_0,theta);

    gsl_vector *D_Max = gsl_vector_calloc(S);               // max dispersal distance for each species
    gsl_matrix *D_spez = gsl_matrix_calloc(D,Z);            // matrix containing the specific distances d_spec=d_mn / d_max
    gsl_matrix *SU = gsl_matrix_calloc(D,Z);                // dispersal success
    gsl_matrix *W = gsl_matrix_calloc(D,Z);                 // dispersal weights
    gsl_matrix *DP = gsl_matrix_calloc(Z,Z);                // distance between patches (P_Dist)
    
    
    if(INPUT_SPP)
    {
        for (i=0; i<S_b; i++) {
            temp5=gsl_vector_get(DmaxPlant,i);
            gsl_vector_set(D_Max, i,temp5);
        }
      
    }
    else
    {
    for (i=0; i<S_b; i++)
    {

        temp5= gsl_rng_uniform(ran)*D_max_max;

        gsl_vector_set(D_Max, i, D_0*temp5);        // Calculate max dispersal distance for basal species

    }
    }
    for (i=S_b; i<S; i++)
        gsl_vector_set(D_Max, i, D_0*D_max_scale*pow(gsl_vector_get(mass,i), eps));  // Calculate max dispersal distance for consumer species

//	 printf("1 D_Max\n");
//	 show_vector(D_Max, S);

    if(INPUT_LANDSCAPE)
        getInputCoords(Loc);
    
    else{
        
        if(gsl_vector_get(Loc, 0) == 0 && gsl_vector_get(Loc, 1) == 0 && gsl_vector_get(Loc, 3) == 0)    // it's probably empty
        {
            // Generate X and Y locations for Random Geometric Graph (RGG)
            for (m=0; m<Z; m++)
            {
                gsl_vector_set(Loc, 0+m*2, gsl_rng_uniform(ran));       // Sets X location of Patch m to random value
                gsl_vector_set(Loc, 1+m*2, gsl_rng_uniform(ran));       // Sets Y location of Patch m to random value
            }
        }
    }
//	 printf("2 Loc\n");
    calc_patch_distances(DP, Loc, RGG_params); 	// Calculate distances d_mn betweeen patches with pythagoras

//    	 printf("3 d_mn\n");

    for (m=0; m<Z; m++)
    {
//			printf("4 \n");
        for (n=m+1; n<Z; n++)
        {
            // Specific Distance - Decides whether two patches are too far away for a species or not
            for (i=0; i<S; i++)
            {
                // Create submatrices
                gsl_matrix_view D_spez_i_mat = gsl_matrix_submatrix(D_spez,i*Z,0,Z,Z);
                gsl_matrix *D_spez_i = &D_spez_i_mat.matrix;
                gsl_matrix_view SU_i_mat = gsl_matrix_submatrix(SU,i*Z,0,Z,Z);
                gsl_matrix *SU_i = &SU_i_mat.matrix;
                gsl_matrix_view W_i_mat = gsl_matrix_submatrix(W,i*Z,0,Z,Z);
                gsl_matrix *W_i = &W_i_mat.matrix;

                d_spec = gsl_matrix_get(DP, m,n)/gsl_vector_get(D_Max,i);   // specific distance d_spec = d_mn/d_max

                if(d_spec <=1)   // d_max_i >= distance between patch m and n d_mn
                {
                    gsl_matrix_set(D_spez_i, m, n, d_spec);
                    gsl_matrix_set(D_spez_i, n, m, d_spec);

                    // Calculate dispersal success s_i,nm = (1-d_i,nm)^theta and save it in matrix SU_S*Z,Z
                    gsl_matrix_set(SU_i, m, n, pow((1-d_spec),theta)); // just setting the success to 1 makes solver get stuck early in simulation
                    gsl_matrix_set(SU_i, n, m, pow((1-d_spec),theta));
                } // else 0 (initial value)
            }
        }

    }
    for(i=0; i<S; i++)
    {
        gsl_matrix_view D_spez_i_mat = gsl_matrix_submatrix(D_spez,i*Z,0,Z,Z);
        gsl_matrix *D_spez_i = &D_spez_i_mat.matrix;
        gsl_matrix_view W_i_mat = gsl_matrix_submatrix(W,i*Z,0,Z,Z);
        gsl_matrix *W_i = &W_i_mat.matrix;
        gsl_matrix_view SU_i_mat = gsl_matrix_submatrix(SU,i*Z,0,Z,Z);
        gsl_matrix *SU_i = &SU_i_mat.matrix;
        gsl_matrix_view SW_i_mat = gsl_matrix_submatrix(SW,i*Z,0,Z,Z);
        gsl_matrix *SW_i = &SW_i_mat.matrix;

        //printf("Matrix W vor bearbeitung: \n");
        //show_matrix(W_i, P, P);

        for(j=0; j<Z; j++)      // über eine Zeile laufen
        {
            sum=0;              // für jede Zeile wieder bei Null anfangen
            for(k=0; k<Z; k++)  // über die Zellen der Zeile j laufen --> iteriert über Patches, die mit j verbunden sind
            {
                if(gsl_matrix_get(D_spez_i, j, k)!=0)
                {
                    sum += 1-gsl_matrix_get(D_spez_i, j, k);
                    gsl_matrix_set(W_i,j,k, (1-gsl_matrix_get(D_spez_i, j, k)));
                }
            }
            gsl_vector_view vrow = gsl_matrix_row(W_i, j);
            gsl_vector *vrow_j = &vrow.vector;

            if(sum!=0)
                gsl_vector_scale(vrow_j, 1/sum);
            // or redo and use gsl_blas_dasum (const gsl vector * x )
        }
        // Calculate success-weighted dispersal matrix elementwise SU*W
        // SW is the matrix we pass to the function (it is our weighted-success-matrix (the end results we want to attain))
        //show_matrix(W_i,Z,Z);
        //show_matrix(SU_i,Z,Z);
        gsl_matrix_memcpy(SW_i, W_i); // copy contents of W_i to SW_i
        if (EMIGRATIONLOSS) {
            gsl_matrix_mul_elements(SW_i, SU_i); // multiply all elements of W (in SW) and SU to get SW (end result)
        }
        
    }
    if (INPUT_SW) {
        getInputSW(SW_input);
        gsl_matrix_memcpy(SW, SW_input);
    }
    calc_spatial_connectance(SW,con,RGG_params);

    gsl_vector_free(D_Max);
    gsl_matrix_free(D_spez);
    gsl_matrix_free(SU);
    gsl_matrix_free(DP);
    gsl_matrix_free(W);

}


//***************************************************************************************
//  calculate the distances between patches, save in matrix
//***************************************************************************************
static void calc_patch_distances(gsl_matrix *P_Dist, gsl_vector *Loc, double RGG_params[])
{

    int m,n;
    double X_1,X_2,Y_1,Y_2,DistX,DistY, d_mn, sum_nn_dist;

    gsl_vector_view X_Vec = gsl_vector_subvector_with_stride(Loc,0,2,Z); // stride = 2, da zwischen den X immer Y Koordinaten sind
    gsl_vector *XVec=&X_Vec.vector;
    gsl_vector_view Y_Vec = gsl_vector_subvector_with_stride(Loc,1,2,Z);
    gsl_vector *YVec=&Y_Vec.vector;

    // Calculate distances d_mn betweeen patches with pythagoras
    for (m=0; m<Z; m++)
    {

        X_1 = gsl_vector_get(XVec,m);                   // take one X-value for substraction
        Y_1 = gsl_vector_get(YVec,m);                   // take one Y-value for substraction
        //printf("Patch %d x=%g y=%g", m, X_1, Y_1);

        //printf("4 \n");

        for (n=m+1; n<Z; n++)
        {
            X_2 = gsl_vector_get(XVec,n);                           // choose another X-Value to be substracted by X1
            DistX = X_1-X_2;                                        // Substraction

            Y_2 = gsl_vector_get(YVec,n);                           // choose another X-Value to be substracted by X1
            DistY = Y_1-Y_2;                                        // Substraction

            //printf("5 \n");

            // Pythagoras for Distances between two Patches
            d_mn = sqrt(pow(DistY,2)+pow(DistX,2));
            gsl_matrix_set(P_Dist, m, n, d_mn);
            gsl_matrix_set(P_Dist, n, m, d_mn);


//            printf("%g\n",d_mn);
        }
    }


    // mean nearest neighbor distance RGG_params[2]
    //alle punkte durchgehen und nächsten nachbarn finden -- mean

//       show_matrix(P_Dist,Z,Z);


    sum_nn_dist = 0;
    gsl_vector *temp = gsl_vector_calloc(Z);
    gsl_vector *temp2 = gsl_vector_calloc(Z*Z);

    for(int i=0; i<Z; i++)
    {

        gsl_vector_view P_Dist_Row = gsl_matrix_row(P_Dist,i); // create matrix view for each row i (length Z)
        gsl_vector *PDist_Row=&P_Dist_Row.vector;

        for(int j=0; j<Z; j++)
            gsl_vector_set(temp2,i*Z+j,gsl_vector_get(PDist_Row,j)); // correct ????

        gsl_vector_set(PDist_Row,i,5);       // set distance to patch itself to X to remove 0 as min dist
        sum_nn_dist += gsl_vector_min(PDist_Row);
        gsl_vector_set(temp,i,gsl_vector_min(PDist_Row)); // correct ???
//       printf("nearest neighbor dist %.9g\n",sum_nn_dist);
        gsl_vector_set(PDist_Row,i,0);       // set distance to patch itself back to 0

    }

//    show_vector(temp2,Z*Z);
//    show_matrix(P_Dist,Z,Z);
    double mean_dist = calc_mean(temp2,0,Z*Z,(Z*Z)-Z);
    double sd_nn_dist = calc_sd(temp,0,Z,Z);
    double sd_dist = calc_sd(temp2,0,Z*Z,(Z*Z)-Z);
//     printf("mean distance: %g\n",sum_dist/Z);
    RGG_params[0] = mean_dist;          // mean distance
    RGG_params[8]= sd_dist; // sd distance include into outpt !!!

    RGG_params[1]= sum_nn_dist/Z; // mean nearest neighbor distance ??
    RGG_params[9]= sd_nn_dist; // sd nearest neighbor distance ??? include into output -- output not correct yet -- finish copying into Cluster/main.cpp

//    printf("test %.9g\n",RGG_params[8]);
//    printf("test %.9g\n",RGG_params[9]);

    gsl_vector_free(temp);
    gsl_vector_free(temp2);
}

//***************************************************************************************
//  calculate spatial connectance, save in RGG_parms[]
//***************************************************************************************
static void calc_spatial_connectance(gsl_matrix *SW, gsl_vector *con, double RGG_params[])
{

    int i, j, k;
    double count2, mean_all, mean_b, mean_c, sd_all, sd_b, sd_c;

    for(i=0; i<S; i++)
    {
        gsl_matrix_view SW_i_mat = gsl_matrix_submatrix(SW,i*Z,0,Z,Z);
        gsl_matrix *SW_i = &SW_i_mat.matrix;

        count2=0;
        for(j=0; j<Z; j++)
        {
            for(k=0; k<Z; k++)
            {
                if(gsl_matrix_get(SW_i,j,k)>0) count2 ++; // number of realized links
            }
        }

        gsl_vector_set(con,i,count2/(pow(Z,2)-Z));     // calculate spatial connectance on species-level con = no.realized.links/no.possible.links (Z^2-Z)
//    printf("%g\n",count2);
//    show_matrix(SW_i,Z,Z);
    }
//show_vector(con,S);

    mean_all = calc_mean(con,0,S,S);
    mean_b = calc_mean(con,0,S_b,S_b);
    mean_c = calc_mean(con,S_b,S,S_c);

    sd_all = calc_sd(con,0,S,S);
    sd_b = calc_sd(con,0,S_b,S_b);
    sd_c = calc_sd(con,S_b,S,S_c);

    RGG_params[2] = mean_all;    // mean spatial connectance
    RGG_params[3] = sd_all;           // sd spatial connectance

    RGG_params[4] = mean_b;   // mean spatial connectance plants
    RGG_params[5] = sd_b;           // sd spatial connectance plants

    RGG_params[6] = mean_c;           // mean spatial connectance consumers
    RGG_params[7] = sd_c;           // sd spatial connectance consumers

//for (i=0;i<8;i++)
//    printf("%.9g\n",RGG_params[i]);

}

double calc_sd(gsl_vector *con, int START, int P, int div)
{
    double sum = 0.0, mean = 0, sd = 0.0;

    int i;


    for(i = START; i < P; ++i)
        sum += gsl_vector_get(con,i);

    mean = sum/div;

    for(i = START; i < P; ++i)
        sd += pow(gsl_vector_get(con,i) - mean, 2);

    return sqrt(sd/div);
}

double calc_mean(gsl_vector *con, int START, int P, int div)
{
    double sum = 0.0;
    int i;


    for(i = START; i < P; ++i)
        sum += gsl_vector_get(con,i);

    return (sum/div);
}




//***************************************************************************************
//  write the adjacency matrix, bodymasses, and biomasses to a file (in matrix/vector-format)
//***************************************************************************************
static void output(gsl_matrix *SW, gsl_vector *mass,  gsl_vector *con, double iniB[], double meanB[], double meanB_tot[], double CV[], double CV_tot[], double params[], int paramslength, double Biomass_Temp[], double RGG_params[], double net_growth_Temp[])
{
    int i,j,k,m;

    FILE *file1,*file2,*file3;

    char str[99];
    char web_out[99] = "web_";
    char mass_out[99] = "masses_";
    char global_out[99] = "global_";

    strcpy(web_out,directory);
    strcpy(mass_out,directory);
    strcpy(global_out,directory);

    strcat(web_out,"web_");
    strcat(mass_out,"mass_");
    strcat(global_out,"global_");

    strcat(web_out,name);
    strcat(mass_out,name);
    strcat(global_out,name);

    strcat(web_out,".out");
    strcat(mass_out,".out");
    strcat(global_out,".out");

    file1 = fopen(web_out,"w");
    file2 = fopen(mass_out,"w");
    file3 = fopen(global_out,"w");


    for(i=0; i<S; i++)
    {
        for(j=0; j<S; j++)
            fprintf(file1,"%g ",params[i*S+j]);
        fprintf(file1,"\n");
    }

    gsl_vector_view P_vec=gsl_vector_view_array(params+2*S*S+2*S,S);     // basal vector (1 for basal species, 0 otherwise)
    gsl_vector *Pvec=&P_vec.vector;




    // global.out
    fprintf(file3, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n","number.of.spp","number.of.plants","number.of.consumers","number.of.patch","rng.seed","max.emigr.rate.plants","shape.emigr.function.plants","max.emigr.rate.consumers","shape.emigr.function.consumers","D_0","theta","eps","mean.patch.dist","sd.patch.dist","mean.nn.dist","sd.nn.dist","mean.con.rgg","sd.con.rgg","mean.con.rgg.plants","sd.con.rgg.plants","mean.con.rgg.consumers","sd.con.rgg.consumers","Ricker");
    fprintf(file3,"%d,%d,%d,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",S,S_b,S_c,Z,seed,emigr_a_S_b,emigr_b_S_b,emigr_a_S_c,emigr_b_S_c,D_0,theta,eps,RGG_params[0],RGG_params[8],RGG_params[1],RGG_params[9],RGG_params[2],RGG_params[3],RGG_params[4],RGG_params[5],RGG_params[6],RGG_params[7],g);

    // mass.out
    fprintf(file2, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n","patch","species","body.mass", "init.biomass", "mean.biomass","biomass.variance","tot.mean.biomass","tot.biomass.variance","if.basal.spp","Biomassess_tend-20.000","Biomassess_tend-10.000","Biomassess_tend","spatial.connectance","Mean.net.growth1","Mean.net.growth2","net.growth.sd");  //local.alpha.var = CoefficientOfVariation

    for(j=0; j<Z; j++)
    {
       

        for(i=0; i<S; i++)
            fprintf(file2,"%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n", j, i, gsl_vector_get(mass,i), iniB[j*S+i],meanB[j*S+i],CV[j*S+i],meanB_tot[i],CV_tot[i],gsl_vector_get(Pvec,i),Biomass_Temp[i+j*S],Biomass_Temp[D+DN+i+j*S],Biomass_Temp[2*(D+DN)+i+j*S],gsl_vector_get(con,i),net_growth_Temp[j*S+i],net_growth_Temp[D+j*S+i],net_growth_Temp[2*D+j*S+i]);

        //nurients -- don't forget to include output !!!
        for(m=(D+j*N); m<(D+j*N+1); m++)
            fprintf(file2,"%d,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%d,%d,%d,%d,%d,%d,%d,%d\n",j, 1000, 0, iniB[m],meanB[m],CV[m],meanB_tot[m],CV_tot[m],2,0,0,0,0,0,0,0);
        for(m=(D+j*N+1); m<(D+j*N+N); m++)
            fprintf(file2,"%d,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%d,%d,%d,%d,%d,%d,%d,%d\n",j, 2000, 0, iniB[m],meanB[m],CV[m],meanB_tot[m],CV_tot[m],2,0,0,0,0,0,0,0);
    }



    if(WRITEPARAMSTOFILE)
        Write_params_to_file(params,paramslength);

//    if(WRITELANDSCAPETOFILE)
//        Write_landscape_to_file(params);
//    Write_SW_matrix_to_file(SW);

    fclose(file1);
    fclose(file2);
    fclose(file3);

    return;
}




//  ********** write the time series of the last run to an output file **********

static void Prepare_timeseries_file(FILE *timeseries)
{
    
    char str[256];
    char timeseries_out[256] = "timeseries_";
    strcpy(timeseries_out,directory);
    strcat(timeseries_out,"timeseries_");
    strcat(timeseries_out,name);
    strcat(timeseries_out,".out");

    timeseries = fopen(timeseries_out,"w");

    fprintf(timeseries,"%s \t","time");

//    for(int i=0; i < D; i++)
//        fprintf(timeseries,"%d \t",i);

    for(int m=0; m < Z; m++){
		for(int i=0; i < S; i++){
         fprintf(timeseries,"B%dP%d\t",i,m);
     }
    }

    fprintf(timeseries,"\n");
    fclose(timeseries);

    return;
}

static void Write_timeseries_to_file(double B[], double t, FILE *timeseries)
{
    //if (t> tend-teval) {
        
    
    char str[256];
    char timeseries_out[256] = "timeseries_";
    strcpy(timeseries_out,directory);
    strcat(timeseries_out,"timeseries_");
    strcat(timeseries_out,name);
    strcat(timeseries_out,".out");

    timeseries = fopen(timeseries_out,"a");

    
    fprintf(timeseries,"%g \t",t);

//    for(int i=0; i < D; i++)
//      fprintf(timeseries,"%.8g \t",B[i]);

 for(int m=0; m < Z; m++){
		for(int i=0; i < S; i++){
            
         fprintf(timeseries,"%.8g \t",B[m*S+i]);
     }
    }

    fprintf(timeseries,"\n");

    fclose(timeseries);
   // }
    return;
}

/*static void Prepare_netrates_file(FILE *netrates) {
    std::ofstream ofs;
    char str[256];
    char netvec_out[256] = "netrates_";
    strcpy(netvec_out,directory);
	strcat(netvec_out,"netrates_");
	strcat(netvec_out,name);
	strcat(netvec_out,".out");

    netrates = fopen(netvec_out,"a");
    fprintf(netrates,"%s \t %s \t %s \t %s \n","patch","species","netgrowthrate","dispersalrate");
    fclose(netrates);

    return;
}


static void Write_netrates_to_file(int p, gsl_vector *netvec,gsl_vector *Emvec_p, FILE *netrates) {
    char str[256];
    char netvec_out[256] = "netrates_";
    strcpy(netvec_out,directory);
	strcat(netvec_out,"netrates_");
	strcat(netvec_out,name);
	strcat(netvec_out,".out");

    netrates = fopen(netvec_out,"a");
    for(int i=0; i < S; i++)
        fprintf(netrates,"%d \t %d \t %.8g \t %.8g \n",p,i,gsl_vector_get(netvec,i),gsl_vector_get(Emvec_p,i));
    fclose(netrates);

    return;
}


static void Prepare_disprates_file(FILE *disprates) {
    char str[256];
    char disp_out[256] = "disprates_";
    strcpy(disp_out,directory);
	strcat(disp_out,"disprates_");
	strcat(disp_out,name);
	strcat(disp_out,".out");

    disprates = fopen(disp_out,"a");
    fprintf(disprates,"%s \ %s \t %s \t %s \n","patch","species","emigration","immigration");
    fclose(disprates);

    return;
}

static void Write_disprates_to_file(int i, int j, gsl_vector *Imvec_i, FILE *disprates) {
    char str[256];
    char disp_out[256] = "disprates_";
    strcpy(disp_out,directory);
	strcat(disp_out,"disprates_");
	strcat(disp_out,name);
	strcat(disp_out,".out");

    disprates = fopen(disp_out,"a");
    fprintf(netrates,"%d \t %d \t %.8g \t %.8g \n",j,i,gsl_vector_get(Imvec_i,j));
    fclose(disprates);

    return;
}*/

static void Write_params_to_file(double params[],int paramslength)
{
    FILE *file1;
    char str[99];
    char params_out[99] = "params_";
    strcpy(params_out,directory);
    strcat(params_out,"params_");
    strcat(params_out,name);
    strcat(params_out,".out");

    file1 = fopen(params_out,"w");

    for (int i=0; i<S*S; i++)
        fprintf(file1, "%s, %g\n", "attack rates",params[i]);
    for (int i=S*S; i<2*S*S; i++)
        fprintf(file1, "%s, %g\n", "handling times",params[i]);
    for (int i=2*S*S; i<2*S*S+S_b; i++)
        fprintf(file1, "%s, %g\n", "plant respiration",params[i]);
    for (int i=2*S*S+S_b; i<2*S*S+S; i++)
        fprintf(file1, "%s, %g\n", "animal respiration",params[i]);
    for (int i=2*S*S+S; i<2*S*S+2*S; i++)
        fprintf(file1, "%s, %g\n", "interference",params[i]);
    for (int i=2*S*S+2*S; i<2*S*S+3*S; i++)
        fprintf(file1, "%s, %g\n", "plant identifier",params[i]);
    for (int i=2*S*S+3*S; i<2*S*S+4*S; i++)
        fprintf(file1, "%s, %g\n", "hill exponent",params[i]);
    for (int i=2*S*S+4*S; i<2*S*S+5*S; i++)
        fprintf(file1, "%s, %g\n", "body masses",params[i]);
    for (int i=2*S*S+5*S; i<2*S*S+5*S+S_b; i++)
        fprintf(file1, "%s, %g\n", "plant assimilarion efficiency",params[i]);
    for (int i=2*S*S+5*S+S_b; i<2*S*S+6*S; i++)
        fprintf(file1, "%s, %g\n", "animal assimilarion efficiency",params[i]);
    for (int i=2*S*S+6*S; i<2*S*S+6*S+N*S_b; i++)
        fprintf(file1, "%s, %g\n", "nutrient uptake half sat. densities",params[i]);
    for (int i=2*S*S+6*S+N*S_b; i<2*S*S+6*S+N*S_b+N; i++)
        fprintf(file1, "%s, %g\n", "nutrient supply concentrations",params[i]);
    for (int i=2*S*S+6*S+N*S_b+N; i<2*S*S+6*S+N*S_b+N+S_b; i++)
        fprintf(file1, "%s, %g\n", "max. nutrient uptake rates",params[i]);
    for (int i=2*S*S+6*S+N*S_b+N+S_b; i<2*S*S+6*S+N*S_b+N+S_b+1; i++)
        fprintf(file1, "%s, %g\n", "content of first nutrient in plant",params[i]);
    for (int i=2*S*S+6*S+N*S_b+N+S_b+1; i<2*S*S+6*S+N*S_b+N+S_b+2; i++)
        fprintf(file1, "%s, %g\n", "content of second nutrient in plant",params[i]);
    for (int i=2*S*S+6*S+N*S_b+N+S_b+2; i<2*S*S+6*S+N*S_b+N+S_b+2+S; i++)
        fprintf(file1, "%s, %g\n", "extinction threshold",params[i]);
    for (int i=2*S*S+7*S+N*S_b+N+S_b+2; i<2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z; i++)
        fprintf(file1, "%s, %g\n", "SW dispersal matrix",params[i]);
    for (int i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z; i<2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z; i++)
        fprintf(file1, "%s, %g\n", "patch location XY",params[i]);
    int i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z;
    fprintf(file1, "%s, %g\n", "max. animal emigration rate",params[i]);
    i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+1;
    fprintf(file1, "%s, %g\n", "shape of animal emigration function",params[i]);
    i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+2;
    fprintf(file1, "%s, %g\n", "max. plant emigration rate",params[i]);
    i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+3;
    fprintf(file1, "%s, %g\n", "shape of plant emigration function",params[i]);
    i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+4;
    fprintf(file1, "%s, %g\n", "RGG epsilon",params[i]);
    i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+5;
    fprintf(file1, "%s, %g\n", "RGG dispersal distance D_0",params[i]);
    i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+6;
    fprintf(file1, "%s, %g\n", "RGG dispersal success theta",params[i]);
    for (int i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+7; i<2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+7+D; i++)
        fprintf(file1, "%s, %g\n", "net growth rate",params[i]);

    fclose(file1);

    /* To compare...
    FILE *file2;
    char t_out[99] = "testparams_";
    strcpy(t_out,directory);
    strcat(t_out,"testparams_");
    strcat(t_out,name);
    strcat(t_out,".out");
    file2 = fopen(t_out,"w");
    for (int i=0; i<paramslength;i++)
        fprintf(file2, "%g\n", params[i]);
    fclose(file2);
    */

    return;
}

static void Write_landscape_to_file(double params[])
{

    FILE *file1;
    char str[99];
    char landsc_out[99] = "landsc_";
    strcpy(landsc_out,directory);
    strcat(landsc_out,"landsc_");
    strcat(landsc_out,name);
    strcat(landsc_out,".out");

    file1 = fopen(landsc_out,"w");
    fprintf(file1, "%s, %s, %s\n","patch", "x", "y");

    for (int i=0; i<Z; i++)
        fprintf(file1, "%d, %g,%g\n", i,params[2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+S*Z+N*Z+0+2*i], params[2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+S*Z+N*Z+1+2*i]);
    fclose(file1);

    return;
}

static void Write_SW_matrix_to_file(gsl_matrix *SW)
{
    FILE *file1;
    char str[99];
    char SW_out[99] = "SW_";
    strcpy(SW_out,directory);
    strcat(SW_out,"SW_");
    strcat(SW_out,name);
    strcat(SW_out,".out");

    file1 = fopen(SW_out,"w");

    for(int i=0; i<S*Z; i++)
    {
        for(int j=0; j<Z; j++)
        {

            if(j==(Z-1))
                fprintf(file1,"%g \n", gsl_matrix_get(SW, i, j));
            else
                fprintf(file1,"%g \t", gsl_matrix_get(SW, i, j));
        }

    }

    fclose(file1);

    return;
}


static void getInputSpp()
{
    
    //std::string s = std::to_string(Web_Id);
    //char const* s=Web_Id;
    //char const *pchar = s.c_str();  //use char const* as target type
    char mass_in[99] = "/BodyMass_";
    strcpy(mass_in,inputdir);
    strcat(mass_in,"/BodyMass_");
    strcat(mass_in,Web_Id);
    strcat(mass_in,".out");
    //~ printf("%s",mass_in);
    
    std::ifstream input(mass_in); //put your program together with this file in the input directory
    std::vector<double> spp;
    std::string line;
    
    if(!input.is_open())
    {
        std::cerr << "There was a problem opening the input file1\n";
        //printf("%s",mass_in);
        std::cerr << mass_in;
        exit(1);
    }
    
    std::string col2; // needed for comma-seperator
    double col1,col3;
    
    while(getline(input,line))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2 >> col3;
        spp.push_back(col3);
    }
    
    S_b=0;
    S_c=0;
    S=0;
    
    
    for(unsigned i=1; i<spp.size(); i++)
    {
        if(spp[i]>0)
        S_b += 1;
        else
        S_c += 1;
    }
    
    S = S_b+S_c;
    
    input.close();  // close input file
    //~ printf("S= %d, S_b = %d, S_c = %d\n",S,S_b,S_c);
    
    return;
}

static void getInputRicker()
{
    
    //std::string s = std::to_string(Web_Id);
    //char const* s=Web_Id;
    //char const *pchar = s.c_str();  //use char const* as target type
    char rick_in[99] = "/params_";
    strcpy(rick_in,inputdir);
    strcat(rick_in,"/params_");
    strcat(rick_in,Web_Id);
    strcat(rick_in,".out");
    //~ printf("%s",rick_in);
    
    std::ifstream input(rick_in); //put your program together with this file in the input directory
    std::vector<double> rick;
    std::string line;
    
    if(!input.is_open())
    {
        std::cerr << "There was a problem opening the input file2\n";
        //printf("%s",rick_in);
        std::cerr << rick_in;
        exit(1);
    }
    
    std::string col2, col4; // needed for comma-seperator
    double col1,col3, col5;
    
    while(getline(input,line))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2 >> col3 >> col4 >> col5;
        rick.push_back(col5);
    }
    
    for(unsigned i=0; i<1;i++)
    {
        g=rick[i+1];
        
        
    }
    

    
    input.close();  // close input file
    //printf("rick=%g.9\n",g);
    
    return;
}


static void getInputMass(gsl_vector *mass, gsl_matrix *Up,gsl_vector *DmaxPlant)
{
    
    //std::string s = std::to_string(Web_Id);
    //char const* s=Web_Id;
    //char const *pchar = s.c_str();  //use char const* as target type
    char mass_in[99] = "/BodyMass_";
    strcpy(mass_in,inputdir);
    strcat(mass_in,"/BodyMass_");
    strcat(mass_in,Web_Id);
    strcat(mass_in,".out");
    //printf(mass_in);
    //~ file1 = fopen(mass_in,"r");
    
    std::ifstream input(mass_in); //put your program together with this file in the input directory
    std::vector<double> masses;
    std::vector<double> dispersaldist;
    std::vector<double> NutUptake1;
    std::string line;
    
    if(!input.is_open())
    {
        std::cerr << "There was a problem opening the input file3\n";
        //printf("%s",mass_in);
        std::cerr << mass_in;
        exit(1);
    }
    
    std::string col2, col4, col6; // needed for comma-seperator
    double col1,col3, col5, col7;
    
    while(getline(input,line))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7;
        masses.push_back(col1);
        dispersaldist.push_back(col5);
        NutUptake1.push_back(col7);
        
    }
    
    for(unsigned i=0; i<(masses.size()-1); i++)
    {
        gsl_vector_set(mass,i,masses[i+1]);
    }
    
    for(unsigned i=0; i<S_b;i++)
    {
        gsl_vector_set(DmaxPlant,i,dispersaldist[i+1]);
        
        gsl_matrix_set(Up,i,0,NutUptake1[i+1]);
        
        
        
    }
    input.close();  // close input file
    //~ show_vector(mass,S);
    
    return;
}

static void getInputLandscape()
{
    
    //std::string s = std::to_string(RGG_Id);
    //char const* s= RGG_Id;
    //char const *pchar = s.c_str();  //use char const* as target type
    char rgg_in[99] = "/landscape_";
    strcpy(rgg_in,inputdir);
    strcat(rgg_in,"/landscape_");
    strcat(rgg_in,RGG_Id);
    strcat(rgg_in,".out");
    //printf("%s",rgg_in);
    
    std::ifstream input(rgg_in); //put your program together with this file in the input directory
    std::vector<double> patch_z;
    std::string line;
    
    if(!input.is_open())
    {
        std::cerr << "There was a problem opening the RGG input file4\n";
        exit(1);
    }
    
    std::string col2; // needed for comma-seperator
    std::string col4; // needed for comma-seperator
    double col1,col3, col5;
    
    while(getline(input,line))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2 >> col3 >> col4 >> col5;
        patch_z.push_back(col1);
    }
    
    Z = 0;
    
    for(unsigned i=1; i<patch_z.size(); i++)
    {
        Z += 1;
    }
    
    input.close();  // close input file
    //printf("Z = %d", Z);
    
    return;
}


static void getInputCoords(gsl_vector *coords_in_v)
{
    
    //std::string s = std::to_string(RGG_Id);
    //char const* s= RGG_Id;
    //char const *pchar = s.c_str();  //use char const* as target type
    char rgg_in[99] = "/landscape_";
    strcpy(rgg_in,inputdir);
    strcat(rgg_in,"/landscape_");
    strcat(rgg_in,RGG_Id);
    strcat(rgg_in,".out");
    //~ printf("%s",rgg_in);
    
    std::ifstream input(rgg_in); //put your program together with this file in the input directory
    std::vector<double> coords_in_X;
    std::vector<double> coords_in_Y;
    std::string line;
    
    if(!input.is_open())
    {
        std::cerr << "There was a problem opening the RGG input file5\n";
        exit(1);
    }
    
    std::string col2; // needed for comma-seperator
    std::string col4; // needed for comma-seperator
    double col1,col3, col5;
    
    while(getline(input,line))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2 >> col3 >> col4 >> col5;
        coords_in_X.push_back(col3);
        coords_in_Y.push_back(col5);
    }
    
    
    int j=0;
    
    for(int i=0; i<Z; i++)
    {
        gsl_vector_set(coords_in_v,0+j*2,coords_in_X[i]);
        gsl_vector_set(coords_in_v,1+j*2,coords_in_Y[i]);
        j+=1;
    }
    
    input.close();  // close input file
    //~ show_vector(coords_in_v,Z*2);
    //~ printf("/n");
    
    return;
}

static void getInputSW(gsl_matrix *SW_input)
{


    char sw_in[99] = "/SW";
    strcpy(sw_in,inputdir);
    strcat(sw_in,"/SW");
    strcat(sw_in,RGG_Id);
    strcat(sw_in,"_");
    strcat(sw_in,Disp);
    strcat(sw_in,".out");
    
    int i, j;


    std::ifstream in(sw_in);
    double d;

    if (!in) {
        std::cerr << "There was a problem opening the SW input file\n";
        exit(1);
    }
    
    
    
    for(i=0; i<S*Z; i++)
    {
        for(j=0; j<Z; j++)
        {
            in >> d;

            gsl_matrix_set(SW_input,i,j,d);
        }
    
    }

    in.close();

    return;
}

static void getInputNutrient(gsl_vector *NS_input)
{
    
    
    char sw_in[99] = "/Nutrient_";
    strcpy(sw_in,inputdir);
    strcat(sw_in,"/Nutrient_");
    strcat(sw_in,NUT);
    strcat(sw_in,".out");
    
    int i, j;
    
    
    std::ifstream in(sw_in);
    double d;
    
    if (!in) {
        std::cerr << "There was a problem opening the Nutrient input file\n";
        exit(1);
    }
    
    
    
  
        for(j=0; j<Z; j++)
        {
            in >> d;
            
            gsl_vector_set(NS_input,j,d);

        }
        
    
    
    in.close();
    
    return;
}

