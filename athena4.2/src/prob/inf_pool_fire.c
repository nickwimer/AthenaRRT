// #include "copyright.h"
/*==============================================================================
 // * FILE: Seeping_Flow.c
 *
 * PURPOSE: 
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "transport.h"
#include "reactions.h"
#include "constants.h"
#include "eos.h"
#include "config.h"

static Real rho;
static Real P0;
static Real lx_burner;
static Real R;
static Real v_in;
static Real v_co;
static Real TFuel;
static Real grav;
static Real XFuel[5];
static Real XAir[5];
static Real XMix[5];
static Real h;
static Real w;
static Real nu_iso;
static Real rhoAir;
static Real rhoFuel;

static Real grav_pot(const Real x1, const Real x2, const Real x3);

#define NO  0
#define YES 1

void walls_ijb(Grid *pGrid);
void hydro_ojb(Grid *pGrid);
void farfield_iib(Grid *pGrid);
void farfield_oib(Grid *pGrid);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(GridS *pGrid, DomainS *pDomain)
{

  int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;
  Real pg;                              /* Hydrostatic Pressure */
  Real M1, M2, M3, E;                   /* Temporary variables to shorten equations */
  int n;
  Real x1, x2, x3;
  Real Y[NSP], rhoX[NSP];
  Real Eterm, EFuel, RFuel;

/* Read problem parameters */
  P0               = par_getd("problem","P0");
  lx_burner        = par_getd("problem","lx_burner");
  v_in             = par_getd("problem","v_in");
  v_co             = par_getd("problem","v_co");
  TFuel            = par_getd("problem","TFuel");
  R                = par_getd("problem","R");
  grav             = par_getd("problem","grav");
  h                = par_getd("grid","x2max");
  w                = par_getd("grid","x1max");
  Gamma            = par_getd("problem","gamma");
  nu_iso           = par_getd("problem","nu_iso");


  // rhoAir = P0/(R*300.0);
  rhoAir = 1.225e-3;
  RFuel = 5.20e6;
  rhoFuel = P0/(RFuel*TFuel);

  ShearViscCoef0  = nu_iso*rhoAir; //make CGS compatible?
  BulkViscCoef0   = nu_iso*rhoAir;
  expTempShearVisc= 0.7e0;
  expTempBulkVisc = 0.7e0;
  Gamma_1 = Gamma - 1.0;
  Gamma_2 = Gamma - 2.0;

  // NSP = 5;

/* Set the diffusion coefficient for the different species */
  for (n=0; n<NSP; ++n) {
    DiffCoef0[n]  = 6.25e-6; // this is for methane, need to find the values for N2, O2, CO2, and H2O
  }

/* Set the fuel and air mass fraction vectors */
  for (n=0; n<NSP; ++n) {
    XFuel[n] = 0.0;
  }
  XFuel[0] = 1.0;

  for (n=0; n<NSP; ++n) {
    XAir[n] = 0.0;
  }
  XAir[1] = 0.233;
  XAir[4] = 1.0 - XAir[1];

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        pGrid->U[k][j][i].d  = rhoAir;
        pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
        for (n=0; n<NSP; ++n) {
          pGrid->U[k][j][i].s[n] = XAir[n]*rhoAir;
        }
        pg  =  grav*rhoAir*x2 + P0;
        pGrid->U[k][j][i].E = pg/Gamma_1;
      }
    }
  }

  // /* Enroll the gravitational potential function */
  StaticGravPot = grav_pot;
  // /* Enroll the custom boundary conditions */
  set_bvals_mhd_fun(left_x2,walls_ijb);
  set_bvals_mhd_fun(right_x2,hydro_ojb);
  set_bvals_mhd_fun(left_x1,farfield_iib);
  set_bvals_mhd_fun(right_x1,farfield_oib);

  iReactionActive = YES;

  return;

}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(GridS *pG, DomainS *pD, FILE *fp)
{
  return;
}

void problem_read_restart(GridS *pG, DomainS *pD, FILE *fp)
{
  int i, is=pG->is, ie = pG->ie;
  int j, js=pG->js, je = pG->je;
  int k, ks=pG->ks, ke = pG->ke;
  Real pg;                              /* Hydrostatic Pressure */
  Real M1, M2, M3, E;                   /* Temporary variables to shorten equations */
  int n;
  Real x1, x2, x3;
  Real Y[NSP], rhoX[NSP];
  Real Eterm, EFuel, RFuel;

/* Read problem parameters */
  P0               = par_getd("problem","P0");
  lx_burner        = par_getd("problem","lx_burner");
  v_in             = par_getd("problem","v_in");
  v_co             = par_getd("problem","v_co");
  TFuel            = par_getd("problem","TFuel");
  R                = par_getd("problem","R");
  grav             = par_getd("problem","grav");
  h                = par_getd("grid","x2max");
  w                = par_getd("grid","x1max");
  Gamma            = par_getd("problem","gamma");
  nu_iso           = par_getd("problem","nu_iso");


  // rhoAir = P0/(R*300.0);
  rhoAir = 1.225e-3;
  RFuel = 5.20e6;
  rhoFuel = P0/(RFuel*TFuel);

  ShearViscCoef0  = nu_iso*rhoAir; //make CGS compatible?
  BulkViscCoef0   = nu_iso*rhoAir;
  expTempShearVisc= 0.7e0;
  expTempBulkVisc = 0.7e0;
  Gamma_1 = Gamma - 1.0;
  Gamma_2 = Gamma - 2.0;

  // NSP = 5;

/* Set the diffusion coefficient for the different species */
  for (n=0; n<NSP; ++n) {
    DiffCoef0[n]  = 6.25e-6; // this is for methane, need to find the values for N2, O2, CO2, and H2O
  }

/* Set the fuel and air mass fraction vectors */
  for (n=0; n<NSP; ++n) {
    XFuel[n] = 0.0;
  }
  XFuel[0] = 1.0;

  for (n=0; n<NSP; ++n) {
    XAir[n] = 0.0;
  }
  XAir[1] = 0.233;
  XAir[4] = 1.0 - XAir[1];

  // /* Enroll the gravitational potential function */
  StaticGravPot = grav_pot;
  // /* Enroll the custom boundary conditions */
  set_bvals_mhd_fun(left_x2,walls_ijb);
  set_bvals_mhd_fun(right_x2,hydro_ojb);
  set_bvals_mhd_fun(left_x1,farfield_iib);
  set_bvals_mhd_fun(right_x1,farfield_oib);

  iReactionActive = YES;

  return;

}

Gasfun_t get_usr_expr(const char *expr){ 
//   if (strcmp(expr,"T")==0) return calc_T;
  return NULL; 
  }

VGFunout_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(GridS *pGrid, DomainS *pDomain)
{
  return;
}

void Userwork_after_loop(GridS *pGrid, DomainS *pDomain)
{
  return;
}

/*------------------------------------------------------------------------------
 * grav_pot: Gravitational potential
 */

static Real grav_pot(const Real x1, const Real x2, const Real x3)
{
  return -grav*x2;
}

void walls_ijb(GridS *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
  int n;
  Real rhod;
  Real M1, M2, M3, pres, E;
  Real pg;
  Real x1,x2,x3;
  Real sburn, eburn;

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  sburn=  (w-lx_burner)/2.0;
  eburn=  sburn + lx_burner;

  // rhoFuel = 0.65e-3;
  // GFuel = 1.32;

  /* Set everthing to wall boundary conditions except for length of burner (inlet)*/
for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pGrid,i,js-j,k,&x1,&x2,&x3);
        if (x1>sburn & x1<eburn) {
          pg = rhoFuel*grav*x2 + P0;
          pGrid->U[k][js-j][i].d  =  rhoFuel;
          pGrid->U[k][js-j][i].M1 =  0.0;
          pGrid->U[k][js-j][i].M2 =  rhoFuel*v_in;
          pGrid->U[k][js-j][i].M3 =  0.0;
          for (n=0; n<NSP; ++n) {
            pGrid->U[k][js-j][i].s[n] = XFuel[n]*rhoFuel;
          }
          pGrid->U[k][js-j][i].E  =  pg/Gamma_1 + 0.5*rhoFuel*(v_in*v_in);
        }
        else { //Air coflow region (make sure to change spongelayer to encourage coflow instead of ambient)
          pg = rhoAir*grav*x2 + P0;
          pGrid->U[k][js-j][i].d  =  rhoAir;
          pGrid->U[k][js-j][i].M1 =  0.0;
          pGrid->U[k][js-j][i].M2 =  rhoAir*v_co;
          pGrid->U[k][js-j][i].M3 =  0.0;
          for (n=0; n<NSP; ++n) {
            pGrid->U[k][js-j][i].s[n] = XAir[n]*rhoAir;
          }
          pGrid->U[k][js-j][i].E  =  pg/Gamma_1 + 0.5*rhoAir*(v_co*v_co);
        }
      }
    }
  }

  return;
}

void hydro_ojb(GridS *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
  int n;
  Real rhod, E, M1, M2, M3;
  Real pg;
  Real x1, x2, x3;
  Real U1, U2, U3;

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pGrid,i,je+j,k,&x1,&x2,&x3);
        rhod  = pGrid->U[k][je][i].d;
        M1    = pGrid->U[k][je][i].M1;
        M2    = pGrid->U[k][je][i].M2;
        M3    = pGrid->U[k][je][i].M3;
        U1    = M1/rhod;
        U2    = M2/rhod;
        U3    = M3/rhod;
        pg    = grav*rhoAir*x2 + P0;
        if (U2 >= 0.0) {  // outflow condition -- let pressure be total and zerogradient velocity
          pGrid->U[k][je+j][i]    = pGrid->U[k][je][i];
          pGrid->U[k][je+j][i].E  = pg/Gamma_1 + 0.5*rhod*(U1*U1 + U2*U2 + U3*U3);
        } else {  //  inflow condition -- modify pressure with dynamic portion and set tangential velocity components to zero
          pGrid->U[k][je+j][i].d  = rhoAir;
          pGrid->U[k][je+j][i].M1 = 0.0;
          pGrid->U[k][je+j][i].M2 = rhoAir*U2; 
          pGrid->U[k][je+j][i].M3 = 0.0; 
          U1    = pGrid->U[k][je+j][i].M1/rhoAir;
          U2    = pGrid->U[k][je+j][i].M2/rhoAir;
          U3    = pGrid->U[k][je+j][i].M3/rhoAir;
          for (n=0; n<NSP; ++n) {
            pGrid->U[k][je+j][i].s[n] = XAir[n]*rhoAir;
          }
          pGrid->U[k][je+j][i].E  = pg/Gamma_1 + 0.5*rhoAir*(U1*U1 + U2*U2 + U3*U3);
        }
      }
    }
  }
  return;
}

void farfield_iib(GridS *pGrid) // this is no longer a hydrostatic boundary and needs to be a coflow boundary
{
  int is = pGrid->is;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,jl,ju; /* i-lower/upper */
  int n;
  Real rhod, E, M1, M2, M3;
  Real x1, x2, x3;
  Real pg;
  Real dudx, dudy, v, u, sig, K, eta;
  Real U1, U2, U3;

  ju = pGrid->je + nghost;
  jl = pGrid->js - nghost;

  for (k=ks; k<=ke; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        cc_pos(pGrid,is-i,j,k,&x1,&x2,&x3);
        rhod  = pGrid->U[k][j][is].d;
        // Pressure boundary condition only
        M1    = pGrid->U[k][j][is].M1;
        M2    = pGrid->U[k][j][is].M2;
        M3    = pGrid->U[k][j][is].M3;
        U1    = M1/rhod;
        U2    = M2/rhod;
        U3    = M3/rhod;
        pg  =  grav*rhoAir*x2 + P0;

        pGrid->U[k][j][is-i].d  = rhoAir;
        pGrid->U[k][j][is-i].M1 = 0.0;
        // pGrid->U[k][j][is-i].M1 = rhoAir*U1;
        pGrid->U[k][j][is-i].M2 = rhoAir*v_co;
        pGrid->U[k][j][is-i].M3 = 0.0;
        for (n=0; n<NSP; ++n) {
          pGrid->U[k][j][is-i].s[n] = XAir[n]*rhoAir;
        }
        pGrid->U[k][j][is-i].E  = pg/Gamma_1 + 0.5*rhoAir*v_co*v_co;
      }
    }
  }
  return;
}

void farfield_oib(GridS *pGrid)
{
  int ie = pGrid->ie;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,jl,ju; /* i-lower/upper */
  int n;
  Real rhod, E, M1, M2, M3;
  Real x1, x2, x3;
  Real pg;
  Real dudx, dudy, v, u, sig, K, eta;
  Real U1, U2, U3;

  ju = pGrid->je + nghost;
  jl = pGrid->js - nghost;

  for (k=ks; k<=ke; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        cc_pos(pGrid,ie+i,j,k,&x1,&x2,&x3);
        rhod  = pGrid->U[k][j][ie].d;
        // Pressure boundary condition only
        M1    = pGrid->U[k][j][ie].M1;
        M2    = pGrid->U[k][j][ie].M2;
        M3    = pGrid->U[k][j][ie].M3;
        U1    = M1/rhod;
        U2    = M2/rhod;
        U3    = M3/rhod;
        pg  =  grav*rhoAir*x2 + P0;

        pGrid->U[k][j][ie+i].d  = rhoAir;
        pGrid->U[k][j][ie+i].M1 = 0.0;
        // pGrid->U[k][j][ie+i].M1 = rhoAir*U1;
        pGrid->U[k][j][ie+i].M2 = rhoAir*v_co;
        pGrid->U[k][j][ie+i].M3 = 0.0;
        for (n=0; n<NSP; ++n) {
          pGrid->U[k][j][ie+i].s[n] = XAir[n]*rhoAir;
        }
        pGrid->U[k][j][ie+i].E  = pg/Gamma_1 + 0.5*rhoAir*v_co*v_co;
      }
    }
  }
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
// 
