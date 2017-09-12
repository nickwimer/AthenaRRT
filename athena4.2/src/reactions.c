#include "copyright.h"
#define REACTIONS_C
/*=============================================================================
 * FILE: reaction.c
 *
 * PURPOSE: Contains functions that provide description of reactive processes.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *    reaction_init()      - sets pointer to an appropriate reaction function
 *    reaction_infinitely_fast() - one-step infinitely_fast chemistry
 *=============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "eos.h"
#include "reactions.h"

#if defined(REACTION) && defined(ISOTHERMAL)
#error Isothermal gas requested together with reactions
#endif

static Real CellVolume;

/*-----------------------------------------------------------------------------*/
/* reaction_init: initialize pointer to an appropriate reaction function and
 *                call the initialization function for the given reaction type
 *                and read (or set to default values) various parameters
 *                controlling the integration.
 *
 * VGDFun_t is a function of type void which takes a Grid and a Domain as 
 * arguments.
 */

VGFun_t reaction_init(GridS *pG, DomainS *pD)
{
  Real dx1 = (pG->Nx1 > 1) ? pG->dx1 : 1.0;
  Real dx2 = (pG->Nx2 > 1) ? pG->dx2 : 1.0;
  Real dx3 = (pG->Nx3 > 1) ? pG->dx3 : 1.0;
  CellVolume = dx1 * dx2 * dx3;

#ifdef REACTION
  maxVariation = par_getd_def("reactions","maxVariation",5.e-2);

#if defined(REACTION_INFINITELY_FAST)
  return reaction_infinitely_fast;
#else
  #error("No reaction model is defined...")
#endif

#endif /* REACTION */

  return NULL;
}

#ifdef REACTION_INFINITELY_FAST
void reaction_infinitely_fast(GridS *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int n, ierror=0, load;
  Real rho, invrho, Ekin, rhoes, Ynow;
  Real xtot, Mtot;
  Real y[NSP], ySave[NSP];
  Real x[NSP], M[NSP], dE, dY;
  Real EburnNow, MburnNow, ERatio, yRatio, Ereac;

  /* Check the reaction activation flag. */
  if (!iReactionActive) {
    dtMinReac = HUGE_NUMBER;
    return;
  }

  EburnNow = 0.0;
  MburnNow = 0.0;
  ERatio   = -HUGE_NUMBER;
  yRatio   = -HUGE_NUMBER;

  /* define molar mass of species (hard coded for CH4)*/
  M[0] = 16.04;   //CH4 [g/mol]
  M[1] = 31.9988;  //O2
  M[2] = 44.01;   //CO2
  M[3] = 18.015;  //H2O
  M[4] = 28.013;  //N2

  /* Loop over the entire grid. */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {

        /* Prepare input */
        rho    = pG->U[k][j][i].d;
        invrho = 1.0 / rho;
        Ekin   = 0.5*(pG->U[k][j][i].M1*pG->U[k][j][i].M1 + pG->U[k][j][i].M2*pG->U[k][j][i].M2 + pG->U[k][j][i].M3*pG->U[k][j][i].M3)*invrho;
        rhoes  = pG->U[k][j][i].E - Ekin;

        xtot = 0.0;
        Mtot = 0.0;
        for (n = 0; n < NSP; ++n) {
    /* Species mass fractions */
          y[n]      = pG->U[k][j][i].s[n] * invrho;
          ySave[n]  = y[n];
          xtot      += (y[n]/M[n]);
        }
    /* Species mole fraction */
        for (n = 0; n < NSP; ++n) {
          x[n]      = (y[n]/M[n])/(xtot);
          Mtot      += x[n] * M[n];
        }
        Ereac = 0.0;
    /* Determine limiting reaction and react species */
        if (x[0] > 2*x[1]) { //O2 is the limiting reactant
          // if (x[1] > 0.1) {
            x[0] = x[0] - 0.5*x[1]; // subtract CH4
            x[2] = x[2] + 0.5*x[1]; // add CO2
            x[3] = x[3] + x[1];     // add H2O
            Ereac= -8.91e12*(0.5*x[1]/Mtot)/4.0;
            x[1] = 0.0;     // subtract O2
          // }
        } else { // CH4 is limiting or they are perfectly matched
          // if (x[0] > 0.1) {
            x[1] = x[1] - 2.0*x[0]; // subtract O2
            x[2] = x[2] + x[0];     // add CO2
            x[3] = x[3] + 2.0*x[0]; // add H2O
            Ereac= -8.91e12*x[0]/Mtot/4.0;
            x[0] = 0.0;     // subtract CH4
          // }
        }
    /* Convert mole fraction back to mass fraction */
        Mtot = 0.0;
        for (n = 0; n < NSP; ++n) {
          Mtot   += x[n] * M[n]; 
        }
        for (n = 0; n < NSP; ++n) {
          y[n]  = x[n] * M[n]/Mtot;
        }

        /* Sum up the fuel consumption rate. Index of species that defines
           "fuel" is 0 for CH4 */
        // MburnNow += pG->U[k][j][i].s[0] - rho * y[0];
        /* Compute changes in mass fractions and energy.
           Update densities of reacting species. */
        dE = Ereac;
        for (n = 0; n < NSP; ++n) {
          dY   = y[n] - ySave[n];
          dY   = fabs(dY);
          Ynow = MAX(0.0,ySave[n]);
          if (dY > yRatio * Ynow) yRatio = dY * (1.0 / Ynow);
    /* Update the Grid composition. Prevent extremely small values, which
       can cause overflow when Grid is saved and later visualized
       offline. */
          if (y[n] < TINY_NUMBER)
            pG->U[k][j][i].s[n] = 0.0;
          else
            pG->U[k][j][i].s[n] = rho*y[n];
        }
        dE *= rho;

        /* Update due to change in sensible internal energy density */
        pG->U[k][j][i].E -= dE;

        /* Sum up the energy release rate. We flip the sign here since
           dY used above to calculate dE is negative since fuel is mostly
     consumed. As a result, dE is also mostly negative. */
        // EburnNow -= dE;

        dE = fabs(dE);
        if (dE > ERatio * rhoes) ERatio = dE * (1.0 / rhoes);

      }
    }
  }

  /* Note, the global quantities returned by this routine are NOT summed over
     processors. Note, also that we need to divide by dt here, since between
     now and the next invocation of dump_history(), timestep will change. */
  /* Multiply global quantities by cell volume. */
  /* Total mass of fuel burned during the current step. */
  // MburnNow *= CellVolume;
  /* Total current fuel consumption rate. */
  // FCR       = MburnNow / pG->dt;
  /* Total burning energy released during the current step. */
  // EburnNow *= CellVolume;
  /* Total current burning rate. Total burning rate is not dE/dt, but the actual
     Ereleased/timestep. In the limit of vanishingly small timesteps they should
     be equal, however in reality it is more accurate to find the burning rate
     estimate based on the actual amount of energy released. */
  // EburnRate = EburnNow / pG->dt;
  /* Finally add current MburnNow and EburnNow to the total Mburn and Eburn. */
  // Mburn += MburnNow;
  // Eburn += EburnNow;

  /* Finally, calculate the next time step. The time step is limited to ensure
     that maximum relative change in fuel mass fraction or internal energy
     does not exceed the specified tolerance, typically 5% per time step. */
  if ((yRatio == -HUGE_NUMBER) && (ERatio == -HUGE_NUMBER)) {
    /* Nothing burned in this grid. */
    dtMinReac = HUGE_NUMBER;
  } else {
    dtMinReac = pG->dt*maxVariation/MAX(yRatio,ERatio);
  }
}
#endif /* REACTION_INFINITELY_FAST */
