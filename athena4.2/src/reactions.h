#include "copyright.h"
#ifndef REACTIONS_H
#define REACTIONS_H
/*******************************************************************************
 * FILE: reactions.h
 *
 * PURPOSE: Global variables used by the routines describing reactive processes.
 *******************************************************************************/

#ifdef REACTION

#include "athena.h"

/*-----------------------------------------------------------------------------*/
/* definitions included in reactions.c */
#ifdef REACTIONS_C

/* Flag indicating whether reaction is activated */
int  iReactionActive;

/* Maximum acceptable relative change in internal energy or individual species
   concentration used to limit the next time step. */
Real maxVariation;
/* Suggested next time step based on the reactions. */
Real dtMinReac;

/* The definition of "fuel" depends on the reaction kinetics - see reactions.c.*/
// Real Mburn=0.0;          /* Accumulator of the total burned fuel mass [g]      */
// Real FCR=0.0;            /* Global fuel consumption rate [g/s] 				  */                
// Real Eburn=0.0;          /* Accumulator of the total energy released [erg]     */
// Real EburnRate=0.0;      /* Global burning rate at the current step [erg/s]    */

#ifdef REACTION_INFINITELY_FAST
int NSP=5;
#endif

/*-----------------------------------------------------------------------------*/
/* definitions included everywhere except reactions.c  */

#else /* REACTIONS_C */

extern int  iReactionActive;
extern Real maxVariation;
extern Real dtMinReac;

// extern Real Mburn;
// extern Real FCR;
// extern Real Eburn;
// extern Real EburnRate;

#ifdef REACTION_INFINITELY_FAST
extern int NSP;
#endif


#endif /* REACTIONS_C */
#endif /* REACTION */
#endif /* REACTIONS_H */
