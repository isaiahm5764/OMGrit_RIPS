/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. For support, post issues to the XBraid Github page.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 ***********************************************************************EHEADER*/

/**
 * Example:       ex-01-optimiztion.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-01-optimization
 *
 * Help with:     This implements a simple optimization iteration for the 
 *                simple ODE example in ex-01 with gradient computation 
 *                in ex-01-adjoint. Read the source for more information.
 *
 * Sample run:    mpirun -np 2 ex-01-optimization
 *
 * Description:   Solves the simple optimization problem 
 * 
 *                   min J
 *                   s.t. u' = lambda u, 
 * 
 *                with initial condition u(0) = 1, design parameter lambda < 0 
 *                and objective function:
 * 
 *                   J = 1/2 * ( 1/T \int_0^T || u ||^2 dt - Target )^2  
 *                      + gamma/2 || lamba ||^2
 * 
 *                for a precomputed 'Target' value. 
 *                Simple steepest-descent optimization is applied using a fixed 
 *                step size for design updates.   
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"

#include "braid_test.h"

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure contains information needed for time-stepping and evaluating the objective funcion */
typedef struct _braid_App_struct
{
   double    design;          /* Store the design variables in the app */
   double    gradient;        /* Store the gradient in the app - should be of same size as design! */
   double    target;          /* Target for tracking type objective function (inverse design) */
   double    gamma;           /* Relaxation parameter for the tracking type objective function */

   int       rank;            /* Rank of the processor */
   double    ntime;           /* Total number of time-steps */
} my_App;

/* Vector structure holds the state variable at one time-step */
typedef struct _braid_Vector_struct
{
   double value;
} my_Vector;

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

   /* Get the design variable from the app */
   double lambda = app->design;

   /* Use backward Euler to propagate solution */
   (u->value) = 1./(1. - lambda* (tstop-tstart))*(u->value);
   
   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   if (t == 0.0) /* Initial condition */
   {
      (u->value) = 1.0;
   }
   else /* All other time points set to arbitrary value */
   {
      (u->value) = 0.456;
   }
   *u_ptr = u;

   return 0;
}

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;

   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->value) = (u->value);
   *v_ptr = v;

   return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u);
   return 0;
}

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
   (y->value) = alpha*(x->value) + beta*(y->value);
   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double dot;

   dot = (u->value)*(u->value);
   *norm_ptr = sqrt(dot);
   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int        index;
   char       filename[255];
   FILE      *file;
   
   braid_AccessStatusGetTIndex(astatus, &index);
   sprintf(filename, "%s.%04d.%03d", "ex-01.out", index, app->rank);
   file = fopen(filename, "w");
   fprintf(file, "%.14e\n", (u->value));
   fflush(file);
   fclose(file);

   return 0;
}

int
my_BufSize(braid_App          app,
           int                *size_ptr,
           braid_BufferStatus bstatus)
{
   *size_ptr = sizeof(double);
   return 0;
}

int
my_BufPack(braid_App          app,
           braid_Vector       u,
           void               *buffer,
           braid_BufferStatus bstatus)
{
   double *dbuffer = buffer;

   dbuffer[0] = (u->value);
   braid_BufferStatusSetSize( bstatus, sizeof(double) );

   return 0;
}

int
my_BufUnpack(braid_App          app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus bstatus)
{
   double    *dbuffer = buffer;
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   (u->value) = dbuffer[0];
   *u_ptr = u;

   return 0;
}


/* Evaluate one term of the time-dependent discretized objective function for
 * vector u.  The result over all time values will be summed by Braid. */
int 
my_ObjectiveT(braid_App              app,
              braid_Vector           u,
              braid_ObjectiveStatus  ostatus,
              double                *objectiveT_ptr)
{
   double objT;
   int    ntime;

   /* Get the total number of time-steps */
   braid_ObjectiveStatusGetNTPoints(ostatus, &ntime);

   /* Evaluate the local objective: 1/N * u(t)**2 */
   objT = 1. / ntime * (u->value) * (u->value);

   *objectiveT_ptr = objT;
   return 0;
}


/* Evaluate the time-independent part of the objective function */
int
my_PostprocessObjective(braid_App   app,
                        double      integral,
                        double     *postprocess_ptr
                        )
{
   double J;

   /* Tracking-type functional */
   J  = 1./2. * pow(integral - app->target,2);
   
   /* Regularization term */
   J += (app->gamma) / 2. * pow(app->design,2);

   *postprocess_ptr = J;

   return 0;
}

int
my_PostprocessObjective_diff(braid_App   app,
                             double      integral,
                             double     *F_bar_ptr
                             )
{
   /* Derivative of tracking type function */
   *F_bar_ptr = integral - app->target;

   /* Derivative of regularization term */
   app->gradient = (app->gamma) * (app->design);

   return 0;
}


/* Transposed partial derivatives of objectiveT times F_bar */
int
my_ObjectiveT_diff(braid_App            app,
                  braid_Vector          u,
                  braid_Vector          u_bar,
                  braid_Real            F_bar,
                  braid_ObjectiveStatus ostatus)
{
   int    ntime;
   double ddu;      /* Derivative wrt u */
   double ddesign;  /* Derivative wrt design */

   /* Get the total number of time-steps */
   braid_ObjectiveStatusGetNTPoints(ostatus, &ntime);

   /* Partial derivative with respect to u times F_bar */
   ddu = 2. / ntime * u->value * F_bar;

   /* Partial derivative with respect to design times F_bar*/
   ddesign = 0.0 * F_bar;

   /* Update u_bar and add to gradient */
   u_bar->value  = ddu;
   app->gradient += ddesign;

   return 0;
}

/* Transposed partial derivatives of the step routine times u_bar */
int
my_Step_diff(braid_App              app,
                braid_Vector        ustop,
                braid_Vector        u,
                braid_Vector        ustop_bar,
                braid_Vector        u_bar,
                braid_StepStatus    status)
{

   /* Get the time step size */
   double tstop, tstart, deltat;
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltat = tstop - tstart;

   double ddu;      /* Derivative wrt u */
   double ddesign;  /* Derivative wrt design */

   /* Get the design from the app */
   double lambda = app->design;

   /* Transposed derivative of step wrt u times u_bar */
   ddu = 1./(1. - lambda * deltat) * (u_bar->value);
 
   /* Transposed derivative of step wrt design times u_bar */
   ddesign = (deltat * (u->value)) / pow(1. - deltat*lambda,2) * (u_bar->value);

   /* Update u_bar and add to gradient */
   u_bar->value      = ddu;              // Make sure to do "=" here, not "+="! 
   app->gradient    += ddesign;


   return 0;
}


/* Set the gradient to zero */
int 
my_ResetGradient(braid_App app)
{
   app->gradient = 0.0;

   return 0;
}


/* Function to allow for the computation of the gradient */
int
gradient_allreduce(braid_App app)
{
   double mygradient = app->gradient;
   double gradient; 

   /* Collect sensitivities from all processors and broadcast it */
   MPI_Allreduce(&mygradient, &gradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   app->gradient = gradient;

   return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;
   double        tstart, tstop;
   int           ntime, rank;
   double        design; 
   double        objective, gamma, target;
   double        gnorm, gtol;
   double        stepsize;
   int           iter, maxiter;
   double        mygradient;

   /* Define time domain: ntime intervals */
   ntime  = 50;
   tstart = 0.0;                             
   tstop  = tstart + ntime/2.;

   /* Set optimization variables */
   design           = -1.0;                  /* Initial design  */
   target           = 1.15231184218078e-01;  /* Precomputed target for inverse design optimization */
   gamma            = 0.0005;                /* Relaxation parameter for the objective function */
   stepsize         = 6.0;                   /* Step size for design updates */
   maxiter          = 100;                   /* Maximum number of optimization iterations */
   gtol             = 1e-6;                  /* Stopping criterion on the gradient norm */


   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->rank        = rank;
   app->design      = design;
   app->gradient    = 0.0;
   app->ntime       = ntime;
   app->target      = target;
   app->gamma       = gamma;

   /* Initialize XBraid */
   braid_Init(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

  /* Initialize adjoint-based gradient computation */
   braid_InitAdjoint( my_ObjectiveT, my_ObjectiveT_diff, my_Step_diff, my_ResetGradient, &core);

   /* Optional: Set the tracking type objective function and derivative */
   braid_SetPostprocessObjective(core, my_PostprocessObjective);
   braid_SetPostprocessObjective_diff(core, my_PostprocessObjective_diff);
  
   /* Set some typical Braid parameters */
   braid_SetMaxLevels(core, 2);            /* Number of time-grid levels */
   braid_SetCFactor(core, -1, 2);          /* Coarsening factor on all levels */
   braid_SetMaxIter(core, 10);             /* Maximum number of state and adjoint iterations */
   braid_SetAccessLevel(core, 1);          /* Access level: only after drive() finishes */
   braid_SetPrintLevel( core, 1);          /* Print level: output level for convergence of state and adjoint */
   braid_SetAbsTol(core, 1e-6);            /* Tolerance on state residual norm */
   braid_SetAbsTolAdjoint(core, 1e-6);     /* Tolerance on adjoint residual norm */


   /* Optimization iteration */
   for (iter = 0; iter < maxiter; iter++)
   {
      /* Run adjoint XBraid to compute objective function and gradient */
      braid_Drive(core);

      /* Get the objective function value */
      braid_GetObjective(core, &objective);

      /* Collect sensitivities from all processors and broadcast the sum */
      mygradient = app->gradient;
      MPI_Allreduce(&mygradient, &app->gradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      /* Compute norm of the gradient */
      gnorm = sqrt((app->gradient)*(app->gradient));

      /* Output */
      if (rank == 0) 
      {
         printf("\n %3d: Objective = %1.8e,  || Gradient || = %1.8e\n", iter, objective, gnorm);
      }

      /* Check optimization convergence */
      if (gnorm < gtol)
      {
         break;
      }

      /* Design update using simple steepest descent method with fixed stepsize */
      app->design -= stepsize * app->gradient;
   }


   /* Print some statistics about the optimization run */
   if (rank == 0)
   {
      if (iter == maxiter)
      {
         printf("\n Max. number of iterations reached! \n\n"); 
      }
      else
      {
         printf("\n");
         printf("  Optimization has converged.\n");
         printf("\n"); 
         printf("  Objective function value = %1.8e\n", objective);
         printf("  Gradient norm            = %1.8e\n", gnorm);
         printf("\n");
         printf("  optimization iterations  = %d\n", iter);
         printf("  max optim iterations     = %d\n", maxiter);
         printf("  gradient norm tolerance  = %1.1e", gtol);
         printf("\n");
      }
   }

   /* Print XBraid statistics */
   braid_PrintStats(core);

   /* Clean up */
   free(app);
   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}
