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
 * Example:       ex-04-serial.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-04-serial
 *
 * 
 * Description:  Solves a simple optimal control problem in time-serial:
 * 
 *                 min   \int_0^1 u_1(t)^2 + u_2(t)^2 + gamma c(t)^2  dt
 * 
 *                  s.t.  d/dt u_1(t) = u_2(t)
 *                        d/dt u_2(t) = -u_2(t) + c(t)
 * 
 *                 with initial condition u_1(0) = 0, u_2(0) = -1
 *                 and piecewise constant control c(t).  
 * 
 *               Implements a steepest-descent optimization iteration
 *               using fixed step size for design updates.   
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

void 
take_step(double* u,         /* state at current time */
          double *uleft,
          double  *design,    /* design at current time */
          double  deltaT,
          int mspace,
          double nu )   /* time step size */
{


   double dx = 1.0/(mspace-1);
   double dt = deltaT;

   int M= mspace;

   double A = ((dt*nu)/(dx*dx)) + (dt/(2*dx));
   double B = 1 - ((2*nu*dt)/(dx*dx));
   double C = (dt*nu)/(dx*dx) - (dt/(2*dx));
   double tmp_u_1 = u[0];
   double tmp_u_2 = u[1];
   double tmp_u_Mm1 = u[mspace-2];
   double tmp_u_M = u[mspace-1];
   
   
   for (int i = 1; i <= mspace - 2; i++)
   {
     u[i] = A*u[i-1] + B*u[i] + C*u[i+1] + deltaT*design[i] -uleft[i];
   }
   
   /* Deal with the u_1 and u_M vectors seperately */
   u[0] = B*tmp_u_1 + C*tmp_u_2 + deltaT*design[0] -uleft[i];
   u[mspace-1] = A*tmp_u_Mm1 + B*tmp_u_M + deltaT*design[mspace-1] -uleft[i];

}

/* Evaluate the time-dependent objective function at the current time */
double 
evalObjectiveT(double* u,        /* state at current time */
               double  *design,   /* design at current time */
               double  deltaT,   /* time step size */
               double  alpha,
               int mspace,
               double *u0)    /* relaxation parameter */
{
   double objectiveT;

   objectiveT  =  0;
   double dx = 1.0/(mspace-1);
   double dt = deltaT;
   for(int i=0; i<mspace; i++)
   {
      objectiveT += ((u[i] - u0[i])*(u[i] - u0[i]) + alpha*design[i]*design[i])*dx;
   }
   objectiveT *= dt*.5;

   return objectiveT;
}


/*--------------------------------------------------------------------------
 * Routines for solving the adjoint equation 
 *--------------------------------------------------------------------------*/

/* Compute the transposed partial derivative of take_step multiplied with current adjoint w */
double
*take_step_diff(double *w, 
               double  deltaT,
               int mspace,
               double nu)
{
   double *ddu = (double*)malloc(mspace*sizeof(double));
   double *gradientT = (double*)malloc(mspace*sizeof(double));
   double dx = 1.0/(mspace-1);
   double dt = deltaT;
   int M= mspace;

   double A = ((dt*nu)/(dx*dx)) + (dt/(2*dx));
    double B = 1 - ((2*nu*dt)/(dx*dx));
    double C = (dt*nu)/(dx*dx) - (dt/(2*dx));
    double tmp_w_1 = w[0];
    double tmp_w_2 = w[1];
    double tmp_w_Mm1 = w[M-2];
    double tmp_w_M = w[M-1];
    
    
    for (int i = 1; i <= M - 2; i++)
    {
       ddu[i] = C*w[i-1] + B*w[i] + A*w[i+1];
    }
    
    /* Deal with the u_1 and u_M vectors seperately */
    ddu[0] = B*tmp_w_1 + A*tmp_w_2;
    ddu[M-1] = C*tmp_w_Mm1 + B*tmp_w_M;


   /* derivative with respect to c  */
   for(int i=0; i<mspace; i++){
      gradientT[i] = deltaT * w[i];
   }

   /* Update adjoint */
   for(int i =0; i<mspace; i++)
   {
      w[i] = ddu[i];
   }
   

   free(ddu);

   return gradientT;
}            

/* Partial derivative of the objective */
double
*evalObjectiveT_diff(double *w, 
                    double *u, 
                    double  *design,
                    double  nu,
                    double  alpha,
                    double  deltaT,
                    int mspace,
                    double *u0)
{
   double *gradientT = (double*)malloc(mspace*sizeof(double));
   double *ddu = (double*)malloc(mspace*sizeof(double));
   double dx = 1.0/(mspace-1);
   double dt = deltaT;

   /* derivative with respect to u */
      for(int i=0; i<mspace; i++)
      {
         ddu[i] = (u[i]-u0[i]) * dx *dt;
      }

   /* Update the adjoint */
   for(int i=0; i<mspace; i++){
      w[i] = ddu[i];
   }

   /* derivative with respect to c */
   for(int i=0; i<mspace; i++)
   {
      gradientT[i] = design[i] * alpha * dx *dt;
   }

   return gradientT;
}

/* Advance an adjoint variable backwards in time 
 * and evaluate local gradient at that time step */
double
*take_adjoint_step(double *w, 
                    double *u, 
                    double  *design,
                    double  nu,
                    double  alpha,
                    double  deltaT,
                    int mspace,
                    double *u0)    /* time step size */
{
   double *gradientT = (double*)malloc(mspace*sizeof(double));

   /* transposed derivatives of evalObjectiveT */
   gradientT = evalObjectiveT_diff(w, u, design, nu, alpha, deltaT, mspace, u0);

   /* transposed derivatives of take_step times w */
   double *gradientT2 = take_step_diff(w, deltaT, mspace, nu);
   for(int i=0; i<mspace; i++)
   {
      gradientT[i] += gradientT2[i];
   }

   return gradientT;
}                 


/*--------------------------------------------------------------------------
 * Utility routines
 *--------------------------------------------------------------------------*/

/* Write state or adjoint vector to file */
void 
write_vec(char*   name,     /* Filename extension (ex-04.out.name) */
          double **u,     /* first components of the vector */
          int     ntime,
          int mspace )   /* total number of time steps (length of vector) */
{
   char  filename[255];
   FILE *file;
   int ts;

   /* Open file for output  */
   sprintf(filename, "advec-grad-serial.out.u");
   file = fopen(filename, "w");
   for(ts = 0; ts < ntime; ts++)
   {
      for(int j=0; j<mspace; j++){
         if(j!=mspace-1){
            fprintf(file, "% 1.14e, ", u[ts][j]);
         }
         else{
            fprintf(file, "% 1.14e\n", u[ts][j]);
         }
      }
   }
   fflush(file);
   fclose(file);
}              

/* Write design vector to file */
void 
write_design_vec(char*   name,     /* Filename extension (ex-04.out.name) */
                 double** v,      /* vector */
                 int     ntime,
                 int mspace )   /* total number of time steps (length of vector) */
{
   char  filename[255];
   FILE *file;
   int ts;

   sprintf(filename, "advec-grad-serial.out.v");
   file = fopen(filename, "w");
   for(ts = 0; ts < ntime; ts++)
   {
      for(int j=0; j<mspace; j++){
         if(j!=mspace-1){
            fprintf(file, "% 1.14e, ", v[ts][j]);
         }
         else{
            fprintf(file, "% 1.14e\n", v[ts][j]);
         }
      }
   }
   fflush(file);
   fclose(file);
}    

/* Compute the squared norm of a vector */
double 
compute_sqnorm(double **vector, 
               int     size, 
               int mspace)

{
   double norm = 0.0;
   int i;

   for(i = 0; i < size; i++) 
   {
      for(int j=0; j<mspace; j++)
      {
         norm += (vector[i][j]) * (vector[i][j]) ;
      }
   }

   return norm;
}



int main (int argc, char *argv[])
{
   double   tstart, tstop, deltaT, start, end, alpha, nu;
   int      ntime, ts, arg_index, mspace;
   int      maxiter,iter;
   double   objective;
   double   gamma;
   double  *u, *w; 
   double   gnorm;
   double   gtol;
   double   stepsize;
   double time;
   start=clock();
   /* Define time domain */
   ntime  = 64;                         /* Total number of time-steps */
   tstart = 0.0;                        /* Beginning of time-domain */
   tstop  = 1.0;                        /* End of time-domain */
   mspace = 16;

   /* Define optimization parameters */
   nu    = 0.01;                    /* Relaxation parameter in the objective function */
   alpha = 0.005;
   stepsize = 50.0;                     /* Step size for design updates */
   maxiter  = 500;                      /* Maximum number of optimization iterations */
   gtol     = 1e-6;                     /* Stopping criterion on the gradient norm */
   
   /* Parse command line */
   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         printf("\n");
         printf(" Solves a simple optimal control problem in time-serial on [0, 1] \n\n");
         printf("  min   \\int_0^1 u_1(t)^2 + u_2(t)^2 + gamma c(t)^2  dt \n\n");
         printf("  s.t.  d/dt u_1(t) = u_2(t) \n");
         printf("        d/dt u_2(t) = -u_2(t) + c(t) \n\n");
         printf("  -ntime <ntime>       : set num points in time\n");
         printf("  -gamma <gamma>       : Relaxation parameter in the objective function \n");
         printf("  -stepsize <stepsize> : Step size for design updates \n");
         printf("  -maxiter <maxiter>   : Maximum number of optimization iterations \n");
         printf("  -gtol <gtol>         : Stopping criterion on the gradient norm \n");
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
         deltaT = (tstop - tstart) / ntime;  /* recompute */
      }
      else if ( strcmp(argv[arg_index], "-alpha") == 0 )
      {
         arg_index++;
         alpha = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu") == 0 )
      {
         arg_index++;
         nu = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mspace") == 0 )
      {
         arg_index++;
         mspace = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-stepsize") == 0 )
      {
         arg_index++;
         stepsize = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-maxiter") == 0 )
      {
         arg_index++;
         maxiter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-gtol") == 0 )
      {
         arg_index++;
         gtol = atof(argv[arg_index++]);
      }
      else
      {
         printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
         return (0);
      }
   }


   deltaT = (tstop - tstart) / ntime;   /* Time-step size */

   /* Initialize the optimization variables */
   u        = (double*) malloc( mspace*sizeof(double) );      /* State at a time-step */
   w        = (double*) malloc( mspace*sizeof(double) );      /* Adjoint at a time-step */
   double *u0 = (double*) malloc(mspace*sizeof(double));
   for(int i=0;i<mspace;i++){
      if(i<mspace/2){
         u0[i]=1.;
      }
      else{
         u0[i]=0.;
      }
   }
   
   double **us = (double **)malloc(ntime * sizeof(double*));
   for(int i = 0; i < ntime; i++) us[i] = (double *)malloc(mspace * sizeof(double));
   
   double **design = (double **)malloc(ntime * sizeof(double*));
   for(int i = 0; i < ntime; i++) design[i] = (double *)malloc(mspace * sizeof(double));

   double **gradient = (double **)malloc(ntime * sizeof(double*));
   for(int i = 0; i < ntime; i++) gradient[i] = (double *)malloc(mspace * sizeof(double));

   for (ts = 0; ts < ntime; ts++)
   {
      for(int j=0; j<mspace; j++){
         us[ts][j]=0.;
         design[ts][j]=0.;
         gradient[ts][j] = 0.;
      }
   }

   /* Prepare optimization output */
   printf("\n#    Objective             || Gradient ||\n");

   /* Optimization iteration */
   for (iter = 0; iter < maxiter; iter++)
   {

      /* Set initial condition for the state */
      for(int i=0;i<mspace;i++){
         if(i<mspace/2){
            u[i]=1.;
         }
         else{
            u[i]=0;
         }
      }

      /* Main forward time-stepping loop */
      objective = 0.0;
      for (ts = 1; ts <= ntime; ts++)
      {
         /* Take the step */
         take_step(u, design[ts-1], deltaT, mspace, nu);

         /* Add to the objective function */
         objective += evalObjectiveT(u, design[ts-1], deltaT, alpha, mspace, u0);

         /* Store the state variable */
         for(int j =0; j<mspace; j++)
         {
            us[ts-1][j] = u[j];
         }
      }

      /* Set terminal condition for the adjoint */
      w[0] = 0.0;
      w[1] = 0.0;

      /* Main backward adjoint time-stepping loop */
      for (ts = ntime; ts >= 1; ts--)
      {
         /* Restore the state variable */
         for(int j =0; j<mspace; j++)
         {
             u[j] = us[ts-1][j];
         };

         /* Take adjoint step and evaluate gradient */
         gradient[ts-1] = take_adjoint_step(w, u, design[ts-1], nu, alpha, deltaT, mspace, u0);
      }

      /* Compute norm of gradient */
      gnorm = compute_sqnorm(gradient, ntime, mspace);
      gnorm = sqrt(gnorm);

      /* Output */
      printf("%3d  %1.14e  %1.14e\n", iter, objective, gnorm);

      /* Check optimization convergence */
      if (gnorm < gtol)
      {
         break;
      }

      /* Design update */
      for(ts = 0; ts < ntime; ts++) 
      {
         for(int j=0; j<mspace; j++)
         {
            design[ts][j] -= stepsize * gradient[ts][j];
         }
         
      }
   }

   if (iter == maxiter)
   {
      printf("\n Max. number of iterations reached.\n\n");

      char    filename[255];
      FILE   *file;
      sprintf(filename, "out/advec-diff-grad.conv.%d.%f", ntime, nu);
      file = fopen(filename, "w");
      fprintf(file, "%f", 0.0);
      fflush(file);
      fclose(file);
   }
   else
   {
      printf("\n Optimization has converged.\n\n");

      char    filename[255];
      FILE   *file;
      sprintf(filename, "out/advec-diff-grad.conv.%d.%f", ntime, nu);
      file = fopen(filename, "w");
      fprintf(file, "%f", 1.0);
      fflush(file);
      fclose(file);
   }

   end=clock();
   time = (double)(end-start)/CLOCKS_PER_SEC;
   printf("Total Run Time: %f s \n", time);
   {
      char    filename[255];
      FILE   *file;
      sprintf(filename, "%s.%d.%d", "out/advec-grad-serial.time", ntime,mspace);
      file = fopen(filename, "w");
      fprintf(file, "%f", time);
      fflush(file);
      fclose(file);
   }

   /* Write final state and design to file */
   write_vec("state", us, ntime, mspace);
   write_design_vec("design", design, ntime, mspace);

   /* Free memory */
   free(design);
   free(us);
   free(u);
   free(gradient);
   free(w);

   return (0);
}
