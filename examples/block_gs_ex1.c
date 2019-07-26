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
 * Example:       advec-diff-omgrit.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-04-adjoint
 *
 * Description:  Solves a simple optimal control problem in time-parallel:
 * 
 *                 min   0.5\int_0^T \int_0^1 (u(x,t)-u0(x))^2+alpha v(x,t)^2 dxdt
 * 
 *                  s.t.  du/dt + du/dx - nu d^2u/dx^2 = v(x,t)
 *                        u(0,t)=u(1,t)=0
 *                                  u(x,0)=u0(x)
 *
 *               Implements a steepest-descent optimization iteration
 *               using fixed step size for design updates.   
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "braid.h"
#include "braid_test.h"
#define PI 3.14159265
#define g(dt,dx) dt/(2*dx)
#define b(dt,dx,nu) nu*dt/(dx*dx)
/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

typedef struct _braid_App_struct
{
   int      myid;        /* Rank of the processor */
   double   gamma;       /* Relaxation parameter for objective function, v(x,t) */
   int      ntime;       /* Total number of time-steps (starting at time 0) */
   int      mspace;      /* Number of space points included in our state vector */
                         /* So including boundaries we have M+2 space points */

   double **w;           /* Adjoint vectors at each time point on my proc */
   double *U0;

} my_App;


/* Define the state vector at one time-step */
typedef struct _braid_Vector_struct
{
   double *values;     /* Holds the R^M state vector (u_1, u_2,...,u_M) */

} my_Vector;

/*--------------------------------------------------------------------------
 * Vector utility routines
 *--------------------------------------------------------------------------*/

void
vec_create(int size, double **vec_ptr)
{
   *vec_ptr = (double*) malloc( size*sizeof(double) );
}

void
vec_destroy(double *vec)
{
   free(vec);
}

/*------------------------------------*/

void
vec_copy(int size, double *invec, double *outvec)
{
   int i;
   for (i = 0; i < size; i++)
   {
      outvec[i] = invec[i];
   }
}

/*------------------------------------*/

void
vec_axpy(int size, double alpha, double *x, double *y)
{
   int i;
   for (i = 0; i < size; i++)
   {
      y[i] = y[i] + alpha*x[i];
   }
}

/*------------------------------------*/

void
vec_scale(int size, double alpha, double *x)
{
   int i;
   for (i = 0; i < size; i++)
   {
      x[i] = alpha*x[i];
   }
}

/*--------------------------------------------------------------------------
 * KKT component routines
 *--------------------------------------------------------------------------*/

/* This is the K=[A B C] matrix. It acts on a vector in R^M */
/* This function requies that M>=3, but this can easily be fixed later */
void
apply_Phi_Inv(double dt, double *u)
{
   u[0]=(1+dt)*u[0]+dt*u[1];
   vec_scale(2, (double)1/(1+dt), u);
}

/*------------------------------------*/

void
apply_PhiAdjoint_Inv(double dt, double *u)
{
   double u0_tmp=u[0];
   u[0]=(1+dt)*u[0];
   u[1]=dt*u0_tmp + u[1];
   vec_scale(2, (double)1/(1+dt), u);
}

/*------------------------------------*/

void
apply_Phi(double dt, double *u)
{
   u[0]=u[0] - dt*u[1];
   u[1]=(1+dt)*u[1];
}

/*------------------------------------*/

void
apply_PhiAdjoint(double dt, double *u)
{
   u[1]=-dt*u[0] + (1+dt)*u[1];
}

/*------------------------------------*/

void
apply_Uinv(double dt, double *u)
{
   u[0]=u[0]/(2*dt);
   u[1]=u[1]/(2*dt);
}

/*------------------------------------*/

void
apply_Vinv(double dt, double gamma, double *u)
{
   u[0]=u[0]/2*gamma*dt;
   u[1]=u[1]/2*gamma*dt;

}

/*------------------------------------*/

void
apply_DAdjoint(double dt, double *u)
{
    u[0]=0.0;
    u[1]=dt*u[1];
}

/*------------------------------------*/

/*--------------------------------------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/  

/*------------------------------------*/

/* This is only called from level 0 */

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;
   int mspace = (app->mspace);

   /* Allocate the vector */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   vec_create(mspace, &(u->values));

   for (int i = 0; i <= mspace-1; i++)
   {
      u->values[i] = ((double)braid_Rand())/braid_RAND_MAX;
   }

   *u_ptr = u;

   return 0;
}

/*------------------------------------*/

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   int mspace = (app->mspace);
   my_Vector *v;

   /* Allocate the vector */
   v = (my_Vector *) malloc(sizeof(my_Vector));
   vec_create(mspace, &(v->values));

   /* Clone the values */
   for (int i = 0; i<= mspace-1; i++)
   {
      v->values[i] = u->values[i];
   }

   *v_ptr = v;

   return 0;
}

/*------------------------------------*/

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u->values);
   free(u);

   return 0;
}

/*------------------------------------*/

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
   int mspace = (app->mspace);
   for (int i = 0; i <= mspace-1; i++)
   {
      (y->values)[i] = alpha*(x->values)[i] + beta*(y->values)[i];
   }

   return 0;
}

/*------------------------------------*/

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   int i;
   double dot = 0.0;
   int mspace = (app->mspace);
   for (i = 0; i <= mspace-1; i++)
   {
      dot += (u->values)[i]*(u->values)[i];
   }
   *norm_ptr = sqrt(dot);

   return 0;
}

/*------------------------------------*/

// ZTODO: Need to compute u from adjoint and it reqires communication

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int   done, index;
   int   mspace = (app->mspace);

   /* Print solution to file if simulation is over */
   braid_AccessStatusGetDone(astatus, &done);

   if (done)
   {
      /* Allocate w array in app (ZTODO: This only works on one proc right now) */
      if ((app->w) == NULL)
      {
         int  ntpoints;
         braid_AccessStatusGetNTPoints(astatus, &ntpoints);
         ntpoints++;  /* ntpoints is really the gupper index */
         (app->w) = (double **) calloc(ntpoints, sizeof(double *));
      }

      braid_AccessStatusGetTIndex(astatus, &index);
      if (app->w[index] != NULL)
      {
         free(app->w[index]);
      }
      vec_create(mspace, &(app->w[index]));
      vec_copy(mspace, (u->values), (app->w[index]));
   }
   return 0;
}

/*------------------------------------*/

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
   *size_ptr = 2*sizeof(double);
   return 0;
}

/*------------------------------------*/

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = buffer;
   int i;
   int mspace = (app->mspace); 

   for(i = 0; i < mspace; i++)
   {
      dbuffer[i] = (u->values)[i];
   }

   braid_BufferStatusSetSize( bstatus,  2*sizeof(double));

   return 0;
}

/*------------------------------------*/

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
{
   my_Vector *u = NULL;
   double    *dbuffer = buffer;
   int i;
   int mspace = (app->mspace); 

   /* Allocate memory */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->values = (double*) malloc( 2*sizeof(double) );

   /* Unpack the buffer */
   for(i = 0; i < mspace; i++)
   {
      (u->values)[i] = dbuffer[i];
   }

   *u_ptr = u;
   return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int
main(int argc, char *argv[])
{
   my_App     *app;
         
   double      tstart, tstop, dt, tol; 
   int         ntime, mspace, maxiter;
   double      gamma, start, end, time, seed;
   /*
   int         max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter;
   int         access_level, print_level;
   double      tol;
   */

   /* Define space domain. Space domain is between 0 and 1, mspace defines the number of steps */
   mspace = 2;
   ntime = 4096;

   /* Define some optimization parameters */
   gamma = .005;            /* parameter in the objective function */
   tol  = 1.0e-6;
   maxiter = 300;
   int arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         printf("\n");
         printf(" Solves the advection-diffusion model problem \n\n");
         printf("  min  1/2 \\int_0^T\\int_0^1 (u(x,t)-ubar(x))^2 + alpha*v(x,t)^2  dxdt \n\n");
         printf("  s.t.  u_t + u_x - nu*u_xx = v(x,t) \n");
         printf("        u(0,t) = u(1,t) = 0 \n\n");
         printf("        u(x,0) = u0(x) \n");
         printf("  -tstop <tstop>          : Upper integration limit for time\n");
         printf("  -ntime <ntime>          : Num points in time\n");
         printf("  -gamma <gamma>          : Constant Parameter in Objective Function  \n");
         printf("  -ml <max_levels>        : Max number of braid levels \n");
         printf("  -mi <maxiter>           : Max iterations \n");
         printf("  -tol <tol>              : Stopping tolerance \n");
         printf("  -seed <seed>            : Seed for initial guess \n");
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tstop") == 0 )
      {
         arg_index++;
         tstop = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-seed") == 0 )
      {
         arg_index++;
         seed = atoi(argv[arg_index++]);
      }        
      else if ( strcmp(argv[arg_index], "-gamma") == 0 )
      {
         arg_index++;
         gamma = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mi") == 0 )
      {
         arg_index++;
         maxiter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tol") == 0 )
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else
      {
         printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
         return (0);
      }
   }

   /* Define time domain and step */
   tstart = 0.0;             /* Beginning of time domain */
   tstop  = 1.0;             /* End of time domain*/
   dt = (tstop-tstart)/ntime; 
    

   /* Set up the app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->ntime    = ntime;
   app->mspace   = mspace;
   app->gamma    = gamma;
   app->w        = NULL;
   /* Set this to whatever u0 is. */
   double *U0 = (double*) malloc( mspace*sizeof(double) );
   double *u_init=(double*) malloc( mspace*sizeof(double) );
   U0[0]=0.0;
   U0[1]=-1.0; 
   
   /* Set the initial guess to be random */
   srand(seed);
   for(int i=0; i<mspace; i++)
   {   
      double random_value;
      random_value=(double)rand()/RAND_MAX*2.0-1.0; 
      u_init[i]=random_value;
   }

   app->U0       = U0;
   

   /* Start the Gauss-Seidel iterations */
   start=clock();
   double norm=0;
   double niters=0;
   double **w = (double **)malloc(ntime * sizeof(double*));
   for(int i = 0; i < ntime; i++) w[i] = (double *)malloc(mspace * sizeof(double));
   
   double *v = (double*) malloc(ntime * sizeof(double) );
   
   double **u = (double **)malloc(ntime * sizeof(double*));
   for(int i = 0; i < ntime; i++) u[i] = (double *)malloc(mspace * sizeof(double));
   
   double **res = (double **)malloc(ntime * sizeof(double*));
   for(int i = 0; i < ntime; i++) res[i] = (double *)malloc(mspace * sizeof(double));
   
   double **res1 = (double **)malloc(ntime * sizeof(double*));
   for(int i = 0; i < ntime; i++) res1[i] = (double *)malloc(mspace * sizeof(double));

   for(int i = 0; i < ntime; i++)
   {
      vec_copy(mspace, u_init, w[i]);
      vec_copy(mspace, u_init, u[i]);      
      vec_copy(mspace, u_init, res[i]);
      vec_copy(mspace, u_init, res1[i]);
      v[i]=u_init[0];
   }
   do
   {
      norm = 0;    
      /**************FORWARD SOLVE*******************/
      /* Solve Lu^(k+1)=g+Dv^(k) */
      u[0][0]=0;
      u[0][1]=dt*v[0];
      vec_axpy(mspace, 1.0, U0, u[0]);
      apply_Phi_Inv(dt, u[0]);
      for (int i = 1; i < ntime; i++)
      {
         u[i][0]=0;
         u[i][1]=dt*v[i];
         vec_axpy(mspace, 1.0, u[i-1], u[i]);
         apply_Phi_Inv(dt, u[i]);
      }

      /* Solve L*w^(k+1)=Uu^(k+1) */
      vec_copy(mspace, u[ntime-1], w[ntime-1]);
      vec_scale(mspace, -2*dt, w[ntime-1]); /* Apply U is the same as scaling by dxdt */
      apply_PhiAdjoint_Inv(dt, w[ntime-1]);
      for(int i = ntime-2; i >= 0; i--)
      {
         vec_copy(mspace, u[i], w[i]);
         vec_scale(mspace, -2*dt, w[i]);
         vec_axpy(mspace, 1.0, w[i+1], w[i]);
         apply_PhiAdjoint_Inv(dt, w[i]);
      }

      /* Solve Vv^(k+1)=Dw^(k+1) */
      for(int i = 0; i < ntime; i++)
      {
         v[i]=w[i][1]/(2*gamma); 
      }
      /****************************RESIDUAL*********************************/
      /* Compute residual for block Lu^(k+1)+Dv^(k+1)-g */
      vec_copy(mspace, u[0], res[0]);
      apply_Phi(dt, res[0]);
      
      res[0][1]=res[0][1]-dt*v[0];

      vec_axpy(mspace, -1.0, U0, res[0]);

      for (int i = 1; i < ntime; i++)
      {
         vec_copy(mspace, u[i], res[i]);
         apply_Phi(dt, res[i]);
         vec_axpy(mspace, -1.0, u[i-1], res[i]);
         
         res[i][1]=res[i][1]-dt*v[i];
      }
      for(int i = 0; i < ntime; i++)
      {
         for (int j = 0; j < mspace; j++)
         {
            norm = norm + res[i][j]*res[i][j];
         }
      }

      /* Compute residual for block Uu^(k+1)+L*w^(k+1)-k */
      for(int i = 0; i < ntime-1; i++)
      {
         vec_copy(mspace, u[i], res[i]);
         vec_scale(mspace, 2*dt, res[i]);

         vec_copy(mspace, w[i], res1[i]);
         apply_PhiAdjoint(dt, res1[i]);
         vec_axpy(mspace, 1.0, res1[i], res[i]);

         vec_copy(mspace, w[i+1], res1[i]);
         vec_axpy(mspace, -1.0, res1[i], res[i]);
      }
      vec_copy(mspace, u[ntime-1], res[ntime-1]);
      vec_scale(mspace, 2*dt, res[ntime-1]);

      vec_copy(mspace, w[ntime-1], res1[ntime-1]);
      apply_PhiAdjoint(dt, res1[ntime-1]);
      vec_axpy(mspace, 1.0, res1[ntime-1], res[ntime-1]);
      for(int i = 0; i < ntime; i++)
      {
         for (int j = 0; j < mspace; j++)
         {
            norm = norm + res[i][j]*res[i][j];
         }
      }

      /* Compute residual for block D*w^(k+1)+Vv^(k+1)-0 */
      for(int i = 0; i< ntime; i++)
      {
         res[i][0]=2*gamma*dt*v[i];
         res[i][0]=res[i][0]-dt*w[i][1];
         norm = norm + res[i][0]*res[i][0];
      }          
      norm = sqrt(norm);        
      /*****************************************/
      niters = niters + 1;
      printf("Residual: ");
      printf("%f", norm);
      printf("\n");
      printf("Iteration number: ");
      printf("%f", niters);
      printf("\n");
   }while(norm > tol && niters < maxiter);

   end=clock();
   time=(double)(end-start)/CLOCKS_PER_SEC;
   printf("The total run time is: ");
   printf("%f", time);
   printf(" seconds\n");
   /* Print out v, w, u and U0 */
  
   /**********************PRINT W OUT**********************/
   char  filename[255];
   FILE *file;
   int   i;
   sprintf(filename, "%s.%03d", "out/block_gs_ex1.out.w", 000);
   file = fopen(filename, "w");
   for (i = 0; i < (app->ntime); i++)
   {
      /*double **w = (app->w);*/
      fprintf(file, "%05d: % 1.14e, % 1.14e\n", (i+1), w[i][0], w[i][1]);
   }
   fflush(file);
   fclose(file);

   /**********************PRINT V OUT**********************/
   sprintf(filename, "%s.%03d", "out/block_gs_ex1.out.v", 000);
   file = fopen(filename, "w");
   for (i = 0; i < (app->ntime); i++)
   {
      fprintf(file, "%05d: % 1.14e\n", (i+1), v[0]);
   }
   fflush(file);
   fclose(file);

   /**********************PRINT U OUT**********************/
   sprintf(filename, "%s.%03d", "out/block_gs_ex1.out.u", 000);
   file = fopen(filename, "w");
   for (i = 0; i < (app->ntime); i++)
   {
      fprintf(file, "%05d: % 1.14e, % 1.14e\n", (i+1), u[i][0], u[i][1]);
   }
   fflush(file);
   fclose(file);

   /**********************PRINT TIME OUT**********************/
   {
      /*
      char    filename[255];
      FILE   *file;
      */
      sprintf(filename, "%s.%d", "out/block_gs_ex1.time", maxiter);
      file = fopen(filename, "w");
      fprintf(file, "%f", time);
      fflush(file);
      fclose(file);
   }   

   /**********************PRINT TOTAL RES OUT**********************/
   {
      /*
      char    filename[255];
      FILE   *file;
      */
      sprintf(filename, "%s.%d", "out/block_gs_ex1.res", maxiter);
      file = fopen(filename, "w");
      fprintf(file, "%f", norm);
      fflush(file);
      fclose(file);
   }      

      /**********************PRINT OUT IF CONVERGE OR DIVERGE**********************/
   {
      sprintf(filename, "%s.%d.%f", "out/block_gs_ex1.conv", ntime, gamma);
      file = fopen(filename, "w");
      if (isinf(norm)||isnan(norm))
      {
         fprintf(file, "%f", 0.0);
      }
      else
      {
         fprintf(file, "%f", 1.0);
      }
     
      fflush(file);
      fclose(file);
   }   

   return (0);
}