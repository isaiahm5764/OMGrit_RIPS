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
 *												u(x,0)=u0(x)
 *
 *               Implements a steepest-descent optimization iteration
 *               using fixed step size for design updates.   
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "braid.h"
#include "braid_test.h"

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

typedef struct _braid_App_struct
{
   int      myid;        /* Rank of the processor */
   double   alpha;       /* Relaxation parameter for objective function, v(x,t) */
   double   nu;          /* Diffusion coefficent, which we take to be large */
   int      ntime;       /* Total number of time-steps (starting at time 0) */
   int      mspace;      /* Number of space points included in our state vector 
                         /* So including boundaries we have M+2 space points */

   double **w;           /* Adjoint vectors at each time point on my proc */

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
apply_Phi(double dt, double dx, double nu, int M, double *u)
{	 
	 /* Define the A,B,C as in the stencil. These can maybe moved to the app struct */
	 /* Also copy the initial values of u_1, u_2, u_M-1, u_M */
	 double A = dt*nu/(dx*dx) + dt/(2*dx);
	 double B = 1 - 2*nu*dt/(dx*dx);
	 double C = dt*nu / (dx*dx) - dt/(2*dx);
	 double tmp_u_1 = u[0];
	 double tmp_u_2 = u[1];
	 double tmp_u_Mm1 = u[M-2];
	 double tmp_u_M = u[M-1];
	 
	 
	 for (int i = 1; i <= M - 2; i++)
	 {
		 u[i] = A*u[i-1] + B*u[i] + C*u[i+1];
	 }
	 
	 /* Deal with the u_1 and u_M vectors seperately */
	 u[0] = B*tmp_u_1 + C*tmp_u_2;
	 u[M-1] = A*tmp_u_Mm1 + B*tmp_u_M;
}

/*------------------------------------*/

void
apply_PhiAdjoint(double dt, double dx, double nu, int M, double *w)
{
	 /* Define the A,B,C as in the stencil. These can maybe moved to the app struct */
	 /* Also copy the initial values of u_1, u_2, u_M-1, u_M */
	 double A = dt*nu/(dx*dx) + dt/(2*dx);
	 double B = 1 - 2*nu*dt/(dx*dx);
	 double C = dt*nu / (dx*dx) - dt/(2*dx);
	 double tmp_w_1 = w[0];
	 double tmp_w_2 = w[1];
	 double tmp_w_Mm1 = w[M-2];
	 double tmp_w_M = w[M-1];
	 
	 
	 for (int i = 1; i <= M - 2; i++)
	 {
		 w[i] = C*w[i-1] + B*w[i] + A*w[i+1];
	 }
	 
	 /* Deal with the u_1 and u_M vectors seperately */
	 w[0] = B*tmp_w_1 + A*tmp_w_2;
	 w[M-1] = C*tmp_w_Mm1 + B*tmp_w_M;
}

/*------------------------------------*/

void
apply_Uinv(double dt, double dx, int M, double *u)
{
   for (int i = 0; i <= M-1; i++)
	 {
		 u[i]/=dx*dt;
	 }
}

/*------------------------------------*/

void
apply_Vinv(double dt, double dx, double alpha, int M, double *v)
{
	for (int i = 0; i <= M-1; i++)
	{
		v[i] /= alpha*dx*dt;
	}
   
}

/*------------------------------------*/

void
apply_D(double dt, int M, double *v)
{
	 for (int i = 0; i <= M-1; i++)
	 {
		 v[i] *= -dt;
	 }
}

/*------------------------------------*/

void
apply_DAdjoint(double dt, int M, double *v)
{
	 for (int i = 0; i <= M-1; i++)
	 {
		 v[i] /= -dt;
	 }
}

/*------------------------------------*/

/*--------------------------------------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/

/* Compute A(u) - f */

int
my_TriResidual(braid_App       app,
               braid_Vector    uleft,
               braid_Vector    uright,
               braid_Vector    f,
               braid_Vector    r,
               braid_Int       homogeneous,
               braid_TriStatus status,
               braid_Vector    u_0)
{
   double  t, tprev, tnext, dt, dx;
   double  nu = (app->nu);
	 double  alpha = (app->alpha);
   double *rtmp, *utmp;
   int     level, index;
   int     mspace = (app->mspace);
   
   braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
   braid_TriStatusGetLevel(status, &level);
   braid_TriStatusGetTIndex(status, &index);

   /* Get the time-step size */
   if (t < tnext)
   {
      dt = tnext - t;
   }
   else
   {
      dt = t - tprev;
   }

   /* Get the space-step size */
   dx = dx/(mspace+2);

   /* Create temporary vectors */
   vec_create(mspace, &rtmp);
   vec_create(mspace, &utmp);

   /* Compute action of center block */

   /* rtmp = U_i^{-1} u */
   vec_copy(mspace, (r->values), utmp);
   apply_Uinv(dt, dx, mspace, utmp);
   vec_copy(mspace, utmp, rtmp);

   /* rtmp = rtmp + D_i^T V_i^{-1} D_i^T u */
   vec_copy(mspace, (r->values), utmp);
   apply_DAdjoint(dt, mspace, utmp);
   apply_Vinv(dt, dx, alpha, mspace, utmp);
   apply_D(dt, mspace, utmp);
   vec_axpy(mspace, 1.0, utmp, rtmp);

   /* rtmp = rtmp + Phi_i U_{i-1}^{-1} Phi_i^T u */
   /* This term is zero at time 0, since Phi_0 = 0 */
   if (uleft != NULL)
   {
      vec_copy(mspace, (r->values), utmp);
      apply_PhiAdjoint(dt, dx, nu, mspace, utmp);
      apply_Uinv(dt, dx, mspace, utmp);
      apply_Phi(dt, dx, nu, mspace, utmp);
      vec_axpy(mspace, 1.0, utmp, rtmp);
   }

   /* Compute action of west block */
   if (uleft != NULL)
   {
      /* rtmp = rtmp - Phi_i U_{i-1}^{-1} uleft */
      vec_copy(mspace, (uleft->values), utmp);
      apply_Uinv(dt, dx, mspace, utmp);
      apply_Phi(dt, dx, nu, mspace, utmp);
      vec_axpy(mspace, -1.0, utmp, rtmp);
   }
   
   /* Compute action of east block */
   if (uright != NULL)
   {
      /* rtmp = rtmp - U_i^{-1} Phi_{i+1}^T uright */
      vec_copy(mspace, (uright->values), utmp);
      apply_PhiAdjoint(dt, dx, nu, mspace, utmp);
      apply_Uinv(dt, dx, mspace, utmp);
      vec_axpy(mspace, -1.0, utmp, rtmp);
   }

   /* Subtract rhs gbar (add g) in non-homogeneous case */
   /* This will also checks if*/
   if ((!homogeneous) && (index == 0))
   {
      /* rtmp = rtmp + g; g = Phi_0 U^0 */
		/* 
      utmp[0] =  0.0;
      utmp[mspace-1] = 0.0;
      apply_Phi(dt, dx, nu, mspace, utmp);
      vec_axpy(mspace, 1.0, utmp, rtmp);
      */

      /* COULD I DO THIS WITHOUT MAKING THE TMP VECTORS */
      for (int i = 0; i <= mspace-1; i++)
      {
         utmp[i]=u0[i];  
      }
      apply_Phi(dt, dx, nu, mspace, utmp);
      vec_axpy(mspace, 1.0, utmp, rtmp);
   }

   /* Subtract rhs f */
   if (f != NULL)
   {
      /* rtmp = rtmp - f */
      vec_axpy(mspace, -1.0, (f->values), rtmp);
   }
   
   /* Copy temporary residual vector into residual */
   vec_copy(mspace, rtmp, (r->values));
   
   /* Destroy temporary vectors */
   vec_destroy(rtmp);
   vec_destroy(utmp);
   
   return 0;
}   

/*------------------------------------*/

/* Solve A(u) = f */

int
my_TriSolve(braid_App       app,
            braid_Vector    uleft,
            braid_Vector    uright,
            braid_Vector    f,
            braid_Vector    u,
            braid_Int       homogeneous,
            braid_TriStatus status)
{
   double  t, tprev, tnext, dt;
   double *utmp, *rtmp;
   int mspace = (app->mspace);
   
   /* Get the time-step size */
   braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
   if (t < tnext)
   {
      dt = tnext - t;
   }
   else
   {
      dt = t - tprev;
   }

   /* Create temporary vector */
   vec_create(mspace, &utmp);

   /* Initialize temporary solution vector */
   vec_copy(mspace, (u->values), utmp);
   
   /* Compute residual */
   my_TriResidual(app, uleft, uright, f, u, homogeneous, status);

   /* Apply center block preconditioner (multiply by \tilde{C}^-1) to -r
    *
    * Using \tilde{C} = | 1/dt            0             |
    *                   |  0    ( 1/dt + dt/(2*gamma) ) |
    */

   /* JTODO: update lines below for general u vectors  */

   rtmp = (u->values);
   if (uleft != NULL)
   {
      rtmp[0] = -rtmp[0]*dt;
      rtmp[1] = -rtmp[1]/(1/dt + dt/(2*gamma));
   }
   else
   {
      /* At the leftmost point, use a different center coefficient approximation */
      rtmp[0] = -rtmp[0]*(2*dt);
      rtmp[1] = -rtmp[1]/(1/(2*dt) + dt/(2*gamma));
   }

   /* Complete residual update */
   vec_axpy(mspace, 1.0, utmp, (u->values));
   
   /* no refinement */
   braid_TriStatusSetRFactor(status, 1);

   return 0;
}   

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
   vec_create(mspace+2, &(u->values));

   u->values[0] = ((double)braid_Rand())/braid_RAND_MAX;
   u->values[1] = ((double)braid_Rand())/braid_RAND_MAX;

   *u_ptr = u;

   return 0;
}

/*------------------------------------*/

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;
   int mspace = (app->mspace); 
   /* Allocate the vector */
   v = (my_Vector *) malloc(sizeof(my_Vector));
   vec_create(mspace+2, &(v->values));

   /* Clone the values */
   v->values[0] = u->values[0];
   v->values[1] = u->values[1];

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

   (y->values)[0] = alpha*(x->values)[0] + beta*(y->values)[0];
   (y->values)[1] = alpha*(x->values)[1] + beta*(y->values)[1];

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

   for (i = 0; i < 2; i++)
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
      vec_create(2, &(app->w[index]));
      vec_copy(2, (u->values), (app->w[index]));
   }

//   {
//      char  filename[255];
//      FILE *file;
//      int  iter;
//      braid_AccessStatusGetIter(astatus, &iter);
//
//      braid_AccessStatusGetTIndex(astatus, &index);
//      sprintf(filename, "%s.%02d.%04d.%03d", "ex-04.out", iter, index, app->myid);
//      file = fopen(filename, "w");
//      fprintf(file, "%1.14e, %1.14e\n", (u->values)[0], (u->values)[1]);
//      fflush(file);
//      fclose(file);
//   }


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

   for(i = 0; i < 2; i++)
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

   /* Allocate memory */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->values = (double*) malloc( 2*sizeof(double) );

   /* Unpack the buffer */
   for(i = 0; i < 2; i++)
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
   braid_Core  core;
   my_App     *app;
         
   double      tstart, tstop, dt; 
   int         rank, ntime, mspace, arg_index;
   double      alpha, nu;
   int         max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter;
   int         access_level, print_level;
   double      tol;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   /* Define space domain. Space domain is between 0 and 1, mspace defines the number of steps */
   mspace = 20;

   /* Define time domain */
   ntime  = 20;              /* Total number of time-steps */
   tstart = 0.0;             /* Beginning of time domain */
   tstop  = 1.0;             /* End of time domain*/

   /* Define some optimization parameters */
   alpha = 0.005;            /* parameter in the objective function */
   nu    = 1;                /* parameter in PDE */

   /* Define some Braid parameters */
   max_levels     = 30;
   min_coarse     = 1;
   nrelax         = 1;
   nrelaxc        = 7;
   maxiter        = 20;
   cfactor        = 2;
   tol            = 1.0e-6;
   access_level   = 1;
   print_level    = 2;

   /* Parse command line */
   arg_index = 1;
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
         printf("  -mspace <mspace>        : Num points in space\n");
         printf("  -nu <nu>                : Constant Parameter in PDE  \n");
         printf("  -alpha <alpha>          : Constant Parameter in Objective Function  \n");
         printf("  -ml <max_levels>        : Max number of braid levels \n");
         printf("  -num  <nrelax>          : Num F-C relaxations\n");
         printf("  -nuc <nrelaxc>          : Num F-C relaxations on coarsest grid\n");
         printf("  -mi <maxiter>           : Max iterations \n");
         printf("  -cf <cfactor>           : Coarsening factor \n");
         printf("  -tol <tol>              : Stopping tolerance \n");
         printf("  -access <access_level>  : Braid access level \n");
         printf("  -print <print_level>    : Braid print level \n");
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
      else if ( strcmp(argv[arg_index], "-mspace") == 0 )
      {
         arg_index++;
         mspace = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ml") == 0 )
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu") == 0 )
      {
         arg_index++;
         nu = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-alpha") == 0 )
      {
         arg_index++;
         alpha = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-num") == 0 )
      {
         arg_index++;
         nrelax = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nuc") == 0 )
      {
         arg_index++;
         nrelaxc = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mi") == 0 )
      {
         arg_index++;
         maxiter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-cf") == 0 )
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tol") == 0 )
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-access") == 0 )
      {
         arg_index++;
         access_level = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-print") == 0 )
      {
         arg_index++;
         print_level = atoi(argv[arg_index++]);
      }
      else
      {
         printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
         return (0);
      }
   }

   /* Solve adjoint equations starting at time point t1=dt, recognizing that
    * braid will label this time point as index=0 instead of 1 */
   dt = (tstop-tstart)/ntime;

   /* Set up the app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->mspace   = mspace;
   app->nu       = nu;
   app->alpha    = alpha;
   app->w        = NULL;

   /* Initialize XBraid */
   braid_InitTriMGRIT(MPI_COMM_WORLD, MPI_COMM_WORLD, dt, tstop, ntime-1, app,
                      my_TriResidual, my_TriSolve, my_Init, my_Clone, my_Free,
                      my_Sum, my_SpatialNorm, my_Access,
                      my_BufSize, my_BufPack, my_BufUnpack, &core);

   /* Set some XBraid(_Adjoint) parameters */
   braid_SetMaxLevels(core, max_levels);
   braid_SetMinCoarse(core, min_coarse);
   braid_SetNRelax(core, -1, nrelax);
   if (max_levels > 1)
   {
      braid_SetNRelax(core, max_levels-1, nrelaxc); /* nrelax on coarsest level */
   }
   braid_SetCFactor(core, -1, cfactor);
   braid_SetAccessLevel(core, access_level);
   braid_SetPrintLevel( core, print_level);       
   braid_SetMaxIter(core, maxiter);
   braid_SetAbsTol(core, tol);

   /* Parallel-in-time TriMGRIT simulation */
   braid_Drive(core);

   if (access_level > 0)
   {
      /* Print adjoint w to file */
      {
         char  filename[255];
         FILE *file;
         int   i;

         sprintf(filename, "%s.%03d", "ex-04.out.w", (app->myid));
         file = fopen(filename, "w");
         for (i = 0; i < (app->ntime); i++)
         {
            double **w = (app->w);

            fprintf(file, "%05d: % 1.14e, % 1.14e\n", (i+1), w[i][0], w[i][1]);
         }
         fflush(file);
         fclose(file);
      }

      /* Compute state u from adjoint w and print to file */
      {
         char    filename[255];
         FILE   *file;
         int     i;
         double *u;

         sprintf(filename, "%s.%03d", "ex-04.out.u", (app->myid));
         file = fopen(filename, "w");
         vec_create(2, &u);
         for (i = 0; i < (app->ntime); i++)
         {
            double **w = (app->w);

            if ((i+1) < (app->ntime))
            {
               vec_copy(2, w[i+1], u);
               apply_PhiAdjoint(dt, u);
               vec_axpy(2, -1.0, w[i], u);
            }
            else
            {
               vec_copy(2, w[i], u);
               vec_scale(2, -1.0, u);
            }
            apply_Uinv(dt, u);

            fprintf(file, "%05d: % 1.14e, % 1.14e\n", (i+1), u[0], u[1]);
         }
         vec_destroy(u);
         fflush(file);
         fclose(file);
      }

      /* Compute control v from adjoint w and print to file */
      {
         char    filename[255];
         FILE   *file;
         int     i;
         double *v;

         sprintf(filename, "%s.%03d", "ex-04.out.v", (app->myid));
         file = fopen(filename, "w");
         vec_create(2, &v);
         for (i = 0; i < (app->ntime); i++)
         {
            double **w = (app->w);

            apply_DAdjoint(dt, w[i], v);
            vec_scale(1, -1.0, v);
            apply_Vinv(dt, (app->gamma), v);

            fprintf(file, "%05d: % 1.14e\n", (i+1), v[0]);
         }
         vec_destroy(v);
         fflush(file);
         fclose(file);
      }
   }

   free(app);
   
   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}