/*BHEADER**********************************************************************
 * Written by Isaiah Meyers, Joseph Munar, Eric Neville, Tom Overman
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
 * Example:       advec-diff-implicit.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make advec-diff-implicit
 *
 * Description:  Solves a linear optimal control problem in time-parallel:
 * 
 *                 min   0.5\int_0^T \int_0^1 (u(x,t)-u0(x))^2+alpha v(x,t)^2 dxdt
 * 
 *                  s.t.  du/dt + du/dx - nu d^2u/dx^2 = v(x,t)
 *                        u(0,t)=u(1,t)=0
 *                                  u(x,0)=u0(x)
 *
 *                 
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
   double   alpha;       /* Relaxation parameter for objective function, v(x,t) */
   double   nu;          /* Diffusion coefficent, which we take to be large */
   int      ntime;       /* Total number of time-steps (starting at time 0) */
   int      mspace;      /* Number of space points included in our state vector */
                         /* So including boundaries we have M+2 space points */

   int      ilower;      /* Lower index for my proc */
   int      iupper;      /* Upper index for my proc */
   int      npoints;     /* Number of time points on my proc */

   double **w;           /* Adjoint vectors at each time point on my proc */
   double *U0;
   double *ai;
   double *li;


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

/* This is the application of A inverse*/

void
apply_Phi(double dt, double dx, double nu, int M, double *u, double *l, double *a)
{   
   /* First solve Lw=u (Lw=f) */
   double *w;
   vec_create(M, &w);
   double *f;
   vec_create(M, &f);
   vec_copy(M, u, f);
   w[0]=f[0];
   for (int i = 1; i < M; i++)
   {
      w[i]=f[i]-l[i-1]*w[i-1];
   }

   /* Now solve Ux=w */ 
   double b = g(dt,dx)-b(dt, dx, nu);
   u[M-1]=w[M-1]/a[M-1];
   for (int i = M-2; i >= 0; i--)
   {
      u[i]=(w[i]-b*u[i+1])/a[i];      
   }
}

/*This is the application of A inverse transpose*/

void
apply_PhiAdjoint(double dt, double dx, double nu, int M, double *u, double *l, double *a)
{
   /* First solve U^Tw=u (U^Tw=f) */
   double *w;

   vec_create(M, &w);
   double *f;
   vec_create(M, &f);
   vec_copy(M, u, f);
   double b = g(dt,dx)-b(dt, dx, nu);
   w[0]=f[0]/a[0];
   for (int i = 1; i < M; i++)
   {
      w[i]=(f[i]-w[i-1]*b)/a[i];
   }

   /* Now solve L^Tx=w */ 
   u[M-1]=w[M-1];
   for (int i = M-2; i >= 0; i--)
   {
      u[i]=w[i]-l[i]*u[i+1];      
   }
}

/*------------------------------------*/

void
apply_A(double dt, double dx, double nu, int M, double *u)
{
   double A = -g(dt,dx)-b(dt,dx,nu);
   double B = 1+2*b(dt,dx,nu);
   double C = g(dt,dx)-b(dt,dx,nu);
   double *uold;
   vec_create(M, &uold);
   vec_copy(M, u, uold);
   u[0]=B*uold[0]+C*uold[1];
   u[M-1]=A*uold[M-2]+B*uold[M-1];
   for(int i = 1; i <= M-2; i++)
   {
      u[i]=A*uold[i-1]+B*uold[i]+C*uold[i+1];
   }
}

/*------------------------------------*/

void
apply_Aadjoint(double dt, double dx, double nu, int M, double *u)
{
   double A = -g(dt,dx)-b(dt,dx,nu);
   double B = 1+2*b(dt,dx,nu);
   double C = g(dt,dx)-b(dt,dx,nu);
   double *uold;
   vec_create(M, &uold);
   vec_copy(M, u, uold);
   u[0]=B*uold[0]+A*uold[1];
   u[M-1]=C*uold[M-2]+B*uold[M-1];
   for(int i = 1; i <= M-2; i++)
   {
      u[i]=C*uold[i-1]+B*uold[i]+A*uold[i+1];
   }
}

/*------------------------------------*/

void
apply_Uinv(double dt, double dx, int M, double *u)
{
   for (int i = 0; i <= M-1; i++)
   {
      u[i] /= dx*dt;
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
apply_D(double dt, double dx, double nu, int M, double *v, double *l, double *a)
{
   for (int i = 0; i <= M-1; i++)
   {
      v[i] *= dt;
   }
}

/*------------------------------------*/

void
apply_DAdjoint(double dt, double dx, double nu, int M, double *v, double *l, double *a)
{
   for (int i = 0; i <= M-1; i++)
   {
      v[i] *= dt;
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
               braid_TriStatus status)
{
   double  t, tprev, tnext, dt, dx;
   double  nu = (app->nu);
   double  alpha = (app->alpha);
   double *rtmp, *utmp;
   int     level, index;
   int     mspace = (app->mspace);
   double *u0 = (app->U0);
   double *li = (app->li);
   double *ai = (app->ai);

   
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
   dx = 1/((double)(mspace+1));

   /* Create temporary vectors */
   vec_create(mspace, &rtmp);
   vec_create(mspace, &utmp);


   /* Compute action of center block */

   /* rtmp = U_i^{-1}AA^T u */
   vec_copy(mspace, (r->values), utmp);
   apply_A(dt, dx, nu, mspace, utmp);
   apply_Aadjoint(dt, dx, nu, mspace, utmp);
   apply_Uinv(dt, dx, mspace, utmp);
   vec_copy(mspace, utmp, rtmp);

   /* rtmp = rtmp + D_i^T V_i^{-1} D_i^T u */
   vec_copy(mspace, (r->values), utmp);
   apply_DAdjoint(dt, dx, nu, mspace, utmp, li, ai);
   apply_Vinv(dt, dx, alpha, mspace, utmp);
   apply_D(dt, dx, nu, mspace, utmp, li, ai);
   vec_axpy(mspace, 1.0, utmp, rtmp);

   /* rtmp = rtmp + Phi_i U_{i-1}^{-1} Phi_i^T u */
   /* This term is zero at time 0, since Phi_0 = 0 */
   if (uleft != NULL)
   {
      vec_copy(mspace, (r->values), utmp);
      apply_Uinv(dt, dx, mspace, utmp);
      vec_axpy(mspace, 1.0, utmp, rtmp);
   }


   /* Compute action of west block */
   if (uleft != NULL)
   {
      vec_copy(mspace, (uleft->values), utmp);
      apply_Aadjoint(dt, dx, nu, mspace, utmp);
      apply_Uinv(dt, dx, mspace, utmp);
      vec_axpy(mspace, -1.0, utmp, rtmp);
   }


   /* Compute action of east block */
   if (uright != NULL)
   {
      vec_copy(mspace, (uright->values), utmp);
      apply_A(dt, dx, nu, mspace, utmp);
      apply_Uinv(dt, dx, mspace, utmp);
      vec_axpy(mspace, -1.0, utmp, rtmp);
   }


   /* No change for index 0 */
   if (!homogeneous)
   {
      vec_copy(mspace, u0, utmp);
      vec_axpy(mspace, 1.0, utmp, rtmp);

      vec_copy(mspace, u0,utmp);
      apply_A(dt, dx, nu, mspace, utmp);
      vec_axpy(mspace, -1.0, utmp, rtmp); 
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
   double  t, tprev, tnext, dt, dx;
   double *utmp, *rtmp;
   int mspace = (app->mspace);
   double nu = (app->nu);
   double *li = (app->li);
   double *ai = (app->ai);
   
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

   /* Get the space-step size */
   dx = 1/((double)(mspace+1));

   /* Create temporary vector */
   vec_create(mspace, &utmp);

   /* Initialize temporary solution vector */
   vec_copy(mspace, (u->values), utmp);
   
   /* Compute residual */
   my_TriResidual(app, uleft, uright, f, u, homogeneous, status);

   /* Apply center block preconditioner (multiply by \tilde{C}^-1) to -r
    *
    * Using [\tilde{C_i}] = 2AA^T
    * 
   */

   rtmp = (u->values);
   vec_scale(mspace, -1.0*dx*dt, rtmp);
   apply_Phi(dt, dx, nu, mspace, rtmp, li, ai);
   apply_PhiAdjoint(dt, dx, nu, mspace, rtmp, li, ai);
   vec_scale(mspace, .5, rtmp);


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
   for (int i = 0; i <= mspace-1; i++)
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
   int   done, index, ii;
   int   mspace = (app->mspace);

   /* Print solution to file if simulation is over */
   braid_AccessStatusGetDone(astatus, &done);
   if (done)
   {
      braid_AccessStatusGetILowerUpper(astatus, &(app->ilower), &(app->iupper));
      (app->npoints) = (app->iupper) - (app->ilower) + 1;

      /* Allocate w array in app */
      if ((app->w) == NULL)
      {
         (app->w) = (double **) calloc((app->npoints), sizeof(double *));
      }

      braid_AccessStatusGetTIndex(astatus, &index);
      ii = index - (app->ilower);
      if (app->w[ii] != NULL)
      {
         free(app->w[ii]);
      }
      vec_create(mspace, &(app->w[ii]));
      vec_copy(mspace, (u->values), (app->w[ii]));
   }


   //  Below prints U, V, and W after selected iterations. This can then be plotted to show how the space-time solution changes after iterations. 

   /*   char  filename[255];
      FILE *file;
      int  iter;
      braid_AccessStatusGetIter(astatus, &iter);
      braid_AccessStatusGetTIndex(astatus, &index);
      / * file format is advec-diff-btcs.out.{iteration #}.{time index} * /
      if(iter%1==0){
         sprintf(filename, "%s.%04d.%04d", "out/advec-diff-btcs.v.out", iter, index);
         file = fopen(filename, "w");
         for(int i = 0; i<mspace; i++){
             if(i<mspace-1){
                fprintf(file, "%1.14e, ", (u->values)[i]);
             }
             else{
                fprintf(file, "%1.14e", (u->values)[i]);
             }
         }
      fflush(file);
      fclose(file);
      }
   */
   return 0;
}

/*------------------------------------*/

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
  int mspace = (app->mspace); 
   *size_ptr = mspace*sizeof(double);
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

   braid_BufferStatusSetSize( bstatus,  mspace*sizeof(double));

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
   u->values = (double*) malloc( mspace*sizeof(double) );

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

   braid_Core  core;
   my_App     *app;
         
   double      tstart, tstop, dt, dx, start, end; 
   int         rank, ntime, mspace, arg_index;
   double      alpha, nu;
   int         max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter;
   int         access_level, print_level;
   double      tol;
   double time;
   

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   /* Define space domain. Space domain is between 0 and 1, mspace defines the number of steps */
   mspace = 8;
   ntime  = 256;

   /* Define some optimization parameters */
   alpha = .005;            /* parameter in the objective function */
   nu    = .3;               /* parameter in PDE */

   /* Define some Braid parameters */
   max_levels     = 3;
   min_coarse     = 1;
   nrelax         = 1;
   nrelaxc        = 30;
   maxiter        = 300;
   cfactor        = 2;
   tol            = 1.0e-6;
   access_level   = 2;
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
         nu = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-alpha") == 0 )
      {
         arg_index++;
         alpha = atof(argv[arg_index++]);
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

   /* Define the space step */
   dx = (double)1/(mspace+1);

   /* Define time domain and step */
   tstart = 0.0;             /* Beginning of time domain */
   tstop  = 1.0;             /* End of time domain*/
   dt     = (tstop-tstart)/ntime; 
    

   /* Set up the app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->mspace   = mspace;
   app->nu       = nu;
   app->alpha    = alpha;
   app->w        = NULL;

   /* Set this to u0 in the problem formulation */
   double *U0 = (double*) malloc( ntime*sizeof(double) );
   for(int i=0; i<mspace/2; i++)
   {
      U0[i] = 1;
   }

   for(int i=mspace/2; i<mspace; i++)
   {
      U0[i] = 0;
   }

   app->U0       = U0;

   /* Find elements of LU decomposition of A */
   double *ai = (double*) malloc( mspace*sizeof(double) );
   double *li = (double*) malloc( (mspace-1)*sizeof(double) );
   ai[0] = 1+2*b(dt,dx,nu);
   for(int i=1; i<mspace; i++)
   {
      li[i-1] = -(b(dt,dx,nu)+g(dt,dx))/ai[i-1];
      ai[i] = ai[0]+(b(dt,dx,nu)-g(dt,dx))*li[i-1];
   }
   app->ai       = ai;
   app->li       = li;

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
   start=clock();
   braid_Drive(core);
   end=clock();

   dx = 1/((double)(mspace+1));;
   
    /* Writes solutions to files */  
    if (access_level > 0)
    {
        /* Print adjoint w to file */
        char  filename[255];
        FILE *file;
        int   i,j, index;

        sprintf(filename, "%s.%03d", "advec-diff-implicit.out.w", (app->myid));
        file = fopen(filename, "w");
        for (i = 0; i < (app->npoints); i++)
        {
            double **w = (app->w);
            index = (app->ilower) + i +1;
            fprintf(file, "%05d: ", index);
            for(j=0; j <mspace; j++)
            {
                if(j==mspace-1)
                {
                    fprintf(file, "% 1.14e", w[i][j]);
                }
                else
                {
                    fprintf(file, "% 1.14e, ", w[i][j]);
                }
            }
            fprintf(file, "\n");
        }
        fflush(file);
        fclose(file);

        double *v;

        double **vs = (double **)malloc(app->npoints * sizeof(double*));
        for(int i = 0; i < app->npoints; i++) vs[i] = (double *)malloc(mspace * sizeof(double));

        /* Compute control v from adjoint w and print to file */
        sprintf(filename, "%s.%03d", "advec-diff-implicit.out.v", (app->myid));
        file = fopen(filename, "w");
        vec_create((app->mspace), &v);
        for (i = 0; i < (app->npoints); i++)
        {
            double **w = (app->w);
            vec_copy(mspace, w[i], v);
            apply_DAdjoint(dt, dx, nu, mspace, v, li, ai);
            apply_Vinv(dt, dx, alpha, mspace,v);

            fprintf(file, "%05d: ", ((app->ilower)+i+1));
            for (j = 0; j < (app->mspace); j++)
            {
                vs[i][j] = v[j];
                if(j==mspace-1)
                {
                    fprintf(file, "% 1.14e", v[j]);
                }
                else
                {
                    fprintf(file, "% 1.14e, ", v[j]);
                }
            }
            fprintf(file, "\n");

        }
        vec_destroy(v);
        fflush(file);
        fclose(file);


        char filename1[255];
        double *us;

        /* Compute state u from adjoint w and print to file */
        sprintf(filename1, "%s.%03d", "advec-diff-implicit.out.u0", (app->myid));
        file = fopen(filename1, "w");
        vec_create(mspace, &us);
        vec_copy(mspace, U0, us);
        for (j = 0; j < mspace; j++)
        {
            if(j!=mspace-1)
            {
                fprintf(file, "% 1.14e, ", us[j]);
            }
            else
            {
                fprintf(file, "% 1.14e", us[j]);
            }
        }
        vec_destroy(us);


        //Calculate U from W and print out to file
        sprintf(filename, "%s.%03d", "advec-diff-implicit.out.u", (app->myid));
        file = fopen(filename, "w");
        double **u = (double **)malloc(app->npoints * sizeof(double*));
        for(int i = 0; i < app->npoints; i++) u[i] = (double *)malloc(mspace * sizeof(double));
        vec_create((app->mspace), &v);
        for (i = 0; i < (app->npoints); i++)
        {
            double **w = (app->w);
            vec_copy(mspace, w[i], v);

            if(i!=app->npoints-1)
            {
                vec_copy(mspace, w[i],u[i]);
                apply_Aadjoint(dt,dx,nu,mspace,u[i]);
                vec_scale(mspace,-1.0/(dx*dt),u[i]);
                vec_axpy(mspace,1.0/(dx*dt),w[i+1],u[i]);
                vec_axpy(mspace,1.0,U0,u[i]);
            }
            else
            {
                vec_copy(mspace, w[i],u[i]);
                apply_Aadjoint(dt,dx,nu,mspace,u[i]);
                vec_scale(mspace,-1.0/(dx*dt),u[i]);
                vec_axpy(mspace,1.0,U0,u[i]);
            }

            fprintf(file, "%05d: ", ((app->ilower)+i+1));
            for (j = 0; j < mspace; j++)
            {
                if(j==mspace-1)
                {
                    fprintf(file, "% 1.14e", u[i][j]);
                }
                else
                {
                    fprintf(file, "% 1.14e, ", u[i][j]);
                }
            }
            fprintf(file, "\n");
        }
        fflush(file);
        fclose(file);

      // Calculates value of objective function
      /*double objective_val=0;

      for(int i=0; i<ntime; i++)
      {
         for(int j=0; j<mspace; j++)
         {
            objective_val += ((u[i][j] - U0[j]) * (u[i][j] - U0[j]) + alpha*vs[i][j]*vs[i][j] ) * dx;
         }
         objective_val *= dt;
      }
      printf("Objective Function Value: %f \n", objective_val);
      */
      
   }

   /* Print runtime to file (for runtime comparisons)*/
   time = (double)(end-start)/CLOCKS_PER_SEC;
   printf("Total Run Time: %f s \n", time);
   {
      char    filename[255];
      FILE   *file;
      
      //Note that this out file appends the number of time steps
      sprintf(filename, "%s.%d", "advec-diff-implicit.time", ntime);

      file = fopen(filename, "w");
      fprintf(file, "%f", time);
      fflush(file);
      fclose(file);
   }

   free(app);
   
   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}