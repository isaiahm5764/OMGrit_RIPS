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
 * Example:       visc-imp-omgrit.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make visc-imp-omgrit
 *
 * Description:  Solves a non-linear optimal control problem in time-parallel:
 * 
 *                 min   0.5\int_0^T \int_0^1 (u(x,t)-u0(x))^2+alpha v(x,t)^2 dxdt
 * 
 *                  s.t.  du/dt + u du/dx - nu d^2u/dx^2 = v(x,t)
 *                        u(0,t)=u(1,t)=0
 *                                  u(x,0)=u0(x)
 *
 *             
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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

void
apply_PhiInv(double dt, double dx, double nu, int M, double *u, double *l, double *a)
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
   double b = -b(dt, dx, nu);
   u[M-1]=w[M-1]/a[M-1];
   for (int i = M-2; i >= 0; i--)
   {
      u[i]=(w[i]-b*u[i+1])/a[i];      
   }
}

/*------------------------------------*/

void
apply_PhiAdjointInv(double dt, double dx, double nu, int M, double *u, double *l, double *a)
{
   /* First solve U^Tw=u (U^Tw=f) */

   /* Need to change this w to some other letter becuase w is already passed as a parameter of this function */
   double *w;

   vec_create(M, &w);
   double *f;
   vec_create(M, &f);
   vec_copy(M, u, f);
   double b = -b(dt, dx, nu);
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

/* Some helper functions declared below*/

/* A in the formulation papers. Just a matrix mult. */
void 
apply_Phi(double dt, int mspace, double nu, double *u)
{
  int M=mspace;
  double dx = 1/((double)(mspace+1));
  double A = -b(dt, dx, nu);
  double B = 1+2*b(dt, dx, nu);
  double C = A;
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

void
apply_PhiAdjoint(double dt, int mspace, double nu, double *w)
{
  int M=mspace;
  double dx = 1/((double)(mspace+1));
  double A = -b(dt, dx, nu);
  double B = 1+2*b(dt, dx, nu);
  double C = A;
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

//applies the B matrix but may use previous time steps (not sure which works) uleft is just the u vector that is in that gamma nonlinear term
void 
apply_B(double dt, int mspace, double nu, double *u, double *uleft){
  u[0] = uleft[1] * u[0]  -uleft[1] * u[1];
  for(int i=1; i<mspace-1; i++)
  {
    u[i] = uleft[i-1]*u[i-1] + u[i] * (uleft[i+1] - uleft[i-1])  -uleft[i+1]*u[i+1];
  }
  u[mspace-1] = uleft[mspace-2] * u[mspace-2] + -uleft[mspace-2] * u[mspace-1];
}

void 
find_gamma(double *u, int mspace){
  double *tmp;
  vec_create(mspace, &tmp);
  vec_copy(mspace, u, tmp);

  u[0] = tmp[0]*tmp[1];
  for(int i=1; i<mspace-1; i++)
  {
    u[i] = tmp[i]*tmp[i+1] - tmp[i]*tmp[i-1];
  }
  u[mspace-1] = -tmp[mspace-1] * tmp[mspace-2];
  vec_destroy(tmp);
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
   double *rtmp, *utmp, *u2tmp, *tmp;
   int     level, index;
   int     mspace = (app->mspace);
   double *u0 = (app->U0);
   double alpha = (app->alpha);

   double *l = (app->li);
   double *a = (app->ai);
   
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
   vec_create(mspace, &u2tmp);
   vec_create(mspace, &tmp);

   /* CENTER BLOCK */

   if(uright != NULL){   
    vec_copy(mspace, r->values, utmp);
    vec_copy(mspace, utmp, tmp);

    vec_axpy(mspace, alpha*dx/dt, utmp, rtmp);

    apply_Phi(dt, mspace, nu, utmp);
    vec_scale(mspace, alpha*dx/dt, utmp);
    apply_PhiAdjoint(dt, mspace, nu, utmp);
    vec_axpy(mspace, 1.0, utmp, rtmp);

    vec_copy(mspace, r->values, utmp);
    vec_axpy(mspace, dx*dt, utmp, rtmp);

    apply_Phi(dt, mspace, nu, utmp);
    apply_B(dt, mspace, nu, utmp, tmp);
    vec_axpy(mspace, alpha/2, utmp, rtmp);

   }
   else{
    vec_copy(mspace, r->values, utmp);
    vec_copy(mspace, utmp, tmp);

    apply_Phi(dt, mspace, nu, utmp);
    vec_scale(mspace, alpha*dx/dt, utmp);
    apply_PhiAdjoint(dt, mspace, nu, utmp);
    vec_axpy(mspace, 1.0, utmp, rtmp);

    vec_copy(mspace, r->values, utmp);
    vec_axpy(mspace, dx*dt, utmp, rtmp);

    apply_Phi(dt, mspace, nu, utmp);
    apply_B(dt, mspace, nu, utmp, tmp);
    vec_axpy(mspace, alpha/2, utmp, rtmp);
   }

   //




   /* WEST BLOCK */

   if(uleft != NULL)
   {
    vec_copy(mspace, r->values, tmp);
    
    vec_copy(mspace, uleft->values, utmp);
    apply_PhiAdjoint(dt, mspace, nu, utmp);
    vec_axpy(mspace, -alpha*dx/dt, utmp ,rtmp);

    vec_copy(mspace, uleft->values, utmp);
    apply_B(dt, mspace, nu, utmp, tmp);
    vec_axpy(mspace, -alpha/2, utmp, rtmp);
    
   }


  /* EAST BLOCK */
   if(uright!=NULL){
   vec_copy(mspace, uright->values, utmp);
    apply_Phi(dt, mspace, nu, utmp);
    vec_axpy(mspace, -alpha*dx/dt, utmp, rtmp);
 }





  /* 3 non-linear terms have not been dealt with by C,W,E */
  /* They are included now.                               */

  /* First of 3 */
   if (uleft == NULL){
    vec_copy(mspace, r->values, tmp);
    vec_copy(mspace, u0, u2tmp);
    apply_PhiInv(dt, dx, nu, mspace, u2tmp, l, a);
    apply_B(dt, mspace, nu, u2tmp, tmp);
    vec_axpy(mspace, -alpha/2, u2tmp, rtmp);
   }

  /* Second of 3 */
   if (uright != NULL){
    vec_copy(mspace, r->values, tmp);
    vec_copy(mspace, uright->values, u2tmp);
    find_gamma(u2tmp, mspace);
    find_gamma(tmp, mspace);
    apply_PhiAdjoint(dt, mspace, nu, tmp);
    vec_axpy(mspace, -1, u2tmp, tmp);
    vec_axpy(mspace, alpha/2, tmp, rtmp);
   }
  else {
    vec_copy(mspace, r->values, tmp);
    find_gamma(tmp, mspace);
    apply_PhiAdjoint(dt, mspace, nu, tmp);
    vec_axpy(mspace, alpha/2, tmp, rtmp);
  }

  /* Third of 3 */
  vec_copy(mspace, r->values, tmp);
  vec_copy(mspace, r->values, u2tmp);
  find_gamma(tmp, mspace);
  apply_B(dt, mspace, nu, tmp, u2tmp);
  vec_axpy(mspace, g(dt,dx)*alpha/2, tmp, rtmp);



   //calculate f and subtract to find residual
   if ((!homogeneous))
   {
      vec_axpy(mspace, -dx*dt, u0, rtmp);

      if (uleft == NULL){
        vec_copy(mspace, u0, u2tmp);
        apply_PhiAdjoint(dt, mspace, nu, u2tmp);
        vec_axpy(mspace, -alpha*dx/dt, u2tmp, rtmp);
      } 
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
   vec_destroy(u2tmp);
   vec_destroy(tmp);
   // free(u0);
   // free(l);
   // free(a);

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
   double alpha = (app->alpha);

   // double *l = (app->li);
   // double *a = (app->ai);
   
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
   vec_create(mspace, &rtmp);

   /* Initialize temporary solution vector */
   vec_copy(mspace, (u->values), utmp);
   
   /* Compute residual */
   my_TriResidual(app, uleft, uright, f, u, homogeneous, status);


   rtmp = (u->values);
   if (uright != NULL){
    rtmp[0] = -rtmp[0]/(dx*dt + (alpha*dx/dt)*(2+4*b(dt,dx,nu)+6*(b(dt,dx,nu))*(b(dt,dx,nu))) + alpha*((1+3*b(dt,dx,nu))*(utmp[1])+g(dt,dx)*(utmp[1]*utmp[1])));
    for(int i=1; i<mspace-1; i++){
      rtmp[i]=-rtmp[i]/(dx*dt + (alpha*dx/dt)*(2+4*b(dt,dx,nu)+6*(b(dt,dx,nu))*(b(dt,dx,nu))) + alpha*((1+3*b(dt,dx,nu))*(utmp[i+1]-utmp[i-1])+g(dt,dx)*(utmp[i+1]*utmp[i+1]+utmp[i-1]*utmp[i-1]-utmp[i-1]*utmp[i+1])));
      // printf("Dividing by %f",(dx*dt + (alpha*dx/dt)*(1+2*b(dt,dx,nu)+6*(b(dt,dx,nu))*(b(dt,dx,nu))) + alpha*((1+3*b(dt,dx,nu))*(utmp[i+1]-utmp[i-1])+2*g(dt,dx)*(utmp[i+1]*utmp[i+1]+utmp[i-1]*utmp[i-1]+utmp[i-1]*utmp[i+1]))));
    }
    rtmp[mspace-1] = -rtmp[mspace-1]/(dx*dt + (alpha*dx/dt)*(2+4*b(dt,dx,nu)+6*(b(dt,dx,nu))*(b(dt,dx,nu))) + alpha*((1+3*b(dt,dx,nu))*(-utmp[mspace-2])+g(dt,dx)*(utmp[mspace-2]*utmp[mspace-2])));
    // printf("\n ========================================================================================================= \n");
   }
   else {
    rtmp[0] = -rtmp[0]/(dx*dt + (alpha*dx/dt)*(1+4*b(dt,dx,nu)+6*(b(dt,dx,nu))*(b(dt,dx,nu))) + alpha*((1+3*b(dt,dx,nu))*(utmp[1])+g(dt,dx)*(utmp[1]*utmp[1])));
    for(int i=1; i<mspace-1; i++){
      rtmp[i]=-rtmp[i]/(dx*dt + (alpha*dx/dt)*(1+4*b(dt,dx,nu)+6*(b(dt,dx,nu))*(b(dt,dx,nu))) + alpha*((1+3*b(dt,dx,nu))*(utmp[i+1]-utmp[i-1])+g(dt,dx)*(utmp[i+1]*utmp[i+1]+utmp[i-1]*utmp[i-1]-utmp[i-1]*utmp[i+1])));
      // printf("Dividing by %f",(dx*dt + (alpha*dx/dt)*(1+2*b(dt,dx,nu)+6*(b(dt,dx,nu))*(b(dt,dx,nu))) + alpha*((1+3*b(dt,dx,nu))*(utmp[i+1]-utmp[i-1])+2*g(dt,dx)*(utmp[i+1]*utmp[i+1]+utmp[i-1]*utmp[i-1]+utmp[i-1]*utmp[i+1]))));
    }
    rtmp[mspace-1] = -rtmp[mspace-1]/(dx*dt + (alpha*dx/dt)*(1+4*b(dt,dx,nu)+6*(b(dt,dx,nu))*(b(dt,dx,nu))) + alpha*((1+3*b(dt,dx,nu))*(-utmp[mspace-2])+g(dt,dx)*(utmp[mspace-2]*utmp[mspace-2])));
    // printf("\n ========================================================================================================= \n");
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

   /* prints U, V, and W after selected iterations. This can then be plotted to show how the space-time solution changes after iterations. */

     char  filename[255];
     FILE *file;
     int  iter;
     braid_AccessStatusGetIter(astatus, &iter);
     braid_AccessStatusGetTIndex(astatus, &index);
     /* file format is advec-diff-btcs.out.{iteration #}.{time index} */
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
         
   double      tstart, tstop, dt, dx; 
   int         rank, ntime, mspace, arg_index;
   double      alpha, nu;
   int         max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter;
   int         access_level, print_level;
   double      tol;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   /* Define space domain. Space domain is between 0 and 1, mspace defines the number of steps */
   mspace = 8;
   ntime = 256;

   /* Define some optimization parameters */
   alpha = .005;            /* parameter in the objective function */
   nu    = 2;                /* parameter in PDE */

   /* Define some Braid parameters */
   max_levels     = 30;
   min_coarse     = 1;
   nrelax         = 30;
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


   /* Define time domain and step */
   tstart = 0.0;             /* Beginning of time domain */
   tstop  = 1.0;             /* End of time domain*/
   dt = (tstop-tstart)/ntime; 
    

   /* Set up the app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->mspace   = mspace;
   app->nu       = nu;
   app->alpha    = alpha;
   app->w        = NULL;

   /* Set this to whatever u0 is. */

   double *U0 = (double*) malloc( ntime*sizeof(double) );
   for(int i=0; i<mspace/2; i++){
      U0[i]=1;
   }

   for(int i=mspace/2; i<mspace; i++)
   {
      U0[i]=0;
   }

   app->U0       = U0;

   dx = 1/((double)(mspace+1));

   double *ai = (double*) malloc( mspace*sizeof(double) );
   double *li = (double*) malloc( (mspace-1)*sizeof(double) );
   double A = -b(dt,dx,nu);
   double B = 1 +2*b(dt,dx,nu);
   double C = A;
   ai[0] = B;
   for(int i=1; i<mspace; i++){
      li[i-1] = A/ai[i-1];
      ai[i] = ai[0]+(-C)*li[i-1];
   }

   app->ai = ai;
   app->li = li;



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

      /* Print w to file, w is actually u, the notation is just carried over from the 
      * optimization case. */

      {
         char  filename[255];
         FILE *file;
         int   i,j;

         sprintf(filename, "%s.%03d", "out/viscous-burgers-2.w", (app->myid));
         file = fopen(filename, "w");
         for (i = 0; i < (app->ntime); i++)
         {
            double **w = (app->w);
            fprintf(file, "%05d: ", (i+1));
            for(j=0; j <mspace; j++){
               if(j==mspace-1){
                  fprintf(file, "% 1.14e", w[i][j]);
               }
               else{
                  fprintf(file, "% 1.14e, ", w[i][j]);
               }
            }
            fprintf(file, "\n");
         }
         fflush(file);
         fclose(file);
      }
      
   }



   free(app);
   
   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}