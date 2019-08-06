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
   double   alpha;       /* Relaxation parameter for objective function, v(x,t) */
   double   nu;          /* Diffusion coefficent, which we take to be large */
   int      ntime;       /* Total number of time-steps (starting at time 0) */
   int      mspace;      /* Number of space points included in our state vector */
                         /* So including boundaries we have M+2 space points */

   double **U;           /* This holds the u, v, and w vectors. U is in R^3MN */
   double *U0;
   double *ai;
   double *li;

   double **u_state;      /* u,v,w_state hold the current values. They should be accessile in TriSolve*/
   double **v_state;      /* Each of these should be of size MN */
   double **w_state;


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

/*------------------------------------*/

void
apply_PhiAdjoint(double dt, double dx, double nu, int M, double *u, double *l, double *a)
{
   /* First solve U^Tw=u (U^Tw=f) */

   /* Need to change this w to some other letter becuase w is already passed as a parameter of this function */
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

void add_Gamma(int mspace, double *u_n, double *v)
{
   v[0]+=u_n[0]*u_n[1];
   for (int i = 0; i<mspace-2;i++)
   {
      v[i]+=u_n[i]*(u_n[i+1]);
   }
   v[mspace-1]+=-u_n[mspace-1]*u_n[mspace-2];
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
   //add all arguments to apply_Phi below based on what Isaiah does
   /* apply_Phi(dt, dx, nu, M, v, l, a); */
    for (int i = 0; i <= M-1; i++)
    {
       v[i] *= dt;
    }
}

/*------------------------------------*/

void
apply_DAdjoint(double dt, double dx, double nu, int M, double *v, double *l, double *a)
{
   //add all arguments to apply_PhiAdjoing based on what Isaiah does
   /* apply_PhiAdjoint(dt, dx, nu, M, v, l, a); */
    for (int i = 0; i <= M-1; i++)
    {
       v[i] *= dt;
    }
}

/*------------------------------------*/

int
my_TriResidual(braid_App       app,
               braid_Vector    uleft,
               braid_Vector    uright,
               braid_Vector    f,
               braid_Vector    r,
               braid_Int       homogeneous,
               braid_TriStatus status)
{
   double  t, tprev, tnext, dt, dx, g;
   double  nu = (app->nu);
   double  alpha = (app->alpha);
   double *rtmp, *uold, *vtmp, *wtmp, *wtmp_right, *utmp;
   int     level, index;
   int     mspace = (app->mspace);
   int     ntime = (app->ntime);
   double **U = (app->U);
   double *u0 = (app->U0);   
   
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

   dx = 1/((double)(mspace+1));

   /*Set the g=gamma term*/
   g=g(dt, dx);

   /***************UPDATE U VARIABLE***************/
   if(index < ntime)
   {
      /*Create the vectors needed for U computation*/
      vec_create(mspace, &rtmp);
      vec_create(mspace, &uold);
      vec_create(mspace, &wtmp);
      vec_create(mspace, &wtmp_right);
      vec_create(mspace, &utmp);

      /*Copy the current values into wtmp, wtmp_right, and copy u into uold*/
      vec_copy(mspace, (r->values), uold);
      vec_copy(mspace, U[index+2*ntime], wtmp);
      vec_copy(mspace, U[index+2*ntime+1], wtmp_right);

      if(index == ntime - 1)
      {
         /*Compute effect of the L^Tw term*/
         vec_copy(mspace, wtmp, rtmp);
         apply_Aadjoint(dt, dx, nu, mspace, rtmp);

         vec_copy(mspace, (r->values), utmp);
         vec_axpy(mspace, -1.0, u0, utmp);
         vec_scale(mspace, dt*dx, utmp);

         vec_axpy(mspace, 1.0, utmp, rtmp);

         /*Compute the effect of the gamma term*/
         rtmp[0]=rtmp[0]-(wtmp[0]-wtmp[1])*uold[1]*g/(dx*dt);
         for(int i = 1; i<mspace-2; i++)
         {
            rtmp[i]=rtmp[i]-(uold[i-1]*wtmp[i-1]+(uold[i+1]-uold[i-1])*wtmp[i]-uold[i+1]*wtmp[i+1])*g/(dx*dt);
         }
         rtmp[mspace-1]=rtmp[mspace-1]-(wtmp[mspace-2]-wtmp[mspace-1])*uold[mspace-2]*g/(dx*dt);
      }
      else
       {
         /*Compute effect of the L^Tw term*/
         vec_copy(mspace, wtmp, rtmp);
         apply_Aadjoint(dt, dx, nu, mspace, rtmp);
         vec_axpy(mspace, -1.0, wtmp_right, rtmp);

         vec_copy(mspace, (r->values), utmp);
         vec_axpy(mspace, -1.0, u0, utmp);
         vec_scale(mspace, dt*dx, utmp);

         vec_axpy(mspace, 1.0, utmp, rtmp);

         /*Compute the effect of the gamma term*/
         rtmp[0]=rtmp[0]-(wtmp[0]-wtmp[1])*uold[1]*g/(dx*dt);
         for(int i = 1; i<mspace-2; i++)
         {
            rtmp[i]=rtmp[i]-(uold[i-1]*wtmp[i-1]+(uold[i+1]-uold[i-1])*wtmp[i]-uold[i+1]*wtmp[i+1])*g/(dx*dt);
         }
         rtmp[mspace-1]=rtmp[mspace-1]-(wtmp[mspace-2]-wtmp[mspace-1])*uold[mspace-2]*g/(dx*dt);
      }
      vec_copy(mspace, rtmp, (r->values)); 
      vec_destroy(rtmp);
      vec_destroy(utmp);
      vec_destroy(wtmp_right);
      vec_destroy(wtmp);
      vec_destroy(uold);     
   }

   /***************SOLVE GRAD V EQUATION***************/
   
   if(index>=ntime && index<2*ntime)
   {
      vec_create(mspace, &rtmp);
      vec_scale(mspace, alpha*dx, rtmp);
      vec_axpy(mspace, -1.0, U[index+ntime], rtmp);
      vec_scale(mspace, dt, rtmp);
      
      vec_copy(mspace, rtmp, (r->values));
      vec_destroy(rtmp);
   }

   /***************SOLVE GRAD W EQUATION***************/

   else
   {
      vec_create(mspace, &rtmp);
      vec_create(mspace, &vtmp);
      vec_create(mspace, &utmp);

      if(index == 2*ntime)
      {
         vec_copy(mspace, U[0], utmp);
         vec_copy(mspace, U[0], rtmp);
         apply_A(dt, dx, nu, mspace, rtmp);

         add_Gamma(mspace, U[0], rtmp);

         vec_axpy(mspace, -1.0, u0, rtmp);
         vec_copy(mspace, U[ntime-index], vtmp);
         vec_scale(mspace, -1.0*dt, vtmp);
         vec_axpy(mspace, 1.0, vtmp, rtmp);
      }
      else
      {
         vec_copy(mspace, U[index], rtmp);
         apply_A(dt, dx, nu, mspace, rtmp);
         vec_copy(mspace, U[index-1], utmp);
         vec_axpy(mspace, -1.0, utmp, rtmp);

         vec_copy(mspace, U[index], utmp);

         add_Gamma(mspace, U[0], rtmp);

         vec_axpy(mspace, -1.0, u0, rtmp);
         vec_copy(mspace, U[ntime-index], vtmp);
         vec_scale(mspace, -1.0*dt, vtmp);
         vec_axpy(mspace, 1.0, vtmp, rtmp);
      }
      vec_copy(mspace, rtmp, (r->values));
      vec_destroy(rtmp);
      vec_destroy(vtmp);
   }

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
   double  t, tprev, tnext, dt, dx, g;
   double *utmp, *uold, *vtmp, *wtmp, *wtmp_right;
   int level, index;
   int mspace = (app->mspace);
   int ntime = (app->ntime);
   double alpha = (app->alpha);
   double nu = (app->nu);
   double **U = (app->U);
   double *u0 = (app->U0);
   
   /* Get the time-step size */
   braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
   braid_TriStatusGetLevel(status, &level);
   braid_TriStatusGetTIndex(status, &index);
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

   /*Set the g=gamma term*/
   g=g(dt, dx);

   /***************UPDATE U VARIABLE***************/
   if(index < ntime)
   {
      /*Create the vectors needed for U computation*/
      
      printf("here 1a\n");
      
      vec_create(mspace, &utmp);
      vec_create(mspace, &uold);
      vec_create(mspace, &wtmp);
      vec_create(mspace, &wtmp_right);

      /*Copy the current values into wtmp, wtmp_right, and copy u into uold*/
      vec_copy(mspace, (u->values), uold);
      

         /*if (1==1)
         {
            printf("%f\n", u0[i]);
         }
         else
         {
            printf("Not null\n");
         }*/
         /*vec_copy(mspace, U[1][1], utmp);*/
         printf("here 2a.1\n");
         /*printf("%f\n", U[0][0]);*/
         double t=U[0][0];
         printf("here 2a.2\n");
      


      vec_copy(mspace, U[index+(2*ntime)], wtmp);
      
      vec_copy(mspace, U[index+(2*ntime)+1], wtmp_right);
      printf("here 2a.3\n");

      

      if((int)index == ntime - 1)
      {
         /*Compute effect of the L^Tw term*/
         printf("here 1b\n");
         vec_copy(mspace, wtmp, utmp);
         apply_Aadjoint(dt, dx, nu, mspace, utmp);
         vec_scale(mspace, -1.0/(dx*dt), utmp);
         vec_axpy(mspace, 1.0, u0, utmp);
         /*Compute the effect of the gamma term*/
         utmp[0]=utmp[0]-(wtmp[0]-wtmp[1])*uold[1]*g/(dx*dt);
         for(int i = 1; i<mspace-2; i++)
         {
            utmp[i]=utmp[i]-(utmp[i-1]*wtmp[i-1]+(uold[i+1]-utmp[i-1])*wtmp[i]-uold[i+1]*wtmp[i+1])*g/(dx*dt);
         }
         utmp[mspace-1]=utmp[mspace-1]-(wtmp[mspace-2]-wtmp[mspace-1])*utmp[mspace-2]*g/(dx*dt);
         printf("here 2b\n");
      }
      else
      {
         /*Compute effect of the L^Tw term*/
         printf("here 1c\n");
         vec_copy(mspace, wtmp, utmp);
         apply_Aadjoint(dt, dx, nu, mspace, utmp);
         vec_axpy(mspace, 1.0, wtmp_right, utmp);
         vec_scale(mspace, -1.0/(dx*dt), utmp);
         vec_axpy(mspace, 1.0, u0, utmp);
         /*Compute the effect of the gamma term*/
         utmp[0]=utmp[0]-(wtmp[0]-wtmp[1])*uold[1]*g/(dx*dt);
         for(int i = 1; i<mspace-2; i++)
         {
            utmp[i]=utmp[i]-(utmp[i-1]*wtmp[i-1]+(uold[i+1]-utmp[i-1])*wtmp[i]-uold[i+1]*wtmp[i+1])*g/(dx*dt);
         }
         utmp[mspace-1]=utmp[mspace-1]-(wtmp[mspace-2]-wtmp[mspace-1])*utmp[mspace-2]*g/(dx*dt);
         printf("here 2c\n");
      }
      printf("here 1d\n");
      vec_copy(mspace, utmp, (u->values)); 
      vec_destroy(utmp);
      vec_destroy(wtmp_right);
      vec_destroy(wtmp);
      vec_destroy(uold);     
      printf("here 2d\n");
   }
   /***************UPDATE V VARIABLE***************/

   if(index>=ntime && index<2*ntime)
   {
      vec_create(mspace, &vtmp);


      if(index == ntime)
      {
         printf("here 1e\n");
         vec_copy(mspace, U[0], vtmp);
         apply_A(dt, dx, nu, mspace, vtmp);
         vec_scale(mspace, 1.0/dt, vtmp);
         vec_axpy(mspace, 1.0, u0, vtmp);
         add_Gamma(mspace, U[0], vtmp);
         printf("here 2e\n");
      }
      else
      {
         printf("here 1f\n");
         vec_copy(mspace, U[index], vtmp);
         apply_A(dt, dx, nu, mspace, vtmp);
         vec_axpy(mspace, -1.0, U[index-1], vtmp);
         vec_scale(mspace, 1.0/dt, vtmp);
         add_Gamma(mspace, U[index], vtmp);
         printf("here 2f\n");
      }
      vec_copy(mspace, vtmp, (u->values));
      vec_destroy(vtmp);
   }

   /***************UPDATE W VARIABLE***************/
   
   else
   {
      printf("here 1g\n");
      vec_create(mspace, &wtmp);
      vec_copy(mspace, U[index+ntime], wtmp);
      vec_scale(mspace, alpha*dx, wtmp);
      vec_copy(mspace, wtmp, (u->values));
      vec_destroy(wtmp);
      printf("here 2g\n");
   }

   /***********************************************/
   /* no refinement */
   braid_TriStatusSetRFactor(status, 1);

   return 0;
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
      if ((app->U) == NULL)
      {
         int  ntpoints;
         braid_AccessStatusGetNTPoints(astatus, &ntpoints);
         ntpoints++;  /* ntpoints is really the gupper index */
         (app->U) = (double **) calloc(ntpoints, sizeof(double *));
      }

      braid_AccessStatusGetTIndex(astatus, &index);
      if (app->U[index] != NULL)
      {
         free(app->U[index]);
      }
      vec_create(mspace, &(app->U[index]));
      vec_copy(mspace, (u->values), (app->U[index]));
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
   braid_Core  core;
   my_App     *app;
         
   double      tstart, tstop, dt; 
   int         rank, ntime, mspace, arg_index;
   double      alpha, nu;
   int         max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter;
   int         access_level, print_level;
   double      tol;

   printf("here 1a.m\n");
   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   printf("here 2a.m\n");
   /* Define space domain. Space domain is between 0 and 1, mspace defines the number of steps */
   mspace = 8;
   ntime = 256;

   /* Define some optimization parameters */
   alpha = .005;            /* parameter in the objective function */
   nu    = 2;                /* parameter in PDE */

   /* Define some Braid parameters */
   max_levels     = 30;
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
   /*dx=(double)1/(mspace+1);*/

   /* Define time domain and step */
   tstart = 0.0;             /* Beginning of time domain */
   tstop  = 1.0;             /* End of time domain*/
   dt = (tstop-tstart)/ntime; 
    
   printf("here 1b.m\n");
   /* Set up the app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->mspace   = mspace;
   app->nu       = nu;
   app->alpha    = alpha;
   app->U        = NULL;
   printf("here 2b.m\n");

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

   printf("here 1c.m\n");
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
   printf("here 2c.m\n");

   printf("here 1d.m\n");
   /* Parallel-in-time TriMGRIT simulation */
   braid_Drive(core);
   printf("here 2d.m\n");

   free(app);
   
   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}