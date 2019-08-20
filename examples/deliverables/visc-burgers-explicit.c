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
 * Example:       visc-burgers-exp-rms.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make visc-burgers-exp-rms
 *
 * Description:  Solves a nonlinear optimal control problem in time-parallel:
 * 
 *                 min   0.5\int_0^T \int_0^1 (u(x,t)-u0(x))^2+alpha v(x,t)^2 dxdt
 * 
 *                  s.t.  du/dt + u*du/dx - nu d^2u/dx^2 = v(x,t)
 *                        u(0,t)=u(1,t)=0
 *                                  u(x,0)=u0(x)
 * Visualization: Use the visualization file provided in the /examples folder. 
 * Run the command python viz.py -help for more info.     
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

  double ***w;          /* Vectors (ui, vi, wi) at each time point (each of these have their own mspace points) */
  double *U0;
  double *ai;
  double *li;

  int ilower;
  int iupper;
  int npoints;

} my_App;


/* Define the state vector at one time-step */
typedef struct _braid_Vector_struct
{
  double **values;     /* Holds the R^M state vector (u_1, u_2,...,u_M) for ui vi and wi*/
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

/* B inverse and C inverse are actually the inverses of (I - gamma*B^n) and (dx*dt + gamma*C^n) 
 * as required in the solve */

void
apply_B_inverse(double dt, double dx, double nu, int M, double *u, double *r)
{   
  /* Find LU decomposition */
  double *ai = (double*) malloc( M*sizeof(double) );
  double *gi = (double*) malloc( (M)*sizeof(double) );
  double *g = (double*) malloc( (M)*sizeof(double) );
  double l = dt/(4*dx);
  ai[0] = 1;
  gi[0] = l*2*u[0]/ai[0];
  for(int i=1; i<M; i++)
  {
    ai[i] = 1-(l*2*u[i])*gi[i-1];
    gi[i] = -l*2*u[i]/ai[i];
  }
   
  g[0] = r[0]/ai[0];
  for(int i=1;i<M;i++)
  {
    g[i] = (r[i] - (l*2*u[i])*g[i-1])/ai[i];
  }

  r[M-1] = g[M-1];
  for(int i=M-2; i>=0; i--)
  {
    r[i] = g[i]-gi[i]*r[i+1];
  }

  free(ai);
  free(gi);
  free(g);
}

/* Applies inverse of (dx*dt - gamma*C^n) */

void
apply_C_inverse(double dt, double dx, double nu, int M, double *u, double *r)
{  
  /* Find LU decomposition */
  double *ai = (double*) malloc( M*sizeof(double) );
  double *li = (double*) malloc( (M-1)*sizeof(double) );
  double *bi = (double*) malloc( (M-1)*sizeof(double) );
  double l = dt/(4*dx);
   
  ai[0] = dx*dt;
  for(int i=1; i<M-1; i++)
  {
    bi[i-1] = -l*(u[i-1]-u[i]);
    li[i-1] = -l*(u[i-1]-u[i])/ai[i-1];
    ai[i] = dx*dt+l*l*(u[i-1]-u[i])*(u[i-1]-u[i])/ai[i-1];
  }


  /* Solve Lw=u (Lw=f) */
  double *w;
  vec_create(M, &w);
  double *f;
  vec_create(M, &f);
  vec_copy(M, r, f);
  w[0]=f[0];
  for (int i = 1; i < M; i++)
  {
    w[i]=f[i]-li[i-1]*w[i-1];
  }

  /* Now solve Ux=w */ 
  r[M-1]=w[M-1]/ai[M-1];
  for (int i = M-2; i >= 0; i--)
  {
    r[i]=(w[i]-bi[i]*r[i+1])/ai[i];      
  }

  vec_destroy(w);
  vec_destroy(f);
  free(ai);
  free(li);
  free(bi);
}

/* Applies H^n which is dxdtI - lambda*C^n */
void
apply_C(double dt, double dx, double nu, int M, double *w, double *r)
{  
  double *rold;
  vec_create(M, &rold);
  vec_copy(M, r, rold);
  double l = dt/(4*dx);
  r[0]=(dx*dt - l*(2*w[0])) * rold[0];
  r[M-1]= (dx*dt - l*(-2*w[M-2])) * rold[M-1];
  for(int i = 1; i <= M-2; i++)
  {
    r[i]=  (dx*dt - l*(-2*w[i-1] +2*w[i+1])) * rold[i];
  }
  vec_destroy(rold);
}

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
  double b = b(dt, dx, nu);
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
  double b = b(dt, dx, nu);
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
  double A = b(dt,dx,nu);
  double B = 1-2*b(dt,dx,nu);
  double C = b(dt,dx,nu);
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
  double A = b(dt,dx,nu);
  double B = 1-2*b(dt,dx,nu);
  double C = b(dt,dx,nu);
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

//applies the B matrix but may use previous time steps. uleft is just the u vector that is in the gamma nonlinear term
void 
apply_B(double dt, int mspace, double nu, double *u, double *uleft)
{
  double *tmp;
  vec_create(mspace, &tmp);
  vec_copy(mspace, u, tmp);
  u[0] = tmp[1] * 2*uleft[0];
  for(int i=1; i<mspace-1; i++)
  {
    u[i] = -2*uleft[i]*tmp[i-1] + 2*uleft[i]*tmp[i+1];
  }
  u[mspace-1] = -2*uleft[mspace-1]*tmp[mspace-1];
}

void 
find_gamma(double *u, int mspace){
  double *tmp;
  vec_create(mspace, &tmp);
  vec_copy(mspace, u, tmp);

  u[0] = -tmp[1]*tmp[1];
  for(int i=1; i<mspace-1; i++)
  {
    u[i] = tmp[i-1]*tmp[i-1] - tmp[i+1]*tmp[i+1];
  }
  u[mspace-1] = -tmp[mspace-2] * tmp[mspace-2];

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
  double  alpha = (app->alpha);
  double *rtmp,*rtmp2,*rtmp3, *utmp, *utmp2;
  int     level, index;
  int     mspace = (app->mspace);
  double *u0 = (app->U0);
  double nu = (app->nu);

   
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
  double l = dt/(4*dx);

  /* Create temporary vectors */
  vec_create(mspace, &rtmp);
  vec_create(mspace, &rtmp2);
  vec_create(mspace, &rtmp3);
  vec_create(mspace, &utmp);
  vec_create(mspace, &utmp2);

  /* Compute residual on first row*/   
  vec_copy(mspace, r->values[0], rtmp);
  vec_scale(mspace, dx*dt, rtmp);
  vec_axpy(mspace, -dx*dt, u0, rtmp);
  vec_axpy(mspace, 1.0, r->values[2], rtmp);
  if(uright!=NULL)
  {
    vec_copy(mspace, uright->values[2], utmp);
    apply_Aadjoint(dt,dx,nu,mspace,utmp);
    vec_axpy(mspace, -1.0, utmp, rtmp);
  }
  if(uleft!=NULL)
  {
    vec_copy(mspace, r->values[2], utmp);
    apply_B(dt,mspace,nu,utmp,uleft->values[0]);
    vec_axpy(mspace, -l, utmp, rtmp);
  }


  /* Compute residual on second row*/
  vec_copy(mspace, r->values[1],rtmp2);
  vec_scale(mspace, alpha*dx*dt, rtmp2);
  vec_axpy(mspace, -dt, r->values[2], rtmp2);


  /* Compute residual on third row*/
  vec_copy(mspace, r->values[0], rtmp3);
  if(uleft!=NULL)
  {
    vec_copy(mspace, uleft->values[0], utmp);
    apply_A(dt,dx,nu,mspace,utmp);
    vec_axpy(mspace, -1.0, utmp, rtmp3);
    vec_copy(mspace, uleft->values[0], utmp);
    find_gamma(utmp,mspace);
    vec_axpy(mspace, -l, utmp, rtmp3);
  }
  else
  {
    vec_copy(mspace, u0, utmp);
    apply_A(dt,dx,nu,mspace,utmp);
    vec_axpy(mspace, -1.0, utmp, rtmp3);
    vec_copy(mspace, u0, utmp);
    find_gamma(utmp,mspace);
    vec_axpy(mspace, -l, utmp, rtmp3);
  }
  vec_axpy(mspace, -dt, r->values[1], rtmp3);

  /* FAS */
  if (f != NULL)
  {
    /* rtmp = rtmp - f */
    vec_axpy(mspace, -1.0, (f->values[0]), rtmp);
    vec_axpy(mspace, -1.0, (f->values[1]), rtmp2);
    vec_axpy(mspace, -1.0, (f->values[2]), rtmp3);
  }


  /* Copy temporary residual vector into residual */
  vec_copy(mspace, rtmp, (r->values[0]));
  vec_copy(mspace, rtmp2, (r->values[1]));
  vec_copy(mspace, rtmp3, (r->values[2]));
   
  /* Destroy temporary vectors */
  vec_destroy(rtmp);
  vec_destroy(rtmp2);
  vec_destroy(rtmp3);
  vec_destroy(utmp);
  vec_destroy(utmp2);
   
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
   double *utmp, *utmp2, *r1, *r2, *r3;
   int mspace = (app->mspace);
   double nu = (app->nu);
   //double *li = (app->li);
   //double *ai = (app->ai);
   double alpha = (app->alpha);

   double *dW, *dU, *dV, *storage1, *storage2, *storage3, *storageU;
   vec_create(mspace, &dW);
   vec_create(mspace, &dU);
   vec_create(mspace, &dV);
   for(int i=0; i<mspace; i++) dW[i]=0.0;
   for(int i=0; i<mspace; i++) dV[i]=0.0;
   for(int i=0; i<mspace; i++) dU[i]=0.0;

   
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
  //index 0 refers to U, 1 to V, 2 to W

  /* Get the space-step size */
  dx = 1/((double)(mspace+1));
  vec_create(mspace, &storage1);
  vec_create(mspace, &storage2);
  vec_create(mspace, &storage3);
  vec_create(mspace, &storageU);
  vec_copy(mspace, (u->values)[0], storage1);
  vec_copy(mspace, (u->values)[1], storage2);
  vec_copy(mspace, (u->values)[2], storage3);
  if(uleft!=NULL)
  {
    vec_copy(mspace, uleft->values[0], storageU);
  }

  /* Create temporary vectors */
  vec_create(mspace, &utmp);
  vec_create(mspace, &utmp2);
  vec_create(mspace, &r1);
  vec_create(mspace, &r2);
  vec_create(mspace, &r3);

  //call our residual routine
  my_TriResidual(app, uleft, uright, f, u, homogeneous, status);

  //copy residual vectors into temp vectors
  vec_copy(mspace, u->values[0], r1);
  vec_copy(mspace, u->values[1], r2);
  vec_copy(mspace, u->values[2], r3);

  /*solve for deltaW*/
  vec_copy(mspace,r1,dW);
  vec_axpy(mspace, -dx*dt, r3, dW);
  vec_axpy(mspace, -dt/alpha, r2, dW);


  //apply c_tilde inverse
  if(uleft==NULL)
  {
    vec_scale(mspace, 1.0/(1+(dt*dt)/alpha), dW);
  }
  else
  {
    vec_scale(mspace, 1.0,dW);
    apply_B_inverse(dt,dx,nu,mspace,storageU,dW);
  }


  //update dU and dV based on dW
  //dV
  vec_copy(mspace, r2, dV);
  vec_axpy(mspace, dt, dW, dV);
  vec_scale(mspace, 1.0/(alpha*dx*dt), dV);

  //dU
  vec_copy(mspace, r3, dU);
  vec_axpy(mspace, dt, dV,dU);

  /* Complete update of solution */
  vec_axpy(mspace, -1.0, dU, storage1);
  vec_axpy(mspace, -1.0, dV, storage2);
  vec_axpy(mspace, -1.0, dW, storage3);

  /* Copy solution back into u */
  vec_copy(mspace, storage1, u->values[0]);
  vec_copy(mspace, storage2, u->values[1]);
  vec_copy(mspace, storage3, u->values[2]);
 
  /* Destroy temporary vectors */
  vec_destroy(r1);
  vec_destroy(r2);
  vec_destroy(r3);
  vec_destroy(utmp);
  vec_destroy(storage1);
  vec_destroy(storage2);
  vec_destroy(storage3);
  vec_destroy(dU);
  vec_destroy(dV);
  vec_destroy(dW);

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
  int mspace = (app->mspace);

  /* Allocate the vector */
  my_Vector *u;
  u= (my_Vector *) malloc(sizeof(my_Vector));
  //one extra component allocated for residual containing 4 entries
  u->values = (double**) malloc(3*sizeof(double*) );
  for(int i=0; i<3; i++) u->values[i] = (double *) malloc(mspace * sizeof(double));

  /* Start with random initial guesses for ui, vi, and wi*/  
  for (int i = 0; i <= mspace-1; i++)
  {
    u->values[0][i] = ((double)braid_Rand())/braid_RAND_MAX;
    u->values[1][i] = ((double)braid_Rand())/braid_RAND_MAX;
    u->values[2][i] = ((double)braid_Rand())/braid_RAND_MAX;
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

  v->values = (double**) malloc(3*sizeof(double*) );
  for(int i=0; i<3; i++) v->values[i] = (double *) malloc(mspace * sizeof(double));

  /* Clone the values */
  for (int i = 0; i <= mspace-1; i++)
  {
    v->values[0][i] = u->values[0][i];
    v->values[1][i] = u->values[1][i];
    v->values[2][i] = u->values[2][i];
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
    (y->values)[0][i] = alpha*(x->values)[0][i] + beta*(y->values)[0][i];
    (y->values)[1][i] = alpha*(x->values)[1][i] + beta*(y->values)[1][i];
    (y->values)[2][i] = alpha*(x->values)[2][i] + beta*(y->values)[2][i];
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
    dot += (u->values)[0][i]*(u->values)[0][i];
    dot += (u->values)[1][i]*(u->values)[1][i];
    dot += (u->values)[2][i]*(u->values)[2][i];
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
  int   done, index, i, j, ii;
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
      (app->w) = (double ***)malloc(ntpoints*sizeof(double**));
      for (i = 0; i < ntpoints; i++)
      {
        app->w[i] = (double **) malloc(3*sizeof(double *));
        for (j = 0; j < 3; j++)
        {
          app->w[i][j] = (double *)malloc(mspace*sizeof(double));
        }
      }
    }
    braid_AccessStatusGetILowerUpper(astatus, &(app->ilower), &(app->iupper));
    (app->npoints) = (app->iupper) - (app->ilower) + 1;
    braid_AccessStatusGetTIndex(astatus, &index);
    ii = index - (app->ilower);
    app->w[ii][0] = u->values[0];
    app->w[ii][1] = u->values[1];
    app->w[ii][2] = u->values[2];
  }

  return 0;
}

/*------------------------------------*/

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
  int mspace = app->mspace;
  *size_ptr = mspace*3*sizeof(double);
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
    dbuffer[i] = (u->values)[0][i];
    dbuffer[mspace + i] = (u->values)[1][i];
    dbuffer[2*mspace + i] = (u->values)[2][i];
  }

  braid_BufferStatusSetSize( bstatus,  mspace*3*sizeof(double));

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
  u= (my_Vector *) malloc(sizeof(my_Vector));
  //one extra component allocated for residual containing 4 entries
  u->values = (double**) malloc(3*sizeof(double*) );
  for(int i=0; i<3; i++) u->values[i] = (double *) malloc(mspace * sizeof(double));

  /* Unpack the buffer */
  for(i = 0; i < mspace; i++)
  {
    (u->values)[0][i] = dbuffer[i];
    (u->values)[1][i] = dbuffer[mspace+i];
    (u->values)[2][i] = dbuffer[2*mspace+i];
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
  start=clock();

  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Define space domain. Space domain is between 0 and 1, mspace defines the number of steps */
  mspace = 16;
  ntime = 256;

  /* Define some optimization parameters */
  alpha = .1;            /* parameter in the objective function */
  nu    = .1;              /* parameter in PDE */

  /* Define some Braid parameters */
  max_levels     = 2;
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
  dx=(double)1/(mspace+1);

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

  /* Set this to u0 in the problem formulation */
  double *U0 = (double*) malloc( ntime*sizeof(double) );
  for(int i=0; i<mspace/2; i++)
  {
    U0[i]=1;
  }

  for(int i=mspace/2; i<mspace; i++)
  {
    U0[i]=0;
  }

  app->U0       = U0;

  /* Find elements of LU decomposition of A */
  double *ai = (double*) malloc( mspace*sizeof(double) );
  double *li = (double*) malloc( (mspace-1)*sizeof(double) );
  ai[0] = 1-2*b(dt,dx,nu);
  for(int i=1; i<mspace; i++)
  {
    li[i-1] = -(-b(dt,dx,nu))/ai[i-1];
    ai[i] = ai[0]+(-b(dt,dx,nu))*li[i-1];
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
  braid_Drive(core);
  end=clock();

   /********* Print out if converge or diverge *********/
   double ***w = (app->w);
   double check = w[1][0][1];
   if(isinf(check)||isnan(check))
   {
    check=0;
   }
   else
   {
    check=1;
   }

   char    filename[255];
   FILE   *file;
   sprintf(filename, "%s.%d.%d.%f.%f.%d", "out/visc-burgers-exp-rms.conv", ntime,mspace,nu,alpha,max_levels);
   file = fopen(filename, "w");
   fprintf(file, "%f", check);
   fflush(file);
   fclose(file);   

   /****************************************************/   


  dx = 1/((double)(mspace+1));;
  
  /* Writes final solution to files */
  if (access_level > 0)
  {
    char  filename[255];
    FILE *file;
    int   i,j,index;

    /* Compute state u from adjoint w and print to file */
    sprintf(filename, "%s.%03d", "out/visc-burgers-exp-rms.out.u", (app->myid));
    file = fopen(filename, "w");
    for (i = 0; i < (app->npoints); i++)
    {
      index = (app->ilower) + i +1;
      fprintf(file, "%05d: ", index);
      for(j=0; j <mspace; j++)
      {
        if(j==mspace-1)
        {
          fprintf(file, "% 1.14e", w[i][0][j]);
        }
        else
        {
          fprintf(file, "% 1.14e, ", w[i][0][j]);
        }
      }
      fprintf(file, "\n");
    }
    fflush(file);
    fclose(file);

    char filename1[255]; 
    double *us;

    sprintf(filename1, "%s.%03d", "out/visc-burgers-exp-rms.out.u0", (app->myid));
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

    /* Compute control v from adjoint w and print to file */
    sprintf(filename, "%s.%03d", "out/visc-burgers-exp-rms.out.v", (app->myid));
    file = fopen(filename, "w");
    for (i = 0; i < (app->npoints); i++)
    {
      double ***w = (app->w);
      index = (app->ilower) + i +1;
      fprintf(file, "%05d: ", index);
      for(j=0; j <mspace; j++)
      {
        if(j==mspace-1)
        {
          fprintf(file, "% 1.14e", w[i][1][j]);
        }
        else
        {
          fprintf(file, "% 1.14e, ", w[i][1][j]);
        }
      }
      fprintf(file, "\n");
    }
    fflush(file);
    fclose(file);

    /* Print adjoint w to file */
    sprintf(filename, "%s.%03d", "out/visc-burgers-exp-rms.out.w", (app->myid));
    file = fopen(filename, "w");
    for (i = 0; i < (app->npoints); i++)
    {
      double ***w = (app->w);
      index = (app->ilower) + i +1;
      fprintf(file, "%05d: ", index);
      for(j=0; j <mspace; j++)
      {
        if(j==mspace-1)
        {
          fprintf(file, "% 1.14e", w[i][2][j]);
        }
        else
        {
          fprintf(file, "% 1.14e, ", w[i][2][j]);
        }
      }
      fprintf(file, "\n");
    }
    fflush(file);
    fclose(file);
  }
  
  /* Print runtime to file */
  time = (double)(end-start)/CLOCKS_PER_SEC;
  printf("Total Run Time: %f s \n", time);
  {
    char    filename[255];
    FILE   *file;
    sprintf(filename, "%s.%d", "out/visc-burgers-exp-rms.time", ntime);
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
