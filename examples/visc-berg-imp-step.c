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
 * Example:       ex-01.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-01
 *
 * Help with:     this is the simplest example available, read the source
 *
 * Sample run:    mpirun -np 2 ex-01
 *
 * Description:   solve the scalar ODE 
 *                   u' = lambda u, 
 *                   with lambda=-1 and y(0) = 1
 *                in a very simplified XBraid setting.
 *                
 *                When run with the default 10 time steps, the solution is:
 *                $ ./ex-01
 *                $ cat ex-01.out.00*
 *                  1.00000000000000e+00
 *                  6.66666666666667e-01
 *                  4.44444444444444e-01
 *                  2.96296296296296e-01
 *                  1.97530864197531e-01
 *                  1.31687242798354e-01
 *                  8.77914951989026e-02
 *                  5.85276634659351e-02
 *                  3.90184423106234e-02
 *                  2.60122948737489e-02
 *                  1.73415299158326e-02
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"
#define g(dt,dx) dt/(2*dx)
#define b(dt,dx,nu) (double)nu*dt/(dx*dx)

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
   int      myid;
   double   alpha;       /* Relaxation parameter for objective function, v(x,t) */
   double   nu;          /* Diffusion coefficent, which we take to be large */
   int      ntime;       /* Total number of time-steps (starting at time 0) */
   int      mspace;      /* Number of space points included in our state vector */
                         /* So including boundaries we have M+2 space points */

   double **w; 
   double *ai;
   double *li;
} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
   double *values;
} my_Vector;

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

void
vec_copy(int size, double *invec, double *outvec)
{
   int i;
   for (i = 0; i < size; i++)
   {
      outvec[i] = invec[i];
   }
}

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
  int index;
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetTIndex(status, &index);
   int mspace = (app->mspace);
   int ntime = (app->ntime);
   double nu = (app->nu);
   double *l = (app->li);
   double *a = (app->ai);
   double dx = 1.0/(mspace-1);
   double dt = tstop - tstart;

   if(index==0){
    double *utmp;
   
    vec_create(mspace, &utmp);
    vec_copy(mspace, u->values, utmp);
    utmp[0] = -1*(g(dt,dx)*u->values[0]*u->values[1] - utmp[0]);
    for(int i=1; i<mspace-1; i++)
    {
      utmp[i] = -1*(g(dt,dx)*u->values[i]*(u->values[i+1] - u->values[i-1]) - utmp[i]);
    }
    utmp[mspace-1] = -1*(g(dt,dx)*u->values[mspace-2]*u->values[mspace-3] - utmp[mspace-1]);
    // double *w;
    //  vec_create(mspace, &w);
    //  double *f;
    //  vec_create(mspace, &f);
    //  vec_copy(mspace, utmp, f);
    //  w[0]=f[0];
    //  for (int i = 1; i < mspace; i++)
    //  {
    //     w[i]=f[i]-l[i-1]*w[i-1];
    //  }

    //  /* Now solve Ux=w */ 
    //  double b = -b(dt,dx,nu);
    //  utmp[mspace-1]=w[mspace-1]/a[mspace-1];
    //  for (int i = mspace-2; i >= 0; i--)
    //  {
    //     utmp[i]=(w[i]-b*utmp[i+1])/a[i];      
    //  }

    for(int i = 1; i <= mspace-2; i++)
     {
      utmp[i] = (1/(1+2*b(dt,dx,nu)))*utmp[i];
      //printf("%d entry divided by %lf \n", i+1,(1+2*b(dt,dx,nu)+g(dt,dx)*(utmp[i+1]-utmp[i-1])));
     }
     utmp[mspace-1] = (1/(1+2*b(dt,dx,nu)))*utmp[mspace-1]; 

       vec_copy(mspace, utmp, u->values);
     }
   else{
     double *utmp;
   
    vec_create(mspace, &utmp);
    vec_copy(mspace, u->values, utmp);
    utmp[0] = g(dt,dx)*u->values[0]*u->values[1] - utmp[0];
    for(int i=1; i<mspace-1; i++)
    {
      utmp[i] = g(dt,dx)*u->values[i]*(u->values[i+1] - u->values[i-1]) - utmp[i];
    }
    utmp[mspace-1] = g(dt,dx)*u->values[mspace-2]*u->values[mspace-3] - utmp[mspace-1];
    // double *w;
    //  vec_create(mspace, &w);
    //  double *f;
    //  vec_create(mspace, &f);
    //  vec_copy(mspace, utmp, f);
    //  w[0]=f[0];
    //  for (int i = 1; i < mspace; i++)
    //  {
    //     w[i]=f[i]-l[i-1]*w[i-1];
    //  }

    //  /* Now solve Ux=w */ 
    //  double b = -b(dt,dx,nu);
    //  utmp[mspace-1]=w[mspace-1]/a[mspace-1];
    //  for (int i = mspace-2; i >= 0; i--)
    //  {
    //     utmp[i]=(w[i]-b*utmp[i+1])/a[i];      
    //  }
    for(int i = 1; i <= mspace-2; i++)
   {
    utmp[i] = (-1/(1+2*b(dt,dx,nu)))*utmp[i];
    //printf("%d entry divided by %lf \n", i+1,(1+2*b(dt,dx,nu)+g(dt,dx)*(utmp[i+1]-utmp[i-1])));
   }
   utmp[mspace-1] = (-1/(1+2*b(dt,dx,nu)))*utmp[mspace-1]; 

     vec_copy(mspace, utmp, u->values);
   }
  
   
   return 0;
}

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
      if(i==0){
         u->values[i] = 0.0;
      }
      else if(i<=mspace/2-1){
        u->values[i] = 1.0;
      }
      else{
        u->values[i] = 0.0;
      }
      
   }

   *u_ptr = u;

   return 0;
   // my_Vector *u;
   // int mspace = (app->mspace);

   // /* Allocate the vector */
   // u = (my_Vector *) malloc(sizeof(my_Vector));
   // vec_create(mspace, &(u->values));

   // for (int i = 0; i <= mspace-1; i++)
   // {
   //    u->values[i] = ((double)braid_Rand())/braid_RAND_MAX;
   // }

   // *u_ptr = u;

   // return 0;
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

   // char  filename[255];
   //   FILE *file;
   //   int  iter;
   //   braid_AccessStatusGetIter(astatus, &iter);
   //   braid_AccessStatusGetTIndex(astatus, &index);
   //   /* file format is advec-diff-btcs.out.{iteration #}.{time index} */
   //   if(iter%1==0){
   //      sprintf(filename, "%s.%04d.%04d", "out/advec-imp-iter.u.out", iter, index);
   //      file = fopen(filename, "w");
   //      for(int i = 0; i<mspace; i++){
   //          if(i<mspace-1){
   //             fprintf(file, "%1.14e, ", (u->values)[i]);
   //          }
   //          else{
   //             fprintf(file, "%1.14e", (u->values)[i]);
   //          }
   //      }
   //   fflush(file);
   //   fclose(file);
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

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;

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

   /* Define some optimization parameters */
   alpha = .005;            /* parameter in the objective function */
   nu    = 2;                /* parameter in PDE */

   /* Define some Braid parameters */
   max_levels     = 5;
   min_coarse     = 1;
   nrelax         = 1;
   nrelaxc        = 7;
   maxiter        = 300;
   cfactor        = 2;
   tol            = 1.0e-6;
   access_level   = 2;
   print_level    = 2;

   max_levels     = 30;
   min_coarse     = 1;
   nrelax         = 5;
   nrelaxc        = 5;

      /* Define time domain */
   tstart = 0.0;             /* Beginning of time domain */
   tstop  = 1.0;             /* End of time domain*/

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


   /* Define the space step for dt computation */
   dx=(double)1/(mspace+1);



   dt = (tstop-tstart)/(double)ntime;   

   /* Set up the app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->mspace   = mspace;
   app->nu       = nu;
   app->alpha    = alpha;
   app->w        = NULL;

   double *ai = (double*) malloc( mspace*sizeof(double) );
   double *li = (double*) malloc( (mspace-1)*sizeof(double) );
   double A = (-nu*dt / (dx*dx));
   double B = 1 + 2*nu*dt/(dx*dx);
   double C = (-nu*dt / (dx*dx));
   ai[0] = B;
   for(int i=1; i<mspace; i++){
      li[i-1] = A/ai[i-1];
      ai[i] = ai[0]+(-C)*li[i-1];
   }

   app->ai = ai;
   app->li = li;
   
   /* initialize XBraid and set options */
   braid_Init(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
   braid_SetMaxLevels(core, max_levels);
   braid_SetMinCoarse(core, min_coarse);
   braid_SetNRelax(core, -1, nrelax);
   if (max_levels > 1)
   {
      braid_SetNRelax(core, max_levels-1, nrelaxc); /* nrelax on coarsest level */
   }
   
   /* Set some typical Braid parameters */
   braid_SetPrintLevel( core, 2);
   braid_SetMaxLevels(core, max_levels);
   braid_SetAbsTol(core, 1.0e-06);
   braid_SetCFactor(core, -1, 2);
   braid_SetMaxIter(core, maxiter);
   braid_SetAccessLevel(core, access_level);
   
   /* Run simulation, and then clean up */
   braid_Drive(core);

   {
         char  filename[255];
         FILE *file;
         int   i,j;

         sprintf(filename, "%s.%03d", "visc-berg-imp-step.out.u", (app->myid));
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

   braid_Destroy(core);
   free(app);
   MPI_Finalize();

   return (0);
}
