#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"
#define PI 3.14159265
#define g(dt,dx) dt/(2*dx)
#define b(dt,dx,nu) nu*dt/(dx*dx)
/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure can contain anything, and be named anything as well */

/* Vector structure can contain anything, and be name anything as well */

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


double
*my_Step(int ntime, int mspace, double nu, double *u, double *uleft, int index)
{
   double dx = 1.0/(mspace-1);
   double dt = 1.0/ntime;

   double *utmp;
   vec_create(mspace, &utmp);
   vec_copy(mspace, u, utmp);

   double A = (nu*dt / (dx*dx));
   double B = 1 - 2*nu*dt/(dx*dx);

   utmp[0] = -b(dt,dx,nu)*u[1]+u[0]*(1+2*b(dt,dx,nu));
   for(int i = 1; i <= mspace-2; i++)
   {
    utmp[i] = -b(dt,dx,nu)*(u[i-1]+u[i+1])+u[i]*(1+2*b(dt,dx,nu));
   }
   utmp[mspace-1] = -b(dt,dx,nu)*(u[mspace-2])+u[mspace-1]*(1+2*b(dt,dx,nu));
   double *u2tmp;
      vec_create(mspace, &u2tmp);
   if (uleft != NULL)
   {
      
      vec_copy(mspace, uleft, u2tmp);
      vec_axpy(mspace, -1.0, u2tmp, utmp);
      utmp[0] = utmp[0]+g(dt,dx)*u2tmp[0]*(utmp[1]);
      for(int i = 1; i <= mspace-2; i++)
      {
        utmp[i] = utmp[i]+g(dt,dx)*u2tmp[i]*(utmp[i+1]-utmp[i-1]);
      }
      utmp[mspace-1] = utmp[mspace-1]+g(dt,dx)*u2tmp[mspace-1]*(utmp[mspace-2]);
   }

   if ((index == 0))
   {
      vec_copy(mspace, u, u2tmp);
      vec_axpy(mspace,-1.0,u2tmp,utmp);
   }

   return utmp;
}

int main (int argc, char *argv[])
{

   double      dt, dx, nu, tstop, tstart; 
   int         ntime, mspace, arg_index;

   /* Define space domain. Space domain is between 0 and 1, mspace defines the number of steps */
   mspace = 16;
   ntime = 512;

   /* Define some optimization parameters */         /* parameter in the objective function */
   nu    = 0.7;                /* parameter in PDE */

   /* Define some Braid parameters */

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
      else if ( strcmp(argv[arg_index], "-nu") == 0 )
      {
         arg_index++;
         nu = atof(argv[arg_index++]);
      }
      else
      {
         printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
         return (0);
      }
   }


   /* Define the space step for dt computation */
   dx=(double)1/(mspace);
   dt = (tstop-tstart)/(double)ntime;

   //vectors below used for implicit solve employed in my_step
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


   double **w = (double **)malloc(ntime * sizeof(double*));
   for(int i = 0; i < ntime; i++) w[i] = (double *)malloc(mspace * sizeof(double));

   for (int i = 0; i <= mspace-1; i++)
   {
      if(i==0){
         w[0][i] = 0.0;
      }
      else if(i<=mspace/2-1){
        w[0][i] = 1;
      }
      else{
        w[0][i] = 0.0;
      }
      
   }
   for(int i=0; i<ntime; i++)
   {
      //my_step returns the residual
      double *rtmp;
      if(i>0){
      rtmp = my_Step(ntime, mspace, nu, w[i], w[i-1],i);
   }
   else{
      rtmp = my_Step(ntime, mspace, nu, w[i], NULL, i);
   }


   //applying c_tilde to the residual
   rtmp[0] = (-1/(1+2*b(dt,dx,nu)))*rtmp[0]; 
   /*printf("1 entry divided by %lf \n", (1+2*b(dt,dx,nu)+g(dt,dx)*(utmp[1])));*/
   for(int i = 1; i <= mspace-2; i++)
   {
    rtmp[i] = (-1/(1+2*b(dt,dx,nu)))*rtmp[i];
    /*printf("%d entry divided by %lf \n", i+1,(1+2*b(dt,dx,nu)+g(dt,dx)*(utmp[i+1]-utmp[i-1])));*/
   }
   rtmp[mspace-1] = (-1/(1+2*b(dt,dx,nu)))*rtmp[mspace-1]; 
   /*printf("%d entry divided by %lf \n\n", mspace, (1+2*b(dt,dx,nu)+g(dt,dx)*(-utmp[mspace-2])));*/

      for(int j=0; j<mspace; j++)
      {
         w[i][j] = w[i][j] + rtmp[j];
      }
      
   }   

   /* Set up the app structure */


   {
         char  filename[255];
         FILE *file;
         int   i,j;

         sprintf(filename, "%s.%03d", "visc-berg-imp.out.u", 000);
         file = fopen(filename, "w");
         for (i = 0; i < ntime; i++)
         {
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
   return (0);
}

