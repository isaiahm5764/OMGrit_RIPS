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

/** \file braid.c
 * \brief Source code for user interface routines.  See braid.h for more information.
 *
 */

#include "_braid.h"
#include "_util.h"

#ifndef DEBUG
#define DEBUG 0
#endif

braid_Int
braid_Drive(braid_Core  core)
{
   MPI_Comm             comm_world      = _braid_CoreElt(core, comm_world);
   braid_Int            myid            = _braid_CoreElt(core, myid_world);
   braid_Int            warm_restart    = _braid_CoreElt(core, warm_restart);
   braid_Int            print_level     = _braid_CoreElt(core, print_level);
   braid_App            app             = _braid_CoreElt(core, app);
   braid_Int            obj_only        = _braid_CoreElt(core, obj_only);
   braid_Int            adjoint         = _braid_CoreElt(core, adjoint);
   braid_Real           tstart0         = _braid_CoreElt(core, tstart);
   braid_Real           tstop0          = _braid_CoreElt(core, tstop);
   braid_Int            ntime0          = _braid_CoreElt(core, ntime);
   braid_Int            nchunks         = _braid_CoreElt(core, nchunks);
   _braid_Grid        **grids           = _braid_CoreElt(core, grids);

   braid_Int      ichunk;
   _braid_Grid   *grid;
   braid_Real    *ta;
   braid_Real     dt_chunk;
   braid_Real     localtime, globaltime;
   braid_Int      level, ilower, iupper, i;

   /* Check for non-supported features */
   _braid_FeatureCheck(core);

   /* Some initial output */
   if (myid == 0 )
   {
      if (!warm_restart && print_level > 0)
      {
         _braid_printf("\n  Braid: Begin simulation, time [%f, %f], %d time steps, %d chunks\n",
                       tstart0, tstop0, ntime0, nchunks);
      }
      if ( adjoint && print_level > 0 )
      {
         if (_braid_CoreElt(core, max_levels) > 1)
         {
            _braid_printf("\n");
            _braid_printf("  Braid:      || r ||      || r_adj ||     Objective\n");
            _braid_printf("  Braid:---------------------------------------------\n");
         }
         else
         {
            _braid_printf("  Braid: Serial time-stepping. \n\n");
         }
      }
   }

   /* Set chunk size */
   _braid_CoreElt(core, ntime)    = (int) (ntime0 / nchunks);
   _braid_CoreElt(core, gupper)   = _braid_CoreElt(core, ntime);
   dt_chunk = (tstop0 - tstart0 ) / nchunks;

   /* Allocate and initialize grids */
   if ( !warm_restart )
   {
      /* Create fine grid */
      _braid_GetDistribution(core, &ilower, &iupper);
      _braid_GridInit(core, 0, ilower, iupper, &grid);

      /* Set t values */
      ta = _braid_GridElt(grid, ta);
      if ( _braid_CoreElt(core, tgrid) != NULL )
      {
         /* Call the user's time grid routine */
         _braid_BaseTimeGrid(core, app, ta, &ilower, &iupper);
      }
      else
      {
         for (i = ilower; i <= iupper; i++)
         {
            ta[i-ilower] = tstart0 + (((braid_Real)i)/_braid_CoreElt(core,ntime))*(dt_chunk);
         }
      }

      /* Create a grid hierarchy */
      _braid_InitHierarchy(core, grid, 0);

      /* Set initial values */
      _braid_InitGuess(core, 0);
   }

   /* Initialize sensitivity computation */
   if ( adjoint )
   {
      if (!warm_restart)
      {
         /* Initialize and allocate the adjoint variables */
         _braid_InitAdjointVars(core, grid);
      }
      else
      {
         /* Prepare for next adjoint iteration in case of warm_restart */
         _braid_CoreElt(core, optim)->sum_user_obj  = 0.0;
         _braid_CoreElt(core, optim)->f_bar         = 0.0;
         if (!obj_only)
         {
            _braid_CoreFcn(core, reset_gradient)(_braid_CoreElt(core, app));
         }
      }

      if ( obj_only )
      {
         _braid_CoreElt(core, record) = 0;
      }
      else
      {
         _braid_CoreElt(core, record) = 1;
      }
   }

   /* Turn on warm_restart, so that further calls to braid_drive() don't
    * initialize the grid and adjoint again. */
   _braid_CoreElt(core, warm_restart) = 1;

   /* Start timer */
   localtime = MPI_Wtime();

   /* Loop over all time chunks */
   for (ichunk = 0; ichunk < nchunks; ichunk++)
   {
      _braid_CoreElt(core, ichunk) = ichunk;

      /* Set start and end time values of current time chunk */
      _braid_CoreElt(core,tstart) = tstart0 + ichunk * dt_chunk;
      _braid_CoreElt(core,tstop)  = _braid_CoreElt(core, tstart) + dt_chunk;

      /* Output */
      if (myid == 0)
      {
         _braid_printf("\n  Braid Chunk %d: [%f, %f], %d time steps\n",
                       ichunk, _braid_CoreElt(core, tstart), _braid_CoreElt(core, tstop),
                       _braid_CoreElt(core,ntime));
      }

      /* Initialize the chunk */
      if ( ichunk > 0 )
      {
         /* Set new initial condition */
         _braid_ChunkSetInitialCondition(core);

         /* Set new time vector ta on all levels */
         for (level = 0; level < _braid_CoreElt(core, nlevels); level++)
         {
            ta = _braid_GridElt(grids[level], ta);
            ilower = _braid_GridElt(grids[level], ilower);
            iupper = _braid_GridElt(grids[level], iupper);
            for (i = ilower-1; i <= iupper+1; i++)
            {
               ta[i-ilower] += dt_chunk ;
            }
         }
      }

      /* Solve this time chunk */
      if ( _braid_CoreElt(core, trimgrit) )
      {
         /* TriMGRIT */
         _braid_TriDrive(core, localtime);
      }
      else
      {
         /* MGRIT */
         _braid_Drive(core, localtime);
      }

      /* Get the time */
      braid_Real mytimediff = MPI_Wtime() - localtime;
      MPI_Allreduce(&mytimediff, &globaltime, 1, braid_MPI_REAL, MPI_MAX, comm_world);
      _braid_CoreElt(core, localtime)  = mytimediff;
      _braid_CoreElt(core, globaltime) = globaltime;

      /* Print statistics for this run */
      if ( (print_level > 1) && (myid == 0) )
      {
         braid_PrintStats(core);
      }
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_Init(MPI_Comm               comm_world,
           MPI_Comm               comm,
           braid_Real             tstart,
           braid_Real             tstop,
           braid_Int              ntime,
           braid_App              app,
           braid_PtFcnStep        step,
           braid_PtFcnInit        init,
           braid_PtFcnClone       clone,
           braid_PtFcnFree        free,
           braid_PtFcnSum         sum,
           braid_PtFcnSpatialNorm spatialnorm,
           braid_PtFcnAccess      access,
           braid_PtFcnBufSize     bufsize,
           braid_PtFcnBufPack     bufpack,
           braid_PtFcnBufUnpack   bufunpack,
           braid_Core            *core_ptr)
{
   _braid_Core           *core;

   /* Braid default values */
   braid_Int              nchunks         = 1;              /* Default: all computation in one chunk */
   braid_Int              cfdefault       = 2;              /* Default coarsening factor */
   braid_Int              nrdefault       = 1;              /* Default number of FC sweeps on each level */
   braid_Int              fmg             = 0;              /* Default fmg (0 is off) */
   braid_Int              nfmg            = -1;             /* Default fmg cycles is -1, indicating all fmg-cycles (if fmg=1) */
   braid_Int              nfmg_Vcyc       = 1;              /* Default num V-cycles at each fmg level is 1 */
   braid_Int              max_iter        = 100;            /* Default max_iter */
   braid_Int              max_levels      = 30;             /* Default max_levels */
   braid_Int              min_coarse      = 2;              /* Default min_coarse */
   braid_Int              seq_soln        = 0;              /* Default initial guess is from user's Init() function */
   braid_Int              print_level     = 2;              /* Default print level */
   braid_Int              io_level        = 1;              /* Default output-to-file level */
   braid_Int              access_level    = 1;              /* Default access level */
   braid_Int              tnorm           = 2;              /* Default temporal norm */
   braid_Real             tol             = 1.0e-09;        /* Default absolute tolerance */
   braid_Int              warm_restart    = 0;              /* Default is no warm restart */
   braid_Int              rtol            = 1;              /* Use relative tolerance */
   braid_Int              skip            = 1;              /* Default skip value, skips all work on first down-cycle */
   braid_Int              max_refinements = 200;            /* Maximum number of F-refinements */
   braid_Int              tpoints_cutoff  = braid_Int_Max;  /* Maximum number of time steps, controls FRefine()*/
   braid_Int              adjoint         = 0;              /* Default adjoint run: Turned off */
   braid_Int              record          = 0;              /* Default action recording: Turned off */
   braid_Int              obj_only        = 0;              /* Default objective only: Turned off */
   braid_Int              verbose_adj     = 0;              /* Default adjoint verbosity Turned off */

   braid_Int              myid_world,  myid;

   MPI_Comm_rank(comm_world, &myid_world);
   MPI_Comm_rank(comm, &myid);

   core = _braid_CTAlloc(_braid_Core, 1);

   _braid_CoreElt(core, comm_world)      = comm_world;
   _braid_CoreElt(core, comm)            = comm;
   _braid_CoreElt(core, myid_world)      = myid_world;
   _braid_CoreElt(core, myid)            = myid;
   _braid_CoreElt(core, tstart)          = tstart;
   _braid_CoreElt(core, tstop)           = tstop;
   _braid_CoreElt(core, ntime)           = ntime;
   _braid_CoreElt(core, nchunks)         = nchunks;
   _braid_CoreElt(core, ichunk)          = 0;
   _braid_CoreElt(core, dt_chunk)        = tstop - tstart;
   _braid_CoreElt(core, app)             = app;

   _braid_CoreElt(core, step)            = step;
   _braid_CoreElt(core, init)            = init;
   _braid_CoreElt(core, sinit)           = NULL;
   _braid_CoreElt(core, clone)           = clone;
   _braid_CoreElt(core, sclone)          = NULL;
   _braid_CoreElt(core, free)            = free;
   _braid_CoreElt(core, sfree)           = NULL;
   _braid_CoreElt(core, sum)             = sum;
   _braid_CoreElt(core, spatialnorm)     = spatialnorm;
   _braid_CoreElt(core, access)          = access;
   _braid_CoreElt(core, bufsize)         = bufsize;
   _braid_CoreElt(core, bufpack)         = bufpack;
   _braid_CoreElt(core, bufunpack)       = bufunpack;
   _braid_CoreElt(core, residual)        = NULL;
   _braid_CoreElt(core, scoarsen)        = NULL;
   _braid_CoreElt(core, srefine)         = NULL;
   _braid_CoreElt(core, tgrid)           = NULL;

   _braid_CoreElt(core, access_level)    = access_level;
   _braid_CoreElt(core, tnorm)           = tnorm;
   _braid_CoreElt(core, print_level)     = print_level;
   _braid_CoreElt(core, io_level)        = io_level;
   _braid_CoreElt(core, max_levels)      = 0; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, min_coarse)      = min_coarse;
   _braid_CoreElt(core, seq_soln)        = seq_soln;
   _braid_CoreElt(core, tol)             = tol;
   _braid_CoreElt(core, rtol)            = rtol;
   _braid_CoreElt(core, warm_restart)    = warm_restart;

   _braid_CoreElt(core, nrels)           = NULL; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, nrdefault)       = nrdefault;

   _braid_CoreElt(core, cfactors)        = NULL; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, cfdefault)       = cfdefault;

   _braid_CoreElt(core, max_iter)        = 0; /* Set with SetMaxIter() below */
   _braid_CoreElt(core, niter)           = 0;
   _braid_CoreElt(core, fmg)             = fmg;
   _braid_CoreElt(core, nfmg)            = nfmg;
   _braid_CoreElt(core, nfmg_Vcyc)       = nfmg_Vcyc;

   _braid_CoreElt(core, storage)         = -1;            /* only store C-points */
   _braid_CoreElt(core, useshell)         = 0;

   _braid_CoreElt(core, gupper)          = ntime;

   _braid_CoreElt(core, refine)          = 0;  /* Time refinement off by default */
   _braid_CoreElt(core, rfactors)        = NULL;
   _braid_CoreElt(core, r_space)         = 0;
   _braid_CoreElt(core, rstopped)        = -1;
   _braid_CoreElt(core, nrefine)         = 0;
   _braid_CoreElt(core, max_refinements) = max_refinements;
   _braid_CoreElt(core, tpoints_cutoff)  = tpoints_cutoff;

   _braid_CoreElt(core, nlevels)         = 0;
   _braid_CoreElt(core, grids)           = NULL; /* Set with SetMaxLevels() below */

   _braid_CoreElt(core, skip)            = skip;

   _braid_CoreElt(core, adjoint)               = adjoint;
   _braid_CoreElt(core, record)                = record;
   _braid_CoreElt(core, obj_only)              = obj_only;
   _braid_CoreElt(core, verbose_adj)           = verbose_adj;
   _braid_CoreElt(core, actionTape)            = NULL;
   _braid_CoreElt(core, userVectorTape)        = NULL;
   _braid_CoreElt(core, barTape)               = NULL;
   _braid_CoreElt(core, optim)                 = NULL;
   _braid_CoreElt(core, objectiveT)            = NULL;
   _braid_CoreElt(core, objT_diff)             = NULL;
   _braid_CoreElt(core, step_diff)             = NULL;
   _braid_CoreElt(core, reset_gradient)        = NULL;
   _braid_CoreElt(core, postprocess_obj)       = NULL;
   _braid_CoreElt(core, postprocess_obj_diff)  = NULL;

   /* Residual history and accuracy tracking for StepStatus*/
   _braid_CoreElt(core, rnorm0)              = braid_INVALID_RNORM;
   _braid_CoreElt(core, rnorms)              = NULL; /* Set with SetMaxIter() below */
   _braid_CoreElt(core, full_rnorm_res)      = NULL;
   _braid_CoreElt(core, full_rnorm0)         = braid_INVALID_RNORM;
   _braid_CoreElt(core, full_rnorms)         = NULL; /* Set with SetMaxIter() below */
   _braid_CoreElt(core, old_fine_tolx)       = -1.0;
   _braid_CoreElt(core, tight_fine_tolx)     = 1;

   braid_SetMaxLevels(core, max_levels);
   braid_SetMaxIter(core, max_iter);

   *core_ptr = core;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_InitAdjoint(braid_PtFcnObjectiveT        objectiveT,
                  braid_PtFcnObjectiveTDiff    objT_diff,
                  braid_PtFcnStepDiff          step_diff,
                  braid_PtFcnResetGradient     reset_gradient,
                  braid_Core                  *core_ptr)
{
   braid_Optim          optim;

   /* Set adjoint flags */
   _braid_CoreElt(*core_ptr, adjoint)      = 1;
   _braid_CoreElt(*core_ptr, record)       = 1;
   _braid_CoreElt(*core_ptr, skip)         = 0;

   /* Define default values */
   braid_Real  tstart_obj     = _braid_CoreElt(*core_ptr, tstart);
   braid_Real  tstop_obj      = _braid_CoreElt(*core_ptr, tstop);
   braid_Real  tol_adj        = 1e-6;
   braid_Int   rtol_adj       = 1;

   /* Allocate memory for the optimization structure */
   optim = (struct _braid_Optimization_struct*) malloc(sizeof(struct _braid_Optimization_struct));

   /* Set optimization variables */
   optim->adjoints       = NULL;    /* will be allocated in InitAdjointVars() */
   optim->tapeinput      = NULL;    /* will be allocated in InitAdjointVars() */
   optim->sendbuffer     = NULL;    /* will be allocated in InitAdjointVars() */
   optim->request        = NULL;    /* will be allocated in InitAdjointVars() */
   optim->objective      = 0.0;
   optim->sum_user_obj   = 0.0;
   optim->f_bar          = 0.0;
   optim->tstart_obj     = tstart_obj;
   optim->tstop_obj      = tstop_obj;
   optim->tol_adj        = tol_adj;
   optim->rtol_adj       = rtol_adj;
   optim->rnorm_adj      = braid_INVALID_RNORM;
   optim->rnorm0_adj     = braid_INVALID_RNORM;
   optim->rnorm          = braid_INVALID_RNORM;
   optim->rnorm0         = braid_INVALID_RNORM;

   /* Store the optim structure in the core */
   _braid_CoreElt(*core_ptr, optim) = optim;

   /* Initialize the tapes */
   _braid_TapeInit( _braid_CoreElt(*core_ptr, actionTape) );
   _braid_TapeInit( _braid_CoreElt(*core_ptr, userVectorTape) );
   _braid_TapeInit( _braid_CoreElt(*core_ptr, barTape) );

   /* Set the user functions */
   _braid_CoreElt(*core_ptr, objectiveT)     = objectiveT;
   _braid_CoreElt(*core_ptr, step_diff)      = step_diff;
   _braid_CoreElt(*core_ptr, objT_diff)      = objT_diff;
   _braid_CoreElt(*core_ptr, reset_gradient) = reset_gradient;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_InitTriMGRIT(MPI_Comm               comm_world,
                   MPI_Comm               comm,
                   braid_Real             tstart,
                   braid_Real             tstop,
                   braid_Int              ntime,
                   braid_App              app,
                   braid_PtFcnTriResidual triresidual,
                   braid_PtFcnTriSolve    trisolve,
                   braid_PtFcnInit        init,
                   braid_PtFcnClone       clone,
                   braid_PtFcnFree        free,
                   braid_PtFcnSum         sum,
                   braid_PtFcnSpatialNorm spatialnorm,
                   braid_PtFcnAccess      access,
                   braid_PtFcnBufSize     bufsize,
                   braid_PtFcnBufPack     bufpack,
                   braid_PtFcnBufUnpack   bufunpack,
                   braid_Core            *core_ptr)
{
   _braid_Core      *core;

   braid_PtFcnStep   step = NULL;

   /* Use normal braid init function (for MGRIT), then modify for TriMGRIT */
   braid_Init(comm_world, comm, tstart, tstop, ntime,
              app, step, init, clone, free, sum, spatialnorm, access,
              bufsize, bufpack, bufunpack, &core);

   /* These are all set to 0 or NULL in the CTAlloc of core in braid_Init() */
   _braid_CoreElt(core, trimgrit)    = 1;
   _braid_CoreElt(core, triresidual) = triresidual;
   _braid_CoreElt(core, trisolve)    = trisolve;

   /* These are the only values currently supported for TriMGRIT */
   _braid_CoreElt(core, nchunks)  = 1;
   _braid_CoreElt(core, seq_soln) = 0;
   _braid_CoreElt(core, fmg)      = 0;
   _braid_CoreElt(core, storage)  = 0;  /* Store all points on all levels */
   _braid_CoreElt(core, useshell) = 0;
   _braid_CoreElt(core, refine)   = 0;  /* Time refinement off (for now) */
   _braid_CoreElt(core, r_space)  = 0;  /* No refinement (for now) */
   _braid_CoreElt(core, skip)     = 0;

   *core_ptr = core;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_Destroy(braid_Core  core)
{
   if (core)
   {
      braid_App               app        = _braid_CoreElt(core, app);
      braid_Int               nlevels    = _braid_CoreElt(core, nlevels);
      _braid_Grid           **grids      = _braid_CoreElt(core, grids);
      braid_Int               level;

      _braid_TFree(_braid_CoreElt(core, nrels));
      _braid_TFree(_braid_CoreElt(core, rnorms));
      _braid_TFree(_braid_CoreElt(core, full_rnorms));
      _braid_TFree(_braid_CoreElt(core, cfactors));
      _braid_TFree(_braid_CoreElt(core, rfactors));
      _braid_TFree(_braid_CoreElt(core, tnorm_a));

      /* Destroy the optimization structure */
      _braid_CoreElt(core, record) = 0;
      if (_braid_CoreElt(core, adjoint))
      {
         _braid_OptimDestroy( core );
         _braid_TFree(_braid_CoreElt(core, optim));
      }

      /* Free last time step, if set */
      if (_braid_CoreElt(core, storage) < 0 )
      {
         if (_braid_GridElt(grids[0], ulast) != NULL)
         {
            _braid_BaseFree(core, app, _braid_GridElt(grids[0], ulast));
         }
      }

      for (level = 0; level < nlevels; level++)
      {
         _braid_GridDestroy(core, grids[level]);
      }

      _braid_TFree(grids);

      _braid_TFree(core);
   }

   if (_braid_printfile != NULL)
   {
      fclose(_braid_printfile);
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_PrintStats(braid_Core  core)
{
   braid_Int     myid          = _braid_CoreElt(core, myid_world);
   braid_Real    tstart        = _braid_CoreElt(core, tstart);
   braid_Real    tstop         = _braid_CoreElt(core, tstop);
   braid_Int     gupper        = _braid_CoreElt(core, gupper);
   braid_Int     nchunks       = _braid_CoreElt(core, nchunks);
   braid_Int     max_levels    = _braid_CoreElt(core, max_levels);
   braid_Int     min_coarse    = _braid_CoreElt(core, min_coarse);
   braid_Int     seq_soln      = _braid_CoreElt(core, seq_soln);
   braid_Int     storage       = _braid_CoreElt(core, storage);
   braid_Real    tol           = _braid_CoreElt(core, tol);
   braid_Int     rtol          = _braid_CoreElt(core, rtol);
   braid_Int    *nrels         = _braid_CoreElt(core, nrels);
   /*braid_Int    *cfactors     = _braid_CoreElt(core, cfactors);*/
   braid_Int     max_iter      = _braid_CoreElt(core, max_iter);
   braid_Int     nrefine       = _braid_CoreElt(core, nrefine);
   braid_Int     niter         = _braid_CoreElt(core, niter);
   braid_Int     nlevels       = _braid_CoreElt(core, nlevels);
   braid_Int     tnorm         = _braid_CoreElt(core, tnorm);
   braid_Int     fmg           = _braid_CoreElt(core, fmg);
   braid_Int     nfmg          = _braid_CoreElt(core, nfmg);
   braid_Int     nfmg_Vcyc     = _braid_CoreElt(core, nfmg_Vcyc);
   braid_Int     access_level  = _braid_CoreElt(core, access_level);
   braid_Int     print_level   = _braid_CoreElt(core, print_level);
   braid_Int     skip          = _braid_CoreElt(core, skip);
   braid_Real    globaltime    = _braid_CoreElt(core, globaltime);
   braid_PtFcnResidual fullres = _braid_CoreElt(core, full_rnorm_res);
   _braid_Grid **grids         = _braid_CoreElt(core, grids);
   braid_Int     adjoint       = _braid_CoreElt(core, adjoint);
   braid_Optim   optim         = _braid_CoreElt(core, optim);

   braid_Real    tol_adj;
   braid_Int     rtol_adj;
   braid_Real    rnorm, rnorm_adj;
   braid_Int     level;

   if (adjoint)
   {
      tol_adj   = optim->tol_adj;
      rtol_adj  = optim->rtol_adj;
      rnorm_adj = optim->rnorm_adj;
   }

   _braid_GetRNorm(core, -1, &rnorm);

   if ( myid == 0 )
   {
      _braid_printf("\n");
      _braid_printf("  Braid Solver Stats:\n");
      _braid_printf("  start time = %e\n", tstart);
      _braid_printf("  stop time  = %e\n", tstop);
      _braid_printf("  time steps = %d\n", gupper);
      _braid_printf("  time chunks = %d\n", nchunks);
      _braid_printf("\n");
      _braid_printf("  use seq soln?         = %d\n", seq_soln);
      _braid_printf("  storage               = %d\n", storage);
      _braid_printf("\n");

      _braid_printf("  max iterations        = %d\n", max_iter);
      _braid_printf("  iterations            = %d\n", niter);
      _braid_printf("\n");

      if ( adjoint )
      {
         _braid_printf("  state   residual norm =  %e", rnorm);
         if ( rtol ) _braid_printf("  (-> rel. stopping tol. = %1.2e)\n", tol);
         else        _braid_printf("  (-> abs. stopping tol. = %1.2e)\n", tol);
         _braid_printf("  adjoint residual norm =  %e", rnorm_adj);
         if (rtol_adj ) _braid_printf("  (-> rel. stopping tol. = %1.2e)\n", tol_adj);
         else           _braid_printf("  (-> abs. stopping tol. = %1.2e)\n", tol_adj);
      }
      else
      {
         _braid_printf("  residual norm         = %e\n", rnorm);
         _braid_printf("  stopping tolerance    = %e\n", tol);
         _braid_printf("  use relative tol?     = %d\n", rtol);
      }

      if (tnorm == 1)
      {
         _braid_printf("                                          --> 1-norm TemporalNorm \n");
      }
      else if (tnorm == 2)
      {
         _braid_printf("                                          --> 2-norm TemporalNorm \n");
      }
      else if (tnorm == 3)
      {
         _braid_printf("                                          --> Inf-norm TemporalNorm \n");
      }
      if (fullres != NULL)
      {
         _braid_GetFullRNorm(core, -1, &rnorm);
         _braid_printf("  Global res 2-norm     = %e\n", rnorm);
      }

      _braid_printf("\n");
      _braid_printf("  use fmg?              = %d\n", fmg);
      if ( fmg )
      {
         _braid_printf("  V-cycles / FMG level  = %d\n", nfmg_Vcyc);
         if ( nfmg != -1 )
         {
            _braid_printf("  number fmg cycles     = %d\n", nfmg);
         }
         else
         {
            _braid_printf("  fmg-cycles for all iteratons\n");
         }
      }
      _braid_printf("  access_level          = %d\n", access_level);
      _braid_printf("  print_level           = %d\n\n", print_level);
      _braid_printf("  max number of levels  = %d\n", max_levels);
      _braid_printf("  min coarse            = %d\n", min_coarse);
      _braid_printf("  number of levels      = %d\n", nlevels);
      _braid_printf("  skip down cycle       = %d\n", skip);
      _braid_printf("  number of refinements = %d\n", nrefine);
      _braid_printf("\n");
      _braid_printf("  level   time-pts   cfactor   nrelax\n");
      for (level = 0; level < nlevels-1; level++)
      {
         _braid_printf("  % 5d  % 8d  % 7d   % 6d\n",
                       level, _braid_GridElt(grids[level], gupper),
                       _braid_GridElt(grids[level], cfactor), nrels[level]);
      }
      /* Print out coarsest level information */
      _braid_printf("  % 5d  % 8d  \n",
                    level, _braid_GridElt(grids[level], gupper) );
      _braid_printf("\n");
      _braid_printf("  wall time = %f\n", globaltime);
      _braid_printf("\n");
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetMaxLevels(braid_Core  core,
                   braid_Int   max_levels)
{
   braid_Int              old_max_levels = _braid_CoreElt(core, max_levels);
   braid_Int             *nrels          = _braid_CoreElt(core, nrels);
   braid_Int             *cfactors       = _braid_CoreElt(core, cfactors);
   _braid_Grid          **grids          = _braid_CoreElt(core, grids);
   braid_Int              level;

   _braid_CoreElt(core, max_levels) = max_levels;

   nrels = _braid_TReAlloc(nrels, braid_Int, max_levels);
   cfactors = _braid_TReAlloc(cfactors, braid_Int, max_levels);
   grids    = _braid_TReAlloc(grids, _braid_Grid *, max_levels);
   for (level = old_max_levels; level < max_levels; level++)
   {
      nrels[level]    = -1;
      cfactors[level] = 0;
      grids[level]    = NULL;
   }
   _braid_CoreElt(core, nrels)    = nrels;
   _braid_CoreElt(core, cfactors) = cfactors;
   _braid_CoreElt(core, grids)    = grids;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetSkip(braid_Core  core,
              braid_Int   skip)
{
   /* Do not set skip=1 if we do sequential integration first */
   if (_braid_CoreElt(core, seq_soln) == 0)
      _braid_CoreElt(core, skip) = skip;
   else if (skip == 1)
      _braid_printf("  Braid: The skip option is not compatible with SeqSoln option\n");

   return _braid_error_flag;
}


/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetMinCoarse(braid_Core  core,
                   braid_Int   min_coarse)
{
   _braid_CoreElt(core, min_coarse) = min_coarse;

   return _braid_error_flag;
}


/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetPrintLevel(braid_Core  core,
                    braid_Int   print_level)
{
   _braid_CoreElt(core, print_level) = print_level;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetFileIOLevel(braid_Core  core,
                     braid_Int   io_level)
{
   _braid_CoreElt(core, io_level) = io_level;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetPrintFile(braid_Core     core,
                   const char    *printfile_name)
{
   braid_Int  myid = _braid_CoreElt(core, myid_world);

   if (myid == 0)
   {
      if ((_braid_printfile = fopen(printfile_name, "w")) == NULL)
      {
         printf("  Braid: Error: can't open output file %s\n", printfile_name);
         exit(1);
      }
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetDefaultPrintFile(braid_Core     core)
{
   const char fname[] = "braid_runtime.out";
   braid_SetPrintFile(core, fname);
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetAccessLevel(braid_Core  core,
                     braid_Int   access_level)
{
   _braid_CoreElt(core, access_level) = access_level;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SplitCommworld(const MPI_Comm  *comm_world,
                     braid_Int       px,
                     MPI_Comm        *comm_x,
                     MPI_Comm        *comm_t)
{
   braid_Int myid, xcolor, tcolor;

   /* Create communicators for the time and space dimensions */
   /* The communicators are based on colors and keys (= myid) */
   MPI_Comm_rank( *comm_world, &myid );
   xcolor = myid / px;
   tcolor = myid % px;

   MPI_Comm_split( *comm_world, xcolor, myid, comm_x );
   MPI_Comm_split( *comm_world, tcolor, myid, comm_t );

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetAbsTol(braid_Core  core,
                braid_Real  tol)
{
   _braid_CoreElt(core, tol)  = tol;
   _braid_CoreElt(core, rtol) = 0;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetRelTol(braid_Core  core,
                braid_Real  tol)
{
   _braid_CoreElt(core, tol)  = tol;
   _braid_CoreElt(core, rtol) = 1;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetNRelax(braid_Core  core,
                braid_Int   level,
                braid_Int   nrelax)
{
   braid_Int  *nrels = _braid_CoreElt(core, nrels);

   if (level < 0)
   {
      /* Set default value */
      _braid_CoreElt(core, nrdefault) = nrelax;
   }
   else
   {
      /* Set factor on specified level */
      nrels[level] = nrelax;
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetCFactor(braid_Core  core,
                 braid_Int   level,
                 braid_Int   cfactor)
{
   braid_Int  *cfactors = _braid_CoreElt(core, cfactors);

   if (level < 0)
   {
      /* Set default value */
      _braid_CoreElt(core, cfdefault) = cfactor;
   }
   else
   {
      /* Set factor on specified level */
      cfactors[level] = cfactor;
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetMaxIter(braid_Core  core,
                 braid_Int   max_iter)
{
   braid_Real  *rnorms      = _braid_CoreElt(core, rnorms);
   braid_Real  *full_rnorms = _braid_CoreElt(core, full_rnorms);
   braid_Int    next_iter   = _braid_CoreElt(core, niter) + 1;
   braid_Int    i;

   /* If rnorms has never been allocated, make sure all entries are initialized */
   if (rnorms == NULL)
   {
      next_iter = 0;
   }

   _braid_CoreElt(core, max_iter) = max_iter;

   rnorms = _braid_TReAlloc(rnorms, braid_Real, max_iter+1);
   for (i = next_iter; i <= max_iter; i++)
   {
      rnorms[i] = braid_INVALID_RNORM;
   }

   /* Allocate even if not using full rnorms (simplifies intialization) */
   full_rnorms = _braid_TReAlloc(full_rnorms, braid_Real, max_iter+1);
   for (i = next_iter; i <= max_iter; i++)
   {
      full_rnorms[i] = braid_INVALID_RNORM;
   }

   _braid_CoreElt(core, rnorms) = rnorms;
   _braid_CoreElt(core, full_rnorms) = full_rnorms;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetRefine(braid_Core  core,
                braid_Int   refine)
{
   _braid_CoreElt(core, refine) = refine;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetMaxRefinements(braid_Core  core,
                        braid_Int   max_refinements)
{
   _braid_CoreElt(core, max_refinements) = max_refinements;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetTPointsCutoff(braid_Core  core,
                       braid_Int   tpoints_cutoff)
{
   _braid_CoreElt(core, tpoints_cutoff) = tpoints_cutoff;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetFMG(braid_Core  core)
{
   _braid_CoreElt(core, fmg) = 1;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetNFMGVcyc(braid_Core  core,
                  braid_Int   nfmg_Vcyc)
{
   _braid_CoreElt(core, nfmg_Vcyc) = nfmg_Vcyc;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetNFMG(braid_Core  core,
              braid_Int   k)
{
   _braid_CoreElt(core, nfmg) = k;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetStorage(braid_Core  core,
                 braid_Int   storage)
{
   _braid_CoreElt(core, storage) = storage;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetTemporalNorm(braid_Core  core,
                      braid_Int   tnorm)
{
   _braid_CoreElt(core, tnorm) = tnorm;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetResidual(braid_Core          core,
                  braid_PtFcnResidual residual)
{
   _braid_CoreElt(core, residual) = residual;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetFullRNormRes(braid_Core          core,
                      braid_PtFcnResidual residual)
{
   _braid_CoreElt(core, full_rnorm_res) = residual;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetTimeGrid(braid_Core          core,
                  braid_PtFcnTimeGrid tgrid
   )
{
   _braid_CoreElt(core, tgrid) = tgrid;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetSpatialCoarsen(braid_Core          core,
                        braid_PtFcnSCoarsen scoarsen)
{
   _braid_CoreElt(core, scoarsen) = scoarsen;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetSpatialRefine(braid_Core         core,
                       braid_PtFcnSRefine srefine)
{
   _braid_CoreElt(core, srefine) = srefine;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetShell(braid_Core          core,
               braid_PtFcnSInit    sinit,
               braid_PtFcnSClone   sclone,
               braid_PtFcnSFree    sfree)
{
   _braid_CoreElt(core, sinit) = sinit;
   _braid_CoreElt(core, sclone) = sclone;
   _braid_CoreElt(core, sfree) = sfree;
   _braid_CoreElt(core, useshell) = 1;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetNumIter(braid_Core   core,
                 braid_Int   *niter_ptr)
{
   *niter_ptr =  _braid_CoreElt(core, niter);
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetRNorms(braid_Core   core,
                braid_Int   *nrequest_ptr,
                braid_Real  *rnorms)
{
   braid_Real  *rnorms_all = _braid_CoreElt(core, rnorms);
   braid_Int    rnorms_len = _braid_CoreElt(core, niter) + 1;

   _braid_GetNEntries(rnorms_all, rnorms_len, nrequest_ptr, rnorms);
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetNLevels(braid_Core  core,
                 braid_Int  *nlevels_ptr)

{
   *nlevels_ptr = _braid_CoreElt(core, nlevels);
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetSpatialAccuracy( braid_StepStatus  status,
                          braid_Real        loose_tol,
                          braid_Real        tight_tol,
                          braid_Real       *tol_ptr )
{
   braid_Int nrequest   = 2;
   braid_Real stol, tol, rnorm, rnorm0, old_fine_tolx;
   braid_Int level;
   braid_Real l_rnorm, l_ltol, l_ttol, l_tol;
   braid_Real *rnorms = (braid_Real *) malloc( 2*sizeof(braid_Real) );

   braid_StepStatusGetTol(status, &tol);
   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetOldFineTolx(status, &old_fine_tolx);

   /* Get the first and then the current residual norms */
   rnorms[0] = -1.0; rnorms[1] = -1.0;
   braid_StepStatusGetRNorms(status, &nrequest, rnorms);
   if((rnorms[0] == -1.0) && (rnorms[1] != -1.0)){
      rnorm0 = rnorms[1];
   }
   else{
      rnorm0 = rnorms[0];
   }
   nrequest = -2;
   braid_StepStatusGetRNorms(status, &nrequest, rnorms);
   if((rnorms[1] == -1.0) && (rnorms[0] != -1.0)){
      rnorm = rnorms[0];
   }
   else{
      rnorm = rnorms[1];
   }


   if ( (level > 0) || (nrequest == 0) || (rnorm0 == -1.0) )
   {
      /* Always return the loose tolerance, if
       * (1) On a coarse grid computation
       * (2) There is no residual history yet (this is the first Braid iteration with skip turned on) */
      *tol_ptr = loose_tol;
   }
   else
   {
      /* Else, do a variable tolerance for the fine grid */
      l_rnorm = -log10(rnorm / rnorm0);
      l_tol   = -log10(tol / rnorm0);
      l_ltol  = -log10(loose_tol);
      l_ttol  = -log10(tight_tol);

      if ( l_rnorm >= (7.0/8.0)*l_tol )
      {
         /* Close to convergence, return tight_tol */
         *tol_ptr = tight_tol;
      }
      else
      {
         /* Do linear interpolation between loose_tol and tight_tol (but with respect to log10) */
         stol = (l_rnorm / l_tol) * (l_ttol - l_ltol) + l_ltol;
         *tol_ptr = pow(10, -stol);

         /* The fine grid tolerance MUST never decrease */
         if ( ((*tol_ptr) > old_fine_tolx) && (old_fine_tolx > 0) )
         {
            *tol_ptr = old_fine_tolx;
         }
      }
   }

   if (level == 0)
   {
      /* Store this fine grid tolerance */
      braid_StepStatusSetOldFineTolx(status, (*tol_ptr));

      /* If we've reached the "tight tolerance", then indicate to Braid that we can halt */
      if ( *tol_ptr == tight_tol )
      {
         braid_StepStatusSetTightFineTolx(status, 1);
      }
      else
      {
         braid_StepStatusSetTightFineTolx(status, 0);
      }
   }

   free(rnorms);
   /* printf( "lev: %d, accuracy: %1.2e, nreq: %d, rnorm: %1.2e, rnorm0: %1.2e, loose: %1.2e, tight: %1.2e, old: %1.2e, braid_tol: %1.2e \n", level, *tol_ptr, nrequest, rnorm, rnorm0, loose_tol, tight_tol, old_fine_tolx, tol); */
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetSeqSoln(braid_Core  core,
                 braid_Int   seq_soln)
{
   /* Skip needs to be 0 if we do a sequential integration first */
   _braid_CoreElt(core, seq_soln) = seq_soln;
   if (seq_soln == 1)
      _braid_CoreElt(core, skip) = 0;

   return _braid_error_flag;
}

/**----------------------------------------------------------------------------
 * Adjoint
 *-----------------------------------------------------------------------------*/

braid_Int
braid_SetTStartObjective(braid_Core core,
                         braid_Real tstart_obj)
{
   if ( !(_braid_CoreElt(core, adjoint)))
   {
      return _braid_error_flag;
   }

   _braid_CoreElt(core, optim->tstart_obj) = tstart_obj;

   /* Sanity check */
   if ( tstart_obj < _braid_CoreElt(core, tstart) )
   {
      _braid_printf("\n  Braid: WARNING: tstart_objective < tstart ! Using default tstart now.\n\n");
      _braid_CoreElt(core, optim->tstart_obj) = _braid_CoreElt(core, tstart);
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetTStopObjective(braid_Core core,
                        braid_Real tstop_obj)
{
   if ( !(_braid_CoreElt(core, adjoint)) )
   {
      return _braid_error_flag;
   }

   _braid_CoreElt(core, optim->tstop_obj) = tstop_obj;

   /* Sanity check */
   if ( tstop_obj > _braid_CoreElt(core, tstop) )
   {
      _braid_printf("\n  Braid: WARNING: tstop_objective > tstop ! Using default tstop now.\n\n");
      _braid_CoreElt(core, optim->tstop_obj) = _braid_CoreElt(core, tstop);
   }


   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetPostprocessObjective(braid_Core                      core,
                              braid_PtFcnPostprocessObjective post_fcn )
{
   if ( !(_braid_CoreElt(core, adjoint)) )
   {
      return _braid_error_flag;
   }

   _braid_CoreElt(core, postprocess_obj) = post_fcn;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetPostprocessObjective_diff(braid_Core                           core,
                                   braid_PtFcnPostprocessObjective_diff post_fcn_diff )
{
   if ( !(_braid_CoreElt(core, adjoint)) )
   {
      return _braid_error_flag;
   }

   _braid_CoreElt(core, postprocess_obj_diff) = post_fcn_diff;
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetAbsTolAdjoint(braid_Core core,
                       braid_Real tol_adj)
{
   if ( !(_braid_CoreElt(core, adjoint)) )
   {
      return _braid_error_flag;
   }

   _braid_CoreElt(core, optim)->tol_adj  = tol_adj;
   _braid_CoreElt(core, optim)->rtol_adj = 0;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetRelTolAdjoint(braid_Core core,
                       braid_Real tol_adj)
{
   if ( !(_braid_CoreElt(core, adjoint)) )
   {
      return _braid_error_flag;
   }

   _braid_CoreElt(core, optim)->tol_adj  = tol_adj;
   _braid_CoreElt(core, optim)->rtol_adj = 1;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetObjective(braid_Core  core,
                   braid_Real *objective_ptr )
{
   if ( !(_braid_CoreElt(core, adjoint)) )
   {
      *objective_ptr = 0.0;
   }
   else
   {
      *objective_ptr = _braid_CoreElt(core, optim)->objective;
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetObjectiveOnly(braid_Core core,
                       braid_Int  boolean)
{
   if ( !(_braid_CoreElt(core, adjoint)) )
   {
      return _braid_error_flag;
   }

   _braid_CoreElt(core, obj_only) = boolean;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetRNormAdjoint(braid_Core  core,
                      braid_Real  *rnorm_adj)
{
   if ( !(_braid_CoreElt(core, adjoint)) )
   {
      return _braid_error_flag;
   }

   *rnorm_adj = _braid_CoreElt(core, optim)->rnorm_adj;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetMyID(braid_Core core,
              braid_Int *myid_ptr)
{
   *myid_ptr = _braid_CoreElt(core, myid);

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

static unsigned long int _braid_rand_next = 1;
braid_Int
braid_Rand(void)
{
   _braid_rand_next = _braid_rand_next * 1103515245 + 12345;
   return (unsigned int) (_braid_rand_next/65536) % braid_RAND_MAX;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetNChunks(braid_Core core,
                 braid_Int  nchunks)
{
   _braid_CoreElt(core, nchunks) = nchunks;

   return _braid_error_flag;
}
