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

#include "_braid.h"
#include "_util.h"

/*----------------------------------------------------------------------------
 * Integrate one time step
 *----------------------------------------------------------------------------*/

braid_Int
_braid_Step(braid_Core         core,
            braid_Int          level,
            braid_Int          index,
            braid_BaseVector   ustop,
            braid_BaseVector   u)
{
   braid_App          app      = _braid_CoreElt(core, app);
   braid_Real         tol      = _braid_CoreElt(core, tol);
   braid_Int          iter     = _braid_CoreElt(core, niter);
   braid_Int          ichunk   = _braid_CoreElt(core, ichunk);
   braid_Int         *rfactors = _braid_CoreElt(core, rfactors);
   _braid_Grid      **grids    = _braid_CoreElt(core, grids);
   braid_StepStatus   status   = (braid_StepStatus)core;
   braid_Int          nrefine  = _braid_CoreElt(core, nrefine);
   braid_Int          gupper   = _braid_CoreElt(core, gupper);
   braid_Int          ilower   = _braid_GridElt(grids[level], ilower);
   braid_Real        *ta       = _braid_GridElt(grids[level], ta);
   braid_BaseVector  *fa       = _braid_GridElt(grids[level], fa);

   braid_Int        ii;

   ii = index-ilower;
   _braid_StepStatusInit(ta[ii-1], ta[ii], index-1, ichunk, tol, iter, level, nrefine, gupper, status);

   /* If ustop is set to NULL, use a default approach for setting it */
   if (ustop == NULL)
   {
      _braid_GetUInit(core, level, index, u, &ustop);
   }

   if (level == 0)
   {
      _braid_BaseStep(core, app,  ustop, NULL, u, level, status);
      rfactors[ii] = _braid_StatusElt(status, rfactor);
      if ( !_braid_CoreElt(core, r_space) && _braid_StatusElt(status, r_space) )
            _braid_CoreElt(core, r_space) = 1;
   }
   else
   {
      if ( _braid_CoreElt(core, residual) == NULL )
      {
         _braid_BaseStep(core, app,  ustop, NULL, u, level, status);
         if(fa[ii] != NULL)
         {
            _braid_BaseSum(core, app,  1.0, fa[ii], 1.0, u);
         }
      }
      else
      {
         _braid_BaseStep(core, app,  ustop, fa[ii], u, level, status);
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Get an initial guess for ustop to use in the step routine (implicit schemes)
 * This vector may just be a shell. User should be able to deal with it
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetUInit(braid_Core         core,
                braid_Int          level,
                braid_Int          index,
                braid_BaseVector   u,
                braid_BaseVector  *ustop_ptr)
{
   _braid_Grid       **grids    = _braid_CoreElt(core, grids);
   braid_Int           ilower   = _braid_GridElt(grids[level], ilower);
   braid_Int           storage  = _braid_CoreElt(core, storage);
   braid_BaseVector   *va       = _braid_GridElt(grids[level], va);
   braid_BaseVector    ustop    = *ustop_ptr;
   braid_Int        ii;

   ii = index-ilower;

   _braid_UGetVectorRef(core, level, index, &ustop);

   /* If ustop is NULL, then storage is only at C-points on this level and this
    * is an F-point.  See the comment block around FRestrict() for the fixed-point
    * logic behind our choices in ustop. */
   if( ustop == NULL)
   {
      if( (level == 0) || ( storage == -2 ) )
      {
         ustop = u;
      }
      else
      {
         ustop = va[ii];
      }
   }

   /* If you have storage at this point, use it, unless you're in compatibility mode (-2). */
   else if( storage == -2 )
   {
      if ( _braid_CoreElt(core, useshell) == 1)
      {
         // Should not happen, ustop is never NULL with useshell option
         // unless there are inconsistent options (i.e. useshell && storage==-2)
         abort();
      }
      ustop = u;
   }

   *ustop_ptr = ustop;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Solve A(u) for time step 'index' on grid 'level'
 *----------------------------------------------------------------------------*/

braid_Int
_braid_TriSolve(braid_Core  core,
                braid_Int   level,
                braid_Int   index)
{
   braid_App          app      = _braid_CoreElt(core, app);
   _braid_Grid      **grids    = _braid_CoreElt(core, grids);
   braid_TriStatus    status   = (braid_TriStatus)core;
   braid_Int          ilower   = _braid_GridElt(grids[level], ilower);
   braid_Real        *ta       = _braid_GridElt(grids[level], ta);
   braid_BaseVector  *fa       = _braid_GridElt(grids[level], fa);

   braid_BaseVector   u, uleft, uright;

   braid_Int          ii = index-ilower, homogeneous = 0;

   /* Update status (core) */
   _braid_StatusElt(status, t)     = ta[ii];
   _braid_StatusElt(status, tprev) = ta[ii-1];
   _braid_StatusElt(status, tnext) = ta[ii+1];
   _braid_StatusElt(status, idx)   = index;
   _braid_StatusElt(status, level) = level;

   /* Solve A(u) */

   _braid_UGetVectorRef(core, level, index-1, &uleft);
   _braid_UGetVectorRef(core, level, index+1, &uright);
   _braid_UGetVectorRef(core, level, index, &u);

   if (level > 0)
   {
      homogeneous = 1;
   }

   if (level == 0)
   {
      /* No FAS rhs */
      _braid_BaseTriSolve(core, app, uleft, uright, NULL, u, homogeneous, status);
   }
   else
   {
      _braid_BaseTriSolve(core, app, uleft, uright, fa[ii], u, homogeneous, status);
   }

   _braid_USetVectorRef(core, level, index, u);

   return _braid_error_flag;
}

