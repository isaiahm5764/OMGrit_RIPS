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
 * Compute residual
 *----------------------------------------------------------------------------*/

braid_Int
_braid_Residual(braid_Core        core,
                braid_Int         level,
                braid_Int         index,
                braid_BaseVector  ustop,
                braid_BaseVector  r)
{
   braid_App        app      = _braid_CoreElt(core, app);
   braid_Real       tol      = _braid_CoreElt(core, tol);
   braid_Int        iter     = _braid_CoreElt(core, niter);
   braid_Int       *rfactors = _braid_CoreElt(core, rfactors);
   _braid_Grid    **grids    = _braid_CoreElt(core, grids);
   braid_StepStatus status   = (braid_StepStatus)core;
   braid_Int        nrefine  = _braid_CoreElt(core, nrefine);
   braid_Int        gupper   = _braid_CoreElt(core, gupper);
   braid_Int        ilower   = _braid_GridElt(grids[level], ilower);
   braid_Int        ichunk   = _braid_CoreElt(core, ichunk);
   braid_Real      *ta       = _braid_GridElt(grids[level], ta);

   braid_BaseVector rstop;
   braid_Int        ii;

   ii = index-ilower;
   _braid_StepStatusInit(ta[ii-1], ta[ii], index-1, ichunk, tol, iter, level, nrefine, gupper, status);
   if ( _braid_CoreElt(core, residual) == NULL )
   {
      /* By default: r = ustop - \Phi(ustart)*/
      _braid_GetUInit(core, level, index, r, &rstop);
      _braid_BaseStep(core, app,  rstop, NULL, r, level, status);
      _braid_BaseSum(core, app,  1.0, ustop, -1.0, r);
      if (level == 0)
      {
         /*TODO Remove this line after modifing the _braid_StatusSetRFactor to set the rfactor in the array directly */
         rfactors[ii] = _braid_StatusElt(status, rfactor);
         /* TODO : Remove these two lines, which are now useless since core==status */
         if ( !_braid_CoreElt(core, r_space) && _braid_StatusElt(status, r_space) )
               _braid_CoreElt(core, r_space) = 1;
      }
   }
   else
   {
      /* Call the user's residual routine */
      _braid_BaseResidual(core, app, ustop, r, status);
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Compute FAS residual = f - residual
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FASResidual(braid_Core        core,
                   braid_Int         level,
                   braid_Int         index,
                   braid_BaseVector  ustop,
                   braid_BaseVector  r)
{
   braid_App          app    = _braid_CoreElt(core, app);
   _braid_Grid      **grids  = _braid_CoreElt(core, grids);
   braid_Int          ilower = _braid_GridElt(grids[level], ilower);
   braid_BaseVector  *fa     = _braid_GridElt(grids[level], fa);

   braid_Int        ii;

   _braid_Residual(core, level, index, ustop, r);
   if (level == 0)
   {
      _braid_BaseSum(core, app,  0.0, r, -1.0, r);
   }
   else
   {
      ii = index-ilower;
      if(fa[ii] == NULL)
      {
         _braid_BaseSum(core, app,  0.0, r, -1.0, r);
      }
      else
      {
         _braid_BaseSum(core, app,  1.0, fa[ii], -1.0, r);
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Compute residual A(u) for time step 'index' on grid 'level'
 *----------------------------------------------------------------------------*/

braid_Int
_braid_TriResidual(braid_Core         core,
                   braid_Int          level,
                   braid_Int          index,
                   braid_Int          fas,          /* FAS residual? */ 
                   braid_BaseVector  *r_ptr)
{
   braid_App          app      = _braid_CoreElt(core, app);
   _braid_Grid      **grids    = _braid_CoreElt(core, grids);
   braid_TriStatus    status   = (braid_TriStatus)core;
   braid_Int          ilower   = _braid_GridElt(grids[level], ilower);
   braid_Real        *ta       = _braid_GridElt(grids[level], ta);
   braid_BaseVector  *fa       = _braid_GridElt(grids[level], fa);

   braid_BaseVector   u, uleft, uright, r;

   braid_Int          ii = index-ilower, homogeneous = 0;

   /* Update status (core) */
   _braid_StatusElt(status, t)     = ta[ii];
   _braid_StatusElt(status, tprev) = ta[ii-1];
   _braid_StatusElt(status, tnext) = ta[ii+1];
   _braid_StatusElt(status, idx)   = index;
   _braid_StatusElt(status, level) = level;
   
   /* Compute residual */

   _braid_UGetVectorRef(core, level, index-1, &uleft);
   _braid_UGetVectorRef(core, level, index+1, &uright);
   _braid_UGetVectorRef(core, level, index, &u);
   _braid_BaseClone(core, app, u, &r);


   if (level > 0)
   {
      homogeneous = 1;
   }

   if ( (level == 0) || (!fas) )
   {
      /* No FAS rhs */
      _braid_BaseTriResidual(core, app, uleft, uright, NULL, r, homogeneous, status);
   }
   else
   {
      _braid_BaseTriResidual(core, app, uleft, uright, fa[ii], r, homogeneous, status);
   }

   *r_ptr = r;

   return _braid_error_flag;
}

