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
 * Returns an index into the local u-vector for grid 'level' at point 'index'.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UGetIndex(braid_Core   core,
                 braid_Int    level,
                 braid_Int    index,
                 braid_Int   *uindex_ptr,
                 braid_Int   *store_flag_ptr)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Int            ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper      = _braid_GridElt(grids[level], iupper);
   braid_Int            clower      = _braid_GridElt(grids[level], clower);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_Int            uindex, ic, iclo, store_flag;

   uindex = -1;
   store_flag = -2;
   if ((index >= ilower) && (index <= iupper))
   {
      if ( _braid_CoreElt(core, useshell) == 1)
      {
         uindex = index-ilower;
         store_flag = 0;
         // If we are not on a fully-stored point
         // then we only have a shell, the store_flag should be -1
         if ( (_braid_CoreElt(core, storage) < 0) ||
              (level < _braid_CoreElt(core, storage)) )
         {
            if ( !_braid_IsCPoint(index, cfactor) )
            {
               store_flag = -1;
            }
         }
      }
      else
      {
         // If on level that only stores C-points
         if ( (_braid_CoreElt(core, storage) < 0) ||
              (level < _braid_CoreElt(core, storage)) )
         {
            if ( _braid_IsCPoint(index, cfactor) )
            {
               _braid_MapFineToCoarse(index, cfactor, ic);
               _braid_MapFineToCoarse(clower, cfactor, iclo);
               uindex = ic-iclo;
               store_flag = 0;
            }
         }
         else
         {
            uindex = index-ilower;
            store_flag = 0;
         }
      }
   }

   *uindex_ptr = uindex;
   *store_flag_ptr = store_flag;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns a reference to the local u-vector on grid 'level' at point 'index'.
 * If the u-vector is not stored, returns NULL. The referred u-vector might
 * just be a shell if that feature is used.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UGetVectorRef(braid_Core         core,
                     braid_Int          level,
                     braid_Int          index,
                     braid_BaseVector  *u_ptr)
{
   _braid_Grid        **grids = _braid_CoreElt(core, grids);
   braid_BaseVector    *ua    = _braid_GridElt(grids[level], ua);
   braid_BaseVector     u     = NULL;
   braid_Int            iu, sflag;

   _braid_UGetIndex(core, level, index, &iu, &sflag);
   if (sflag>-2) // We have a full point or a shell (iu>=0)
   {
      u = ua[iu];
   }

   *u_ptr = u;

   return _braid_error_flag;
}


/*----------------------------------------------------------------------------
 * Returns a reference to the vector at the last time step.
 * Return NULL if it is not stored on this processor.
 *---------------------------------------------------------------------------*/
braid_Int
_braid_UGetLast(braid_Core        core,
                braid_BaseVector *u_ptr)
{
   _braid_Grid        **grids = _braid_CoreElt(core, grids);
   int                 ntime = _braid_CoreElt(core, ntime);
   braid_BaseVector    ulast;

   /* Get last time step */
   if (_braid_CoreElt(core, storage) < 0 )
   {
      ulast  = _braid_GridElt(grids[0], ulast);
   }
   else
   {
     _braid_UGetVectorRef(core, 0, ntime, &ulast);
   }

  *u_ptr = ulast;

   return _braid_error_flag;
}


/*----------------------------------------------------------------------------
 * Stores a reference to the local u-vector on grid 'level' at point 'index'.
 * If the shellvector feature is used, the u-vector might be emptied so that
 * only the shell is stored. Otherwise, if the u-vector is not stored, nothing
 * is done.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_USetVectorRef(braid_Core        core,
                     braid_Int         level,
                     braid_Int         index,
                     braid_BaseVector  u)
{
   _braid_Grid        **grids = _braid_CoreElt(core, grids);
   braid_BaseVector    *ua    = _braid_GridElt(grids[level], ua);
   braid_Int            iu, sflag;

   _braid_UGetIndex(core, level, index, &iu, &sflag);
   // If sflag ==0, we have a full point, if sflag == -1, we have a shell
   if (sflag == 0)
   {
      ua[iu] = u;
   }
   else if (sflag == -1)
   {
      braid_App    app = _braid_CoreElt(core, app);
      _braid_BaseSFree(core,  app, u);
      ua[iu] = u;
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns a copy of the u-vector on grid 'level' at point 'index'.  If 'index'
 * is my "receive index" (as set by UCommInit(), for example), the u-vector will
 * be received from a neighbor processor.  If the u-vector is not stored, NULL
 * is returned. The copy might just be a shell if this feature is used.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UGetVector(braid_Core         core,
                  braid_Int          level,
                  braid_Int          index,
                  braid_BaseVector  *u_ptr)
{
   braid_App            app         = _braid_CoreElt(core, app);
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_BaseVector    *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            recv_index  = _braid_GridElt(grids[level], recv_index);
   _braid_CommHandle   *recv_handle = _braid_GridElt(grids[level], recv_handle);
   braid_BaseVector     u           = NULL;
   braid_Int            iu, sflag;

   if (index == recv_index)
   {
      /* If a recv was initiated, receive u value from neighbor processor */
      if (recv_index > -1)
      {
         _braid_CommWait(core, &recv_handle);
         _braid_GridElt(grids[level], recv_index)  = -1;
         _braid_GridElt(grids[level], recv_handle) = recv_handle;
         u = ua[-1];
      }
   }
   else
   {
      _braid_UGetIndex(core, level, index, &iu, &sflag);
      if (sflag == 0)
      {
         _braid_BaseClone(core, app, ua[iu], &u);
      }
      else if (sflag == -1)
      {
         // In this case, sclone != NULL
         _braid_BaseSClone(core, app, ua[iu], &u);
      }
   }

   *u_ptr = u;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Stores the u-vector on grid 'level' at point 'index'.  If 'index' is my "send
 * index", a send is initiated to a neighbor processor.  If 'move' is true, the
 * u-vector is moved into core storage instead of copied.  If the u-vector is
 * not stored, nothing is done or only the shell is copied/moved when the shellvector
 * feature is used.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_USetVector(braid_Core        core,
                  braid_Int         level,
                  braid_Int         index,
                  braid_BaseVector  u,
                  braid_Int         move)
{
   braid_App            app         = _braid_CoreElt(core, app);
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_BaseVector    *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            send_index  = _braid_GridElt(grids[level], send_index);
   _braid_CommHandle   *send_handle = _braid_GridElt(grids[level], send_handle);
   braid_Int            iu, sflag;

   if (index == send_index)
   {
      /* Post send to neighbor processor */
      _braid_CommSendInit(core, level, index, u, &send_handle);
      _braid_GridElt(grids[level], send_index)  = -1;
      _braid_GridElt(grids[level], send_handle) = send_handle;
   }

   _braid_UGetIndex(core, level, index, &iu, &sflag);
   if (sflag == 0) // We have a full point
   {
      if (ua[iu] != NULL)
      {
         _braid_BaseFree(core, app,  ua[iu]);
      }
      if (move)
      {
         ua[iu] = u;                                   /* move the vector */
      }
      else
      {
         _braid_BaseClone(core, app,  u, &ua[iu]); /* copy the vector */
      }
   }
   else if (sflag == -1) // We have a shell
   {
      if (ua[iu] != NULL)
      {
         _braid_BaseFree(core, app,  ua[iu]);
      }
      if (move)
      {
         // We are on an F-point, with shellvector option. We only keep the shell.
         _braid_BaseSFree(core,  app, u);
         ua[iu] = u;                                   /* move the vector */
      }
      else
      {
         _braid_BaseSClone(core,  app, u, &ua[iu]); /* copy the vector */
      }
   }
   else if (move) // We store nothing
   {
      _braid_BaseFree(core, app,  u);              /* free the vector */
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Basic communication (from the left, to the right)
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UCommInitBasic(braid_Core  core,
                      braid_Int   level,
                      braid_Int   recv_msg,
                      braid_Int   send_msg,
                      braid_Int   send_now)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Int            ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper      = _braid_GridElt(grids[level], iupper);
   braid_BaseVector    *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            recv_index  = -1;
   braid_Int            send_index  = -1;
   _braid_CommHandle   *recv_handle = NULL;
   _braid_CommHandle   *send_handle = NULL;
   braid_Int            iu, sflag;

   /* Post receive */
   if (recv_msg)
   {
      recv_index = ilower-1;
      _braid_CommRecvInit(core, level, recv_index, &ua[-1], &recv_handle);
   }

   /* Post send */
   if (send_msg)
   {
      send_index = iupper;
      if (send_now)
      {
         _braid_UGetIndex(core, level, send_index, &iu, &sflag);
         if (sflag < 0)
         {
            // We should never get here : we do not communicate shells...
            abort();
         }
         _braid_CommSendInit(core, level, send_index, ua[iu], &send_handle);
         send_index = -1;
      }
   }

   _braid_GridElt(grids[level], recv_index)  = recv_index ;
   _braid_GridElt(grids[level], send_index)  = send_index ;
   _braid_GridElt(grids[level], recv_handle) = recv_handle;
   _braid_GridElt(grids[level], send_handle) = send_handle;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Working on all intervals
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UCommInit(braid_Core  core,
                 braid_Int   level)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Int            ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper      = _braid_GridElt(grids[level], iupper);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_BaseVector    *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            recv_index  = -1;
   braid_Int            send_index  = -1;
   _braid_CommHandle   *recv_handle = NULL;
   _braid_CommHandle   *send_handle = NULL;
   braid_Int            iu, sflag;
   
   /* Note that this routine works for the case of all points being C-points, 
    * i.e., cfactor = 1.  A send and receive are always posted. */

   if (ilower <= iupper)
   {
      /* Post receive */
      _braid_CommRecvInit(core, level, ilower-1, &ua[-1], &recv_handle);
      recv_index = ilower-1;
      
      /* Only post send if iupper is a C-point, otherwise compute and send later */
      if ( _braid_IsCPoint(iupper, cfactor) )
      {
         _braid_UGetIndex(core, level, iupper, &iu, &sflag);
         _braid_CommSendInit(core, level, iupper, ua[iu], &send_handle);
         send_index = -1;
      }
      else
      {
         send_index = iupper;
      }
   }

   _braid_GridElt(grids[level], recv_index)  = recv_index ;
   _braid_GridElt(grids[level], send_index)  = send_index ;
   _braid_GridElt(grids[level], recv_handle) = recv_handle;
   _braid_GridElt(grids[level], send_handle) = send_handle;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Working only on F-pt intervals
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UCommInitF(braid_Core  core,
                  braid_Int   level)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Int            ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper      = _braid_GridElt(grids[level], iupper);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_BaseVector    *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            recv_index  = -1;
   braid_Int            send_index  = -1;
   _braid_CommHandle   *recv_handle = NULL;
   _braid_CommHandle   *send_handle = NULL;
   braid_Int            iu, sflag;

   if (ilower <= iupper)
   {
      /* Only post receive if ilower is an F-point */
      if ( _braid_IsFPoint(ilower, cfactor) )
      {
         _braid_CommRecvInit(core, level, ilower-1, &ua[-1], &recv_handle);
         recv_index = ilower-1;
      }

      /* Only post send if iupper is a C-point and iupper+1 is an F-point.  This
       * check allows for the case of cfactor=1, i.e., all C-points.  Otherwise, 
       * if iupper+1 is an F-point, set send_index, so that when that point is 
       * computed later, it is sent. */
      if ( _braid_IsCPoint(iupper, cfactor) && _braid_IsFPoint(iupper+1, cfactor))
      {
         _braid_UGetIndex(core, level, iupper, &iu, &sflag);
         _braid_CommSendInit(core, level, iupper, ua[iu], &send_handle);
         send_index = -1;
      }
      else if ( _braid_IsFPoint(iupper+1, cfactor) )
      {
         send_index = iupper;
      }
   }

   _braid_GridElt(grids[level], recv_index)  = recv_index ;
   _braid_GridElt(grids[level], send_index)  = send_index ;
   _braid_GridElt(grids[level], recv_handle) = recv_handle;
   _braid_GridElt(grids[level], send_handle) = send_handle;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Finish up communication
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UCommWait(braid_Core  core,
                 braid_Int   level)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   _braid_CommHandle   *recv_handle = _braid_GridElt(grids[level], recv_handle);
   _braid_CommHandle   *send_handle = _braid_GridElt(grids[level], send_handle);

   _braid_CommWait(core, &recv_handle);
   _braid_CommWait(core, &send_handle);
   _braid_GridElt(grids[level], recv_index)  = -1;
   _braid_GridElt(grids[level], send_index)  = -1;
   _braid_GridElt(grids[level], recv_handle) = recv_handle;
   _braid_GridElt(grids[level], send_handle) = send_handle;

   return _braid_error_flag;
}

