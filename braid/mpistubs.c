/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision$
 ***********************************************************************EHEADER*/

/******************************************************************************
 * MPI stubs to generate serial codes without mpi
 *****************************************************************************/

#include "_braid.h"

#ifdef braid_SEQUENTIAL

MPI_Comm
MPI_Comm_f2c( int comm )
{
   return (MPI_Comm) comm;
}

MPI_Comm
MPI_Comm_c2f( int comm )
{
   return (MPI_Comm) comm;
}

int
MPI_Init( int   *argc,
                char      ***argv )
{
   return(0);
}

int
MPI_Finalize( )
{
   return(0);
}

int
MPI_Abort( MPI_Comm comm,
                 int      errorcode )
{
   return(0);
}

double
MPI_Wtime( )
{
   return(0.0);
}

double
MPI_Wtick( )
{
   return(0.0);
}

int
MPI_Barrier( MPI_Comm comm )
{
   return(0);
}

int
MPI_Comm_create( MPI_Comm   comm,
                       MPI_Group  group,
                       MPI_Comm  *newcomm )
{
   return(0);
}

int
MPI_Comm_dup( MPI_Comm  comm,
                    MPI_Comm *newcomm )
{
   return(0);
}

int
MPI_Comm_size( MPI_Comm  comm,
                     int      *size )
{ 
   *size = 1;
   return(0);
}

int
MPI_Comm_rank( MPI_Comm  comm,
                     int      *rank )
{ 
   *rank = 0;
   return(0);
}

int
MPI_Comm_free( MPI_Comm *comm )
{
   return 0;
}

int
MPI_Comm_group( MPI_Comm   comm,
                      MPI_Group *group )
{
   return(0);
}

int
MPI_Comm_split( MPI_Comm  comm,
                      int       n,
                      int       m,
                      MPI_Comm *comms )
{
   return(0);
}

int
MPI_Group_incl( MPI_Group  group,
                      int        n,
                      int       *ranks,
                      MPI_Group *newgroup )
{
   return(0);
}

int
MPI_Group_free( MPI_Group *group )
{
   return 0;
}

int
MPI_Address( void           *location,
                   MPI_Aint *address )
{
   return(0);
}

int
MPI_Get_count( MPI_Status   *status,
                     MPI_Datatype  datatype,
                     int          *count )
{
   return(0);
}

int
MPI_Alltoall( void               *sendbuf,
                    int           sendcount,
                    MPI_Datatype  sendtype,
                    void               *recvbuf,
                    int           recvcount,
                    MPI_Datatype  recvtype,
                    MPI_Comm      comm )
{
   return(0);
}

int
MPI_Allgather( void               *sendbuf,
                     int           sendcount,
                     MPI_Datatype  sendtype,
                     void               *recvbuf,
                     int           recvcount,
                     MPI_Datatype  recvtype,
                     MPI_Comm      comm ) 
{
   int i;

   switch (sendtype)
   {
      case MPI_INT:
      {
         int *crecvbuf = (int *)recvbuf;
         int *csendbuf = (int *)sendbuf;
         for (i = 0; i < sendcount; i++)
         {
	    crecvbuf[i] = csendbuf[i];
         }
      } 
      break;

      case MPI_DOUBLE:
      {
         double *crecvbuf = (double *)recvbuf;
         double *csendbuf = (double *)sendbuf;
         for (i = 0; i < sendcount; i++)
         {
	    crecvbuf[i] = csendbuf[i];
         }
      } 
      break;

      case MPI_CHAR:
      {
         char *crecvbuf = (char *)recvbuf;
         char *csendbuf = (char *)sendbuf;
         for (i = 0; i < sendcount; i++)
         {
	    crecvbuf[i] = csendbuf[i];
         }
      } 
      break;

      case MPI_BYTE:
      {
         memcpy(recvbuf, sendbuf, sendcount);
      } 
      break;

      case MPI_REAL:
      {
         double *crecvbuf = (double *)recvbuf;
         double *csendbuf = (double *)sendbuf;
         for (i = 0; i < sendcount; i++)
         {
	    crecvbuf[i] = csendbuf[i];
         }
      } 
      break;

      case MPI_FLOAT:
      {
         float *crecvbuf = (float *)recvbuf;
         float *csendbuf = (float *)sendbuf;
         for (i = 0; i < sendcount; i++)
         {
	    crecvbuf[i] = csendbuf[i];
         }
      } 
      break;

   }

   return(0);
}

int
MPI_Allgatherv( void               *sendbuf,
                      int           sendcount,
                      MPI_Datatype  sendtype,
                      void               *recvbuf,
                      int          *recvcounts,
                      int          *displs, 
                      MPI_Datatype  recvtype,
                      MPI_Comm      comm ) 
{ 
   return ( MPI_Allgather(sendbuf, sendcount, sendtype,
                                recvbuf, *recvcounts, recvtype, comm) );
}

int
MPI_Gather( void               *sendbuf,
                  int           sendcount,
                  MPI_Datatype  sendtype,
                  void               *recvbuf,
                  int           recvcount,
                  MPI_Datatype  recvtype,
                  int           root,
                  MPI_Comm      comm )
{
   return ( MPI_Allgather(sendbuf, sendcount, sendtype,
                                recvbuf, recvcount, recvtype, comm) );
}

int
MPI_Gatherv( void              *sendbuf,
                  int           sendcount,
                  MPI_Datatype  sendtype,
                  void               *recvbuf,
                  int          *recvcounts,
                  int          *displs,
                  MPI_Datatype  recvtype,
                  int           root,
                  MPI_Comm      comm )
{
   return ( MPI_Allgather(sendbuf, sendcount, sendtype,
                                recvbuf, *recvcounts, recvtype, comm) );
}

int
MPI_Scatter( void               *sendbuf,
                   int           sendcount,
                   MPI_Datatype  sendtype,
                   void               *recvbuf,
                   int           recvcount,
                   MPI_Datatype  recvtype,
                   int           root,
                   MPI_Comm      comm )
{
   return ( MPI_Allgather(sendbuf, sendcount, sendtype,
                                recvbuf, recvcount, recvtype, comm) );
}

int
MPI_Scatterv( void               *sendbuf,
                   int           *sendcounts,
                   int           *displs,
                   MPI_Datatype   sendtype,
                   void                *recvbuf,
                   int            recvcount,
                   MPI_Datatype   recvtype,
                   int            root,
                   MPI_Comm       comm )
{
   return ( MPI_Allgather(sendbuf, *sendcounts, sendtype,
                                recvbuf, recvcount, recvtype, comm) );
}

int
MPI_Bcast( void               *buffer,
                 int           count,
                 MPI_Datatype  datatype,
                 int           root,
                 MPI_Comm      comm ) 
{ 
   return(0);
}

int
MPI_Send( void               *buf,
                int           count,
                MPI_Datatype  datatype,
                int           dest,
                int           tag,
                MPI_Comm      comm ) 
{ 
   return(0);
}

int
MPI_Recv( void               *buf,
                int           count,
                MPI_Datatype  datatype,
                int           source,
                int           tag,
                MPI_Comm      comm,
                MPI_Status   *status )
{ 
   return(0);
}

int
MPI_Isend( void               *buf,
                 int           count,
                 MPI_Datatype  datatype,
                 int           dest,
                 int           tag,
                 MPI_Comm      comm,
                 MPI_Request  *request )
{ 
   return(0);
}

int
MPI_Irecv( void               *buf,
                 int           count,
                 MPI_Datatype  datatype,
                 int           source,
                 int           tag,
                 MPI_Comm      comm,
                 MPI_Request  *request )
{ 
   return(0);
}

int
MPI_Send_init( void               *buf,
                     int           count,
                     MPI_Datatype  datatype,
                     int           dest,
                     int           tag, 
                     MPI_Comm      comm,
                     MPI_Request  *request )
{
   return 0;
}

int
MPI_Recv_init( void               *buf,
                     int           count,
                     MPI_Datatype  datatype,
                     int           dest,
                     int           tag, 
                     MPI_Comm      comm,
                     MPI_Request  *request )
{
   return 0;
}

int
MPI_Irsend( void               *buf,
                  int           count,
                  MPI_Datatype  datatype,
                  int           dest,
                  int           tag, 
                  MPI_Comm      comm,
                  MPI_Request  *request )
{
   return 0;
}

int
MPI_Startall( int          count,
                    MPI_Request *array_of_requests )
{
   return 0;
}

int
MPI_Probe( int         source,
                 int         tag,
                 MPI_Comm    comm,
                 MPI_Status *status )
{
   return 0;
}

int
MPI_Iprobe( int         source,
                  int         tag,
                  MPI_Comm    comm,
                  int        *flag,
                  MPI_Status *status )
{
   return 0;
}

int
MPI_Test( MPI_Request *request,
                int         *flag,
                MPI_Status  *status )
{
   *flag = 1;
   return(0);
}

int
MPI_Testall( int          count,
                   MPI_Request *array_of_requests,
                   int         *flag,
                   MPI_Status  *array_of_statuses )
{
   *flag = 1;
   return(0);
}

int
MPI_Wait( MPI_Request *request,
                MPI_Status  *status )
{
   return(0);
}

int
MPI_Waitall( int          count,
                   MPI_Request *array_of_requests,
                   MPI_Status  *array_of_statuses )
{
   return(0);
}

int
MPI_Waitany( int          count,
                   MPI_Request *array_of_requests,
                   int         *index,
                   MPI_Status  *status )
{
   return(0);
}

int
MPI_Allreduce( void              *sendbuf,
                     void              *recvbuf,
                     int          count,
                     MPI_Datatype datatype,
                     MPI_Op       op,
                     MPI_Comm     comm )
{ 
   int i;
   
   switch (datatype)
   {
      case MPI_INT:
      {
         int *crecvbuf = (int *)recvbuf;
         int *csendbuf = (int *)sendbuf;
         for (i = 0; i < count; i++)
         {
            crecvbuf[i] = csendbuf[i];
         }
         
      } 
      break;

      case MPI_DOUBLE:
      {
         double *crecvbuf = (double *)recvbuf;
         double *csendbuf = (double *)sendbuf;
         for (i = 0; i < count; i++)
         {
            crecvbuf[i] = csendbuf[i];
         }
      } 
      break;

      case MPI_CHAR:
      {
         char *crecvbuf = (char *)recvbuf;
         char *csendbuf = (char *)sendbuf;
         for (i = 0; i < count; i++)
         {
            crecvbuf[i] = csendbuf[i];
         }
      } 
      break;

      case MPI_BYTE:
      {
         memcpy(recvbuf, sendbuf, count);
      } 
      break;

      case MPI_REAL:
      {
         double *crecvbuf = (double *)recvbuf;
         double *csendbuf = (double *)sendbuf;
         for (i = 0; i < count; i++)
         {
            crecvbuf[i] = csendbuf[i];
         }
      } 
      break;

   }

   return 0;
}

int
MPI_Reduce( void               *sendbuf,
                  void               *recvbuf,
                  int           count,
                  MPI_Datatype  datatype,
                  MPI_Op        op,
                  int           root,
                  MPI_Comm      comm )
{ 
   MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
   return 0;
}

int
MPI_Scan( void               *sendbuf,
                void               *recvbuf,
                int           count,
                MPI_Datatype  datatype,
                MPI_Op        op,
                MPI_Comm      comm )
{ 
   MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
   return 0;
}

int
MPI_Request_free( MPI_Request *request )
{
   return 0;
}

int
MPI_Type_contiguous( int           count,
                           MPI_Datatype  oldtype,
                           MPI_Datatype *newtype )
{
   return(0);
}

int
MPI_Type_vector( int           count,
                       int           blocklength,
                       int           stride,
                       MPI_Datatype  oldtype,
                       MPI_Datatype *newtype )
{
   return(0);
}

int
MPI_Type_hvector( int           count,
                        int           blocklength,
                        MPI_Aint      stride,
                        MPI_Datatype  oldtype,
                        MPI_Datatype *newtype )
{
   return(0);
}

int
MPI_Type_struct( int           count,
                       int          *array_of_blocklengths,
                       MPI_Aint     *array_of_displacements,
                       MPI_Datatype *array_of_types,
                       MPI_Datatype *newtype )
{
   return(0);
}

int
MPI_Type_commit( MPI_Datatype *datatype )
{
   return(0);
}

int
MPI_Type_free( MPI_Datatype *datatype )
{
   return(0);
}

int
MPI_Op_free( MPI_Op *op )
{
   return(0);
}

#endif
