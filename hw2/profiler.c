
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _EXTERN_C_
#ifdef __cplusplus
#define _EXTERN_C_ extern "C"
#else /* __cplusplus */
#define _EXTERN_C_
#endif /* __cplusplus */
#endif /* _EXTERN_C_ */

#ifdef MPICH_HAS_C2F
_EXTERN_C_ void *MPIR_ToPointer(int);
#endif // MPICH_HAS_C2F

#ifdef PIC
/* For shared libraries, declare these weak and figure out which one was linked
   based on which init wrapper was called.  See mpi_init wrappers.  */
#pragma weak pmpi_init
#pragma weak PMPI_INIT
#pragma weak pmpi_init_
#pragma weak pmpi_init__
#endif /* PIC */

_EXTERN_C_ void pmpi_init(MPI_Fint *ierr);
_EXTERN_C_ void PMPI_INIT(MPI_Fint *ierr);
_EXTERN_C_ void pmpi_init_(MPI_Fint *ierr);
_EXTERN_C_ void pmpi_init__(MPI_Fint *ierr);

/* **************************************************************** *
 * Global per-node data structures and constants                    *
 * **************************************************************** */

#define USCORE "_"

/* ================== Constants for MPI Operations ================== */
#define _MPI_INIT_       0
#define _MPI_SEND_       1
#define _MPI_RECV_       2
#define _MPI_ISEND_      3
#define _MPI_IRECV_      4
#define _MPI_BARRIER_    5 
#define _MPI_SCATTER_    6
#define _MPI_GATHER_     7
#define _MPI_REDUCE_     8
#define _MPI_ALLREDUCE_  9
#define _MPI_WAIT_      10
#define _MPI_WAITALL_   11
#define _MPI_ALLTOALL_  12
#define _MPI_FINALIZE_  13

typedef struct edgeMetaData
{
    char *key;
    int weight, bytes;

} edgeMetaData;

typedef struct vertexMetaData
{
    int mpiOp;
    int fromRank, toRank;
    int tag, opSeq;
    char *key;
    char *parent;
    char *inTreeChild;
    char *intraTreeChild;
    edgeMetaData inTreeEdge;
    edgeMetaData intraTreeEdge;

} vertexMetaData;

int myRank, numNodes;

void generateKey(char *key, int op, int toRank, int tag, int opSeq)
{
    char *baseOp;
    switch(op)
    {
        case _MPI_INIT_         : baseOp = "MPI_Init"; break;
        case _MPI_SEND_         : baseOp = "MPI_Send"; break;
        case _MPI_RECV_         : baseOp = "MPI_Recv"; break;
        case _MPI_ISEND_        : baseOp = "MPI_Isend"; break;
        case _MPI_IRECV_        : baseOp = "MPI_Irecv"; break;
        case _MPI_BARRIER_      : baseOp = "MPI_Barrier"; break;
        case _MPI_SCATTER_      : baseOp = "MPI_Scatter"; break;
        case _MPI_GATHER_       : baseOp = "MPI_Gather"; break;
        case _MPI_ALLREDUCE_    : baseOp = "MPI_Allreduce"; break;
        case _MPI_WAIT_         : baseOp = "MPI_Wait"; break;
        case _MPI_WAITALL_      : baseOp = "MPI_Waitall"; break;
        case _MPI_ALLTOALL_     : baseOp = "MPI_Alltoall"; break;
        case _MPI_FINALIZE_     : baseOp = "MPI_Finalize"; break;
    }

    strcat(key, baseOp); strcat(key, USCORE);
    strcat(key, itoa(myRank)); strcat(key, USCORE);    
    strcat(key, itoa(toRank)); strcat(key, USCORE);
    strcat(key, itoa(tag)); strcat(key, USCORE);
    strcat(key, itoa(opSeq));
}

/* ================== C Wrappers for MPI_Barrier ================== */
_EXTERN_C_ int PMPI_Barrier(MPI_Comm arg_0);
_EXTERN_C_ int MPI_Barrier(MPI_Comm arg_0) 
{ 
    int _wrap_py_return_val = 0;
 
  {
    _wrap_py_return_val = PMPI_Barrier(arg_0);
  }
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Alltoall ================== */
_EXTERN_C_ int PMPI_Alltoall(void *arg_0, int arg_1, MPI_Datatype arg_2, 
                             void *arg_3, int arg_4, MPI_Datatype arg_5, 
                             MPI_Comm arg_6);
_EXTERN_C_ int MPI_Alltoall(void *arg_0, int arg_1, MPI_Datatype arg_2, 
                            void *arg_3, int arg_4, MPI_Datatype arg_5, 
                            MPI_Comm arg_6) 
{ 
    int _wrap_py_return_val = 0;
 
  {
    _wrap_py_return_val = PMPI_Alltoall(arg_0, arg_1, arg_2, arg_3, arg_4, 
                                        arg_5, arg_6);
  }
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Scatter ================== */
_EXTERN_C_ int PMPI_Scatter(void *arg_0, int arg_1, MPI_Datatype arg_2, 
                            void *arg_3, int arg_4, MPI_Datatype arg_5, 
                            int arg_6, MPI_Comm arg_7);
_EXTERN_C_ int MPI_Scatter(void *arg_0, int arg_1, MPI_Datatype arg_2, 
                           void *arg_3, int arg_4, MPI_Datatype arg_5, 
                           int arg_6, MPI_Comm arg_7) 
{ 
    int _wrap_py_return_val = 0;
 
    {
      _wrap_py_return_val = PMPI_Scatter(arg_0, arg_1, arg_2, arg_3, arg_4, 
                                         arg_5, arg_6, arg_7);
    }
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Gather ================== */
_EXTERN_C_ int PMPI_Gather(void *arg_0, int arg_1, MPI_Datatype arg_2, 
                           void *arg_3, int arg_4, MPI_Datatype arg_5, 
                           int arg_6, MPI_Comm arg_7);
_EXTERN_C_ int MPI_Gather(void *arg_0, int arg_1, MPI_Datatype arg_2, 
                          void *arg_3, int arg_4, MPI_Datatype arg_5, 
                          int arg_6, MPI_Comm arg_7) 
{ 
    int _wrap_py_return_val = 0;
 
    {
      _wrap_py_return_val = PMPI_Gather(arg_0, arg_1, arg_2, arg_3, arg_4, 
                                        arg_5, arg_6, arg_7);
    }
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Reduce ================== */
_EXTERN_C_ int PMPI_Reduce(void *arg_0, void *arg_1, int arg_2, 
                           MPI_Datatype arg_3, MPI_Op arg_4, 
                           int arg_5, MPI_Comm arg_6);
_EXTERN_C_ int MPI_Reduce(void *arg_0, void *arg_1, int arg_2, 
                          MPI_Datatype arg_3, MPI_Op arg_4, 
                          int arg_5, MPI_Comm arg_6) 
{ 
    int _wrap_py_return_val = 0;
 
    {
      _wrap_py_return_val = PMPI_Reduce(arg_0, arg_1, arg_2, arg_3, arg_4, 
                                        arg_5, arg_6);
    }
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Allreduce ================== */
_EXTERN_C_ int PMPI_Allreduce(void *arg_0, void *arg_1, int arg_2, 
                              MPI_Datatype arg_3, MPI_Op arg_4, 
                              MPI_Comm arg_5);
_EXTERN_C_ int MPI_Allreduce(void *arg_0, void *arg_1, int arg_2, 
                             MPI_Datatype arg_3, MPI_Op arg_4, 
                             MPI_Comm arg_5) 
{ 
    int _wrap_py_return_val = 0;
 
    {
      _wrap_py_return_val = PMPI_Allreduce(arg_0, arg_1, arg_2, arg_3, 
                                           arg_4, arg_5);
    }
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Send ================== */
_EXTERN_C_ int PMPI_Send(void *arg_0, int arg_1, MPI_Datatype arg_2, 
                         int arg_3, int arg_4, MPI_Comm arg_5);
_EXTERN_C_ int MPI_Send(void *arg_0, int arg_1, MPI_Datatype arg_2, 
                        int arg_3, int arg_4, MPI_Comm arg_5) 
{ 
    int _wrap_py_return_val = 0;
 
    {
      _wrap_py_return_val = PMPI_Send(arg_0, arg_1, arg_2, arg_3, arg_4, arg_5);
    }
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Isend ================== */
_EXTERN_C_ int PMPI_Isend(void *arg_0, int arg_1, MPI_Datatype arg_2, 
                          int arg_3, int arg_4, MPI_Comm arg_5, 
                          MPI_Request *arg_6);
_EXTERN_C_ int MPI_Isend(void *arg_0, int arg_1, MPI_Datatype arg_2, 
                         int arg_3, int arg_4, MPI_Comm arg_5, 
                         MPI_Request *arg_6) 
{ 
    int _wrap_py_return_val = 0;
 
    {
      _wrap_py_return_val = PMPI_Isend(arg_0, arg_1, arg_2, arg_3, arg_4, 
                                       arg_5, arg_6);
    }
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Recv ================== */
_EXTERN_C_ int PMPI_Recv(void *arg_0, int arg_1, MPI_Datatype arg_2, int arg_3,
                         int arg_4, MPI_Comm arg_5, MPI_Status *arg_6);
_EXTERN_C_ int MPI_Recv(void *arg_0, int arg_1, MPI_Datatype arg_2, int arg_3,
                        int arg_4, MPI_Comm arg_5, MPI_Status *arg_6) 
{ 
    int _wrap_py_return_val = 0;
 
    {
      _wrap_py_return_val = PMPI_Recv(arg_0, arg_1, arg_2, arg_3, arg_4, 
                                      arg_5, arg_6);
    }
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Irecv ================== */
_EXTERN_C_ int PMPI_Irecv(void *arg_0, int arg_1, MPI_Datatype arg_2,int arg_3,
                          int arg_4, MPI_Comm arg_5, MPI_Request *arg_6);
_EXTERN_C_ int MPI_Irecv(void *arg_0, int arg_1, MPI_Datatype arg_2, int arg_3,
                         int arg_4, MPI_Comm arg_5, MPI_Request *arg_6) 
{ 
    int _wrap_py_return_val = 0;
 
    {
        _wrap_py_return_val = PMPI_Irecv(arg_0, arg_1, arg_2, arg_3, arg_4, 
                                         arg_5, arg_6);
    }
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Wait ================== */
_EXTERN_C_ int PMPI_Wait(MPI_Request *arg_0, MPI_Status *arg_1);
_EXTERN_C_ int MPI_Wait(MPI_Request *arg_0, MPI_Status *arg_1) 
{ 
    int _wrap_py_return_val = 0;
 
    {
      _wrap_py_return_val = PMPI_Wait(arg_0, arg_1);
    }
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Waitall ================== */
_EXTERN_C_ int PMPI_Waitall(int arg_0, MPI_Request *arg_1, MPI_Status *arg_2);
_EXTERN_C_ int MPI_Waitall(int arg_0, MPI_Request *arg_1, MPI_Status *arg_2) 
{ 
    int _wrap_py_return_val = 0;
 
    {
      _wrap_py_return_val = PMPI_Waitall(arg_0, arg_1, arg_2);
    }
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Init ================== */
_EXTERN_C_ int PMPI_Init(int *arg_0, char ***arg_1);
_EXTERN_C_ int MPI_Init(int *arg_0, char ***arg_1) { 
    int _wrap_py_return_val = 0;
 
{
  _wrap_py_return_val = PMPI_Init(arg_0, arg_1);
}
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Finalize ================== */
_EXTERN_C_ int PMPI_Finalize();
_EXTERN_C_ int MPI_Finalize() { 
    int _wrap_py_return_val = 0;
 
{
  _wrap_py_return_val = PMPI_Finalize();
}
    return _wrap_py_return_val;
}


