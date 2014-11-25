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

/* ================== Global definitions ================== */

#define PROTOCOL_TYPE "PROTOCOL_TYPE"
#define MIRRORED 0
#define PARALLEL 1
#define ALIVE 1
#define DEAD  0
#define BARRIER_TAG 0xb001

int protoType;
int trueRank, trueNumNodes;
int fakeRank, fakeNumNodes;
int barrierCount = 0, isLeaderNode = 0;
int *alive, *scheduleDeath;

int My_Barrier(MPI_Comm comm)
{
    int err = 0, i, size = 0, firstAlive = -1, nextAlive = -1;
    for(i = 0; i < trueNumNodes; ++i)
    {
        if(alive[i] == ALIVE)
        {   
            if(firstAlive != -1 && nextAlive == -1)
                nextAlive = i;
            if(firstAlive == -1)
                firstAlive = i;
            ++size;
        }
    }

    if(trueRank != firstAlive)
    {
        err = PMPI_Send(NULL, 0, MPI_BYTE, firstAlive, BARRIER_TAG, comm);
        if(err !=  MPI_SUCCESS)
            return err;
        err = PMPI_Recv(NULL, 0, MPI_BYTE, firstAlive, BARRIER_TAG, comm, 
                        MPI_STATUS_IGNORE);
        if(err !=  MPI_SUCCESS)
            return err;
    }
    else
    {
        for(i = nextAlive; i < trueNumNodes; ++i)
        {
            if(alive[i])
            {
                err = PMPI_Recv(NULL, 0, MPI_BYTE, i, BARRIER_TAG, comm, 
                                MPI_STATUS_IGNORE);
                if(err !=  MPI_SUCCESS)
                    return err;
            }
        }

        for(i = nextAlive; i < trueNumNodes; ++i)
        {
            if(alive[i])
            {
                err = PMPI_Send(NULL, 0, MPI_BYTE, i, BARRIER_TAG, comm);
                if(err !=  MPI_SUCCESS)
                    return err;            
            }
        }
    }

    return MPI_SUCCESS;
}

/* ================== C Wrappers for MPI_Barrier ================== */
_EXTERN_C_ int PMPI_Barrier(MPI_Comm arg_0);
_EXTERN_C_ int MPI_Barrier(MPI_Comm arg_0) { 
    int _wrap_py_return_val = 0, i;

    barrierCount++;
    for(i = 0; i < trueNumNodes; ++i)
        if(scheduleDeath[i] == DEAD)
            alive[i] = DEAD;

    printf("\nAlive vector rank#%d ", trueRank);
    for(i = 0; i < trueNumNodes; ++i)
        printf("%d ", alive[i]);
    printf("\n");

    if(alive[trueRank] == DEAD)
    {
        printf("\n%d Calling MPI_Finalize", trueRank);
        MPI_Finalize();
        exit(0);
    }
    else
    {
        printf("\n%d Barrier", trueRank);
        _wrap_py_return_val = My_Barrier(arg_0);
        printf("\n%d out of barrier", trueRank);
    }
    
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Pcontrol ================== */
_EXTERN_C_ int PMPI_Pcontrol(int level, ...);
_EXTERN_C_ int MPI_Pcontrol(int level, ...) {
    int _wrap_py_return_val = 0; int i;

{
  _wrap_py_return_val = PMPI_Pcontrol(level);
}
    if(level < trueNumNodes)
        scheduleDeath[level] = DEAD;

    printf("\nSchedule Death vector rank#%d ", trueRank);
    for(i = 0; i < trueNumNodes; ++i)
        printf("%d ", scheduleDeath[i]);
    printf("\n");

    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Send ================== */
_EXTERN_C_ int PMPI_Send(void *buf, int cnt, MPI_Datatype datatype, int dest, 
                         int tag, MPI_Comm comm);
_EXTERN_C_ int MPI_Send(void *buf, int cnt, MPI_Datatype datatype, int dest, 
                         int tag, MPI_Comm comm) 
{
    int _wrap_py_return_val = 0;
    if(!isLeaderNode)
        dest = dest + fakeNumNodes;

    // MIRRORED PROTOCOL
    if(protoType == MIRRORED)
    {
        if(isLeaderNode)
        {
            // Send to next node first
            if(alive[dest])
            {
                printf("\nRank#%d MPI_Send %d->%d tag %d", trueRank, 
                   trueRank, dest, tag);
                _wrap_py_return_val = PMPI_Send(buf, cnt, datatype, dest, 
                                                tag, comm);
            }

            //Get replica node rank and send if alive
            dest = dest + fakeNumNodes;
            
            if(alive[dest])
            {
                printf("\nRank#%d MPI_Send %d->%d tag %d", trueRank, 
                    trueRank, dest, tag);
                _wrap_py_return_val = PMPI_Send(buf, cnt, datatype, dest, 
                                                tag, comm);
            }
        }
        else
        {
            // Send to next node first
            if(alive[dest])
            {
                printf("\nRank#%d MPI_Send %d->%d tag %d", trueRank, trueRank, 
                    dest, tag + 100);
                _wrap_py_return_val = PMPI_Send(buf, cnt, datatype, dest, 
                                                tag + 100, comm);
            }

            //Get replica node rank and send if alive
            dest = dest - fakeNumNodes;
            
            if(alive[dest])
            {
                printf("\nRank#%d MPI_Send %d->%d tag %d", trueRank, trueRank, 
                    dest, tag + 100);
                _wrap_py_return_val = PMPI_Send(buf, cnt, datatype, dest, 
                                                tag + 100, comm);    
            }
        }
    }
    // PARALLEL PROTOCOL
    else
    {

    }
{
  //_wrap_py_return_val = PMPI_Send(buf, cnt, datatype, dest, tag, comm);
}
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Recv ================== */
_EXTERN_C_ int PMPI_Recv(void *buf, int cnt, MPI_Datatype datatype, int source,
                         int tag, MPI_Comm comm, MPI_Status *status);
_EXTERN_C_ int MPI_Recv(void *buf, int cnt, MPI_Datatype datatype, int source, 
                         int tag, MPI_Comm comm, MPI_Status *status) 
{
    int _wrap_py_return_val = 0;
    if(!isLeaderNode)
        source = source + fakeNumNodes;

    // MIRRORED PROTOCOL
    if(protoType == MIRRORED)
    {
        if(isLeaderNode)
        {
            
            // Send to next node first
            if(alive[source])
            {
                printf("\nRank#%d MPI_Recv %d->%d tag %d", trueRank, 
                    source, trueRank, tag);
                _wrap_py_return_val = PMPI_Recv(buf, cnt, datatype, source, 
                                                tag, comm, status);
            }

            //Get replica node rank and send if alive
            source = source + fakeNumNodes;
            
            if(alive[source])
            {
                printf("\nRank#%d MPI_Recv %d->%d tag %d", trueRank, source, 
                    trueRank, tag + 100);
                _wrap_py_return_val = PMPI_Recv(buf, cnt, datatype, source, 
                                                tag + 100, comm, status);
            }
        }
        else
        {
            // Send to next node first
            if(alive[source])
            {
                printf("\nRank#%d MPI_Recv %d->%d tag %d", trueRank, source, 
                    trueRank, tag + 100);

                _wrap_py_return_val = PMPI_Recv(buf, cnt, datatype, source, 
                                                tag + 100, comm, status);
            }

            //Get replica node rank and send if alive
            source = source - fakeNumNodes;
            
            if(alive[source])
            {
                printf("\nRank#%d MPI_Recv %d->%d tag %d", trueRank, source, 
                    trueRank, tag);
                _wrap_py_return_val = PMPI_Recv(buf, cnt, datatype, source, 
                                                tag, comm, status);    
            }
        }
    }
    // PARALLEL PROTOCOL
    else
    {

    }
{
//_wrap_py_return_val = PMPI_Recv(buf, cnt, datatype, source, tag, 
//                                      comm, status);
}
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Comm_size ================== */
_EXTERN_C_ int PMPI_Comm_size(MPI_Comm arg_0, int *arg_1);
_EXTERN_C_ int MPI_Comm_size(MPI_Comm arg_0, int *arg_1) { 
    int _wrap_py_return_val = 0;
 
    {
    //  _wrap_py_return_val = PMPI_Comm_size(arg_0, arg_1);
    }
    //return _wrap_py_return_val;
    *arg_1 = trueNumNodes / 2;
    return 0;
}

/* ================== C Wrappers for MPI_Comm_rank ================== */
_EXTERN_C_ int PMPI_Comm_rank(MPI_Comm arg_0, int *arg_1);
_EXTERN_C_ int MPI_Comm_rank(MPI_Comm arg_0, int *arg_1) { 
    int _wrap_py_return_val = 0;
 
    {
    //  _wrap_py_return_val = PMPI_Comm_rank(arg_0, arg_1);
    }   
    //return _wrap_py_return_val;
    *arg_1 = trueRank % (trueNumNodes / 2);
    return 0;
}

/* ================== C Wrappers for MPI_Init ================== */
_EXTERN_C_ int PMPI_Init(int *arg_0, char ***arg_1);
_EXTERN_C_ int MPI_Init(int *arg_0, char ***arg_1) { 
    int _wrap_py_return_val = 0; int i;
    
    char *bufProto = (char*)malloc(10);

    {
      _wrap_py_return_val = PMPI_Init(arg_0, arg_1);
    }

    PMPI_Comm_rank(MPI_COMM_WORLD, &trueRank);
    PMPI_Comm_size(MPI_COMM_WORLD, &trueNumNodes);

    MPI_Comm_rank(MPI_COMM_WORLD, &fakeRank);
    MPI_Comm_size(MPI_COMM_WORLD, &fakeNumNodes);

    if(trueRank == 0)
    {
        bufProto = getenv(PROTOCOL_TYPE);
        switch(strcmp(bufProto, "MIRRORED"))
        {
            case 0 : protoType = MIRRORED; break;
            case 1 : protoType = PARALLEL; break;
        }
    }

    alive = (int*) malloc (trueNumNodes * sizeof(int));
    scheduleDeath = (int*) malloc (trueNumNodes * sizeof(int));
    for(i = 0; i < trueNumNodes; ++i)
    {
        alive[i] = ALIVE;
        scheduleDeath[i] = ALIVE;
    }

    isLeaderNode = trueRank < fakeNumNodes ? 1 : 0;

    MPI_Bcast(&protoType, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Finalize ================== */
_EXTERN_C_ int PMPI_Finalize();
_EXTERN_C_ int MPI_Finalize() { 
    int _wrap_py_return_val = 0;
 
    printf("\n%d says bye bye\n", trueRank);
{
  _wrap_py_return_val = PMPI_Finalize();
}
    free(scheduleDeath);
    free(alive);
    return _wrap_py_return_val;
}