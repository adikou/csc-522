#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

#define PLUS "+"
// White space characters fail at tokenizer 
#define FOUR_SPACES " "

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
#define _NUM_MPI_OPS_   14

typedef struct vertexMetaData
{
    char *key;
    char *interTreeChild;
} vertexMetaData;

typedef struct adjList
{
    int dest;
    int weight;
    struct adjList* next;
} adjList;
 
 
adjList* newAdjListNode(int dest, int weight)
{
    adjList* new = (adjList*) malloc(sizeof(adjList));
    new->dest = dest;
    new->weight = weight;
    new->next = NULL;
    return new;
} 

typedef struct graphVertex
{
    char *key;
    char *sendTarget;
    char *opName;
    int fromRank, toRank;
    int tag, opSeq;
    int inTreeWeight, interTreeWeight;
    int bytes;
    adjList *root;
        
} graphVertex;

typedef struct graph
{
    int numGraphVertices;
    graphVertex *vertexListArray;
} graph; 

graph* createGraph(int numGraphVertices)
{
    int i;

    graph* Graph = (graph*) malloc(sizeof(graph));
    Graph->numGraphVertices = numGraphVertices;
 
    Graph->vertexListArray = (graphVertex*) malloc(numGraphVertices * sizeof(graphVertex));
 
    for (i = 0; i < numGraphVertices; ++i)
        Graph->vertexListArray[i].root = NULL;
 
    return Graph;
}
 
void addEdge(graph* Graph, int src, int dest, int weight)
{
    adjList* new = newAdjListNode(dest, weight);
    new->next = Graph->vertexListArray[src].root;
    Graph->vertexListArray[src].root = new; 
}

void printGraph(graph *Graph)
{
    int i;
    for (i = 0; i < Graph->numGraphVertices; ++i)
    {
        adjList *cur = Graph->vertexListArray[i].root;
        printf("\nAdjacency list of vertex %d %s\nroot ", i, Graph->vertexListArray[i].key);
        while (cur)
        {
            printf("-> %s %d", Graph->vertexListArray[cur->dest].key, 
                   cur->weight);
            cur = cur->next;
        }
        printf("\n");
    }
}
 
int myRank, numNodes;
FILE *fin, *fout;
char *baseFileName = "tmp", *fileName;

double startTime, endTime;

int **opSeqCount, totalOps = 0, *numVertices, numCollectives = 0;
graphVertex **keys;
graph *Graph;
int numGraphVertices, *rankOffsetInMatrix;

char *lastParentKey;

/* ================== End of Constants for MPI Operations ================== */

int **allocateOpSeqCount()
{
    int i;
    int *vals, **temp;

    // allocate values
    vals = (int *) malloc (numNodes * _NUM_MPI_OPS_ * sizeof(int));

    // allocate vector of pointers
    temp = (int **) malloc (_NUM_MPI_OPS_ * sizeof(int*));

    for(i=0; i < _NUM_MPI_OPS_; i++)
        temp[i] = &(vals[i * numNodes]);

    return temp;
}

void initOpSeqCount()
{
    int i, j;
    for(i = 0; i < _NUM_MPI_OPS_; i++)
        for(j = 0; j < numNodes; j++)
            opSeqCount[i][j] = 0;
}

void countCollectives()
{
    int j;
    for (j = 0; j < numVertices[0]; ++j)
    {
        switch(strcmp(keys[0][j].opName, "MPI_Barrier"))
        {
            case 0 : numCollectives++; break;
        }
    }
}

void decomposeKey(graphVertex *vertex)
{
    char *buf, i = 0;
    char *string = (char*)malloc(strlen(vertex->key) + 1);
    strcpy(string, vertex->key);
    buf = strtok(string,"+");
    while(buf)
    {
        switch(i)
        {
            case 0 : vertex->opName = (char*)malloc(strlen(buf) + 1);
                     strcpy(vertex->opName, buf);
                     break;
            case 1 : vertex->fromRank = atoi(buf);
                     break;
            case 2 : vertex->toRank = atoi(buf);
                     break;
            case 3 : vertex->tag = atoi(buf);
                     break;
            case 4 : vertex->opSeq = atoi(buf);
                     break;
        }
        buf = strtok(NULL,"+");
        i++;
    }
}

void generateKey(char *key, int op, int fromRank, int toRank, 
                 int tag, int opSeq)
{
    char *baseOp, *buf;
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

    strcpy(key, baseOp); strncat(key, PLUS, sizeof(PLUS));
    buf = malloc(2);
    sprintf(buf, "%d", fromRank); strncat(key, buf, sizeof(buf)); 
    strncat(key, PLUS, sizeof(PLUS));    
    sprintf(buf, "%d", toRank); strncat(key, buf, sizeof(buf)); 
    strncat(key, PLUS, sizeof(PLUS));
    sprintf(buf, "%d", tag); strncat(key, buf, sizeof(buf)); 
    strncat(key, PLUS, sizeof(PLUS));
    sprintf(buf, "%d", opSeq); strncat(key, buf, sizeof(buf));
}

void writeMetaData(int mpiOp, int fromRank, int toRank, int tag, 
                   int bytes, int opSeq)
{
    vertexMetaData tmp;
    char *buf, *bufStr;

    if(lastParentKey == NULL)
        lastParentKey = "NULL";

    switch(mpiOp)
    {
        case _MPI_INIT_ : 
                        fout = fopen(fileName, "w");
                        startTime = MPI_Wtime();
                        if(myRank == 0)
                        {
                            tmp.key = malloc(30);
                            generateKey(tmp.key, mpiOp, fromRank, toRank, 
                                        tag, opSeq);
                            fprintf(fout, "%s", tmp.key);
                            fprintf(fout, "%s", FOUR_SPACES);
                            buf = malloc(10);
                            sprintf(buf, "%d", 0);
                            fprintf(fout, "%s", buf);
                            fprintf(fout, "%s","\n");
                        }
                        fclose(fout);

                        break;
        case _MPI_SEND_ : 
                        fout = fopen(fileName, "a");
                        endTime = MPI_Wtime();

                        tmp.key = malloc(50);
                        tmp.interTreeChild = malloc(50);
                        generateKey(tmp.key, mpiOp, fromRank, 
                                    toRank, tag, opSeq);
                        lastParentKey = tmp.key;
                        fprintf(fout, "%s", tmp.key);
                        fprintf(fout, "%s", FOUR_SPACES);
                        buf = malloc(10);
                        sprintf(buf, "%d", (int)round((endTime - startTime)*1000)+1);
                        fprintf(fout, "%s", buf);
                        fprintf(fout, "%s", FOUR_SPACES);
                        generateKey(tmp.interTreeChild, _MPI_RECV_, fromRank, 
                                    toRank, tag, opSeq);
                        fprintf(fout, "%s", tmp.interTreeChild);
                        fprintf(fout, "%s", FOUR_SPACES);
                        sprintf(buf, "%d", 0);
                        fprintf(fout, "%s", buf);
                        fprintf(fout, "%s", FOUR_SPACES);
                        sprintf(buf, "%d", bytes);
                        fprintf(fout, "%s", buf);
                        fprintf(fout, "%s", "\n");

                        fclose(fout);
                        startTime = MPI_Wtime();
                        break;
        case _MPI_ISEND_ : break;
        case _MPI_RECV_ : 
                        fout = fopen(fileName, "a");
                        endTime = MPI_Wtime();

                        tmp.key = malloc(50);
                        generateKey(tmp.key, mpiOp, toRank, 
                                    fromRank, tag, opSeq);
                        fprintf(fout, "%s", tmp.key);
                        fprintf(fout, "%s", FOUR_SPACES);
                        
                        buf = malloc(10);
                        sprintf(buf, "%d", (int)round((endTime - startTime)*1000)+1);
                        fprintf(fout, "%s", buf);
                        fprintf(fout, "%s", "\n");

                        fclose(fout);
                        startTime = MPI_Wtime();
 
                        break;
        case _MPI_IRECV_ : 
                        fout = fopen(fileName, "a");
                        endTime = MPI_Wtime();
                        
                        tmp.key = malloc(50);
                        generateKey(tmp.key, mpiOp, toRank, 
                                    fromRank, tag, opSeq);
                        fprintf(fout, "%s", tmp.key);
                        fprintf(fout, "%s", FOUR_SPACES);
                        
                        buf = malloc(10);
                        sprintf(buf, "%d", (int)round((endTime - startTime)*1000)+1);
                        fprintf(fout, "%s", buf);
                        fprintf(fout, "%s", "\n");

                        fclose(fout);
                        startTime = MPI_Wtime();
                        break;
        case _MPI_BARRIER_ : 
                        fout = fopen(fileName, "a");
                        endTime = MPI_Wtime();
                        tmp.key = malloc(50);
                        generateKey(tmp.key, mpiOp, fromRank, 
                                    toRank, tag, opSeq);
                        fprintf(fout, "%s", tmp.key);
                        fprintf(fout, "%s", FOUR_SPACES);   
                        buf = malloc(10);
                        sprintf(buf, "%d", (int)round((endTime - startTime)*1000)+1);
                        fprintf(fout, "%s", buf);
                        fprintf(fout, "%s", FOUR_SPACES);
                        fprintf(fout, "%s", "\n");
                        
                        fclose(fout);
                        startTime = MPI_Wtime();
                        break;
        case _MPI_SCATTER_      : break;
        case _MPI_GATHER_       : break;
        case _MPI_ALLREDUCE_    : break;
        case _MPI_WAIT_         : break;
        case _MPI_WAITALL_      : break;
        case _MPI_ALLTOALL_     : break;
        case _MPI_FINALIZE_ : 
                        fout = fopen(fileName, "a");
                        if(myRank == 0)
                        {
                            endTime = MPI_Wtime();
                            tmp.key = malloc(30);
                            generateKey(tmp.key, mpiOp, fromRank, 
                                        toRank, tag, opSeq);

                            lastParentKey = tmp.key;
                            fprintf(fout, "%s", tmp.key);
                            fprintf(fout, "%s", FOUR_SPACES);
                            buf = malloc(10);
                            sprintf(buf, "%d", (int)round((endTime - startTime)*1000)+1);
                            fprintf(fout, "%s", buf);
                            fprintf(fout, "%s", "\n");
                        }
                        buf = malloc(10); 
                        sprintf(buf, "%d", totalOps);
                        bufStr = malloc(strlen("totalOps=") + strlen(buf) + 1);
                        bufStr = "totalOps=";
                        fprintf(fout, "%s", bufStr);
                        fprintf(fout, "%s", buf);
                        fclose(fout);
                        break;
    }

}

int idInGraph(graph *Graph, char* key, int whichRank)
{
    int i, flag = -1;
    for(i = rankOffsetInMatrix[whichRank]; i < rankOffsetInMatrix[whichRank+1]; ++i)
    {
        if(strcmp(Graph->vertexListArray[i].key, key)==0)
        {
            flag = i;
            break;
        }
    }
    return flag ;
}

int isCollective(char *key)
{
    if(strcmp(key, "MPI_Barrier")==0)
        return 1;
}


/* ================== C Wrappers for MPI_Barrier ================== */
_EXTERN_C_ int PMPI_Barrier(MPI_Comm arg_0);
_EXTERN_C_ int MPI_Barrier(MPI_Comm arg_0) 
{ 
    int _wrap_py_return_val = 0;
    totalOps++;
    
    writeMetaData(_MPI_BARRIER_, 0, 0, -1, 0, 
                  ++opSeqCount[_MPI_BARRIER_][myRank]); 
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
_EXTERN_C_ int PMPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, 
                         int tag, MPI_Comm comm);
_EXTERN_C_ int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, 
                         int tag, MPI_Comm comm) 
{ 
    int _wrap_py_return_val = 0;
    totalOps++;
    {
      _wrap_py_return_val = PMPI_Send(buf, count, datatype, dest, tag, comm);
    }
    writeMetaData(_MPI_SEND_, myRank, dest, tag, count, 
                  ++opSeqCount[_MPI_SEND_][dest]);
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
_EXTERN_C_ int PMPI_Recv(void *buf, int count, MPI_Datatype datatype, int dest, 
                         int tag, MPI_Comm comm, MPI_Status *status);
_EXTERN_C_ int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int dest, 
                         int tag, MPI_Comm comm, MPI_Status *status) 
{ 
    int _wrap_py_return_val = 0;
    totalOps++;
    writeMetaData(_MPI_RECV_, myRank, dest, tag, count, 
                  ++opSeqCount[_MPI_RECV_][dest]);
    {
      _wrap_py_return_val = PMPI_Recv(buf, count, datatype, dest, tag, 
                                      comm, status);
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
_EXTERN_C_ int PMPI_Init(int *argc, char ***argv);
_EXTERN_C_ int MPI_Init(int *argc, char ***argv) 
{ 
    int _wrap_py_return_val = 0;
    char *buf;
    int i, j;

    fileName = malloc(strlen(baseFileName) + 4);
    buf = malloc(2);

    strcpy(fileName, baseFileName);
    opSeqCount = allocateOpSeqCount();
    initOpSeqCount();
    {
      _wrap_py_return_val = PMPI_Init(argc, argv);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); 
    MPI_Comm_size(MPI_COMM_WORLD, &numNodes);

    sprintf(buf, "%d", myRank); strncat(fileName, buf, sizeof(buf));
    strncat(fileName, ".txt", 4);
    writeMetaData(_MPI_INIT_, myRank, myRank, -1, 0, 0);

    if(myRank == 0)
    {
        totalOps++;
        numVertices = (int*) malloc(numNodes * sizeof(int));
        rankOffsetInMatrix = (int*) malloc(numNodes + 1 * sizeof(int));
        keys = (graphVertex**) malloc (numNodes * sizeof(graphVertex*));
    }

    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Finalize ================== */
_EXTERN_C_ int PMPI_Finalize();
_EXTERN_C_ int MPI_Finalize() 
{ 
    int _wrap_py_return_val = 0;
    char *buf = malloc(5), *tmpFileName;
    char ch, *bufKey, *line = NULL, *bufLine = NULL;
    ssize_t read; 
    size_t len = 0;
    int idInMatrix;
    int countTillEquals = 0, i, j, k, t, totalVerticesSoFar = 0;
    int firstNode, lastNode;
    FILE *fp; 
    graphVertex *vals = NULL;

    if(myRank == 0)
        totalOps++;

    writeMetaData(_MPI_FINALIZE_, myRank, myRank, -1, 0, 0);
    PMPI_Barrier(MPI_COMM_WORLD);
    bufKey = (char*)malloc(50);
    if(myRank == 0)
    {
        for(i = 0; i < numNodes; i++)
        {
            tmpFileName = malloc(strlen("tmp") + 5);
            strcpy(tmpFileName, "tmp");
            sprintf(buf, "%d", i); strncat(tmpFileName, buf, sizeof(buf));
            strncat(tmpFileName, ".txt", sizeof(".txt"));      

            fin = fopen(tmpFileName, "r");
            fseek(fin, -1, SEEK_END);
            countTillEquals = 0;
            ch = fgetc(fin);
            while(ch != '=')
            {
                fseek(fin, -2, SEEK_CUR);
                ch = fgetc(fin);
                countTillEquals++;
            }

            for(j = 0; j < countTillEquals; j++)
                buf[j] = fgetc(fin);
            buf[j] = '\0';
            
            numVertices[i] = atoi(buf);
            totalVerticesSoFar += numVertices[i];
            keys[i] = (graphVertex*)realloc(vals, numVertices[i] * sizeof(graphVertex));

            fclose(fin);
            fin = fopen(tmpFileName, "r");
            
            for(j = 0; j < numVertices[i]; j++)
            {
                keys[i][j].key = NULL;
                read = getline(&line, &len, fin);
                bufLine = (char*)malloc(strlen(line)+1);
                strcpy(bufLine, line);
                if(read != -1)
                {
                    bufKey = strtok(line," \n");
                    keys[i][j].key = (char*)malloc(strlen(bufKey) + 1);
                    strcpy(keys[i][j].key, bufKey);

                    k = -1;
                    while(bufKey)
                    {
                        bufKey = strtok (NULL," \n");
                        if(bufKey)
                            if(k == -1)
                                keys[i][j].inTreeWeight = atoi(bufKey);
                        if(bufKey)
                        {
                            switch(k)
                            {
                                case 0: keys[i][j].sendTarget = 
                                    (char*)malloc(strlen(bufKey)+1);
                                        strcpy(keys[i][j].sendTarget, bufKey);
                                        break;
                                case 1: keys[i][j].interTreeWeight = atoi(bufKey);
                                        break;
                                case 2: keys[i][j].bytes = atoi(bufKey);
                                        break;
                            }
                        }
                        k++;
                    }
                    decomposeKey(&keys[i][j]);
                }
            }
            fclose(fin);
        }
        countCollectives();

        numGraphVertices = numVertices[0];
        rankOffsetInMatrix[0] = 0;
        for (i = 1; i < numNodes; ++i)
        {
            numGraphVertices += numVertices[i] - numCollectives; 
            rankOffsetInMatrix[i] = numGraphVertices - 
                                    (numVertices[i] - numCollectives);
        }
        rankOffsetInMatrix[numNodes] = numGraphVertices;

        Graph = createGraph(numGraphVertices);        
        j = 0;

        for(i = 0; i < numNodes; i++)
        {
            //Link all Init nodes to span-n ranks and n-ranks to Finalize
            //if(i)
            //addEdge(Graph, rankOffsetInMatrix[i+1] - 1, numVertices[0] - 1, 1);

            //Connect all edges in one thread of one rank
            for(k = 0; k < numVertices[i]; k++)
            {
                if(i == 0)
                {
                    Graph->vertexListArray[j].key = (char*)malloc(strlen(keys[i][k].key));
                    strcpy(Graph->vertexListArray[j].key, keys[i][k].key);
                    j++;
                }
                else
                {
                    if(strcmp(keys[i][k].opName, "MPI_Barrier") != 0)
                    {
                        Graph->vertexListArray[j].key = (char*)malloc(strlen(keys[i][k].key));
                        strcpy(Graph->vertexListArray[j].key, keys[i][k].key);
                        j++;

                    }
                }
                
            }
            
        }
        for(i = 0; i< numNodes; ++i)
        {
            for(t = 0; t < numVertices[i]-1; t++)
            {
                // Current Single op Send/Isend - One Span to Recv also
                
                    if(strcmp(keys[i][t].opName, "MPI_Send") == 0)
                    {
                        addEdge(Graph, idInGraph(Graph,keys[i][t].key,i), 
                                idInGraph(Graph, keys[i][t].sendTarget,keys[i][t].toRank), 
                                keys[i][t].interTreeWeight);
                        if(strcmp(keys[i][t+1].opName, "MPI_Barrier")==0)
                        {
                            addEdge(Graph, idInGraph(Graph,keys[i][t].key,i),
                                idInGraph(Graph,keys[i][t+1].key,keys[i][t+1].fromRank),
                                keys[i][t+1].inTreeWeight);
                        }
                        else addEdge(Graph, idInGraph(Graph,keys[i][t].key,i),
                                idInGraph(Graph,keys[i][t+1].key,i),
                                keys[i][t+1].inTreeWeight);
                    }

                // Current Collective op - Span N
                    else if(strcmp(keys[i][t].opName, "MPI_Barrier") == 0)
                    { 
                        if(strcmp(keys[i][t+1].opName, "MPI_Barrier")==0)
                        {
                            addEdge(Graph, idInGraph(Graph,keys[i][t].key,i),
                                idInGraph(Graph,keys[i][t+1].key,keys[i][t+1].fromRank),
                                keys[i][t+1].inTreeWeight);
                        }
                        else addEdge(Graph, idInGraph(Graph,keys[i][t].key,0),
                                idInGraph(Graph,keys[i][t+1].key,i),
                                keys[i][t+1].inTreeWeight);
                    }
                    // Non-Send single op - fans out to 1/N numNodes
                    else
                    {
                        if(strcmp(keys[i][t+1].opName, "MPI_Barrier")==0)
                        {
                            addEdge(Graph, idInGraph(Graph,keys[i][t].key,i),
                                    idInGraph(Graph,keys[i][t+1].key,keys[i][t+1].fromRank),
                                    keys[i][t+1].inTreeWeight);
                        }
                        else addEdge(Graph, idInGraph(Graph,keys[i][t].key,i),
                                    idInGraph(Graph,keys[i][t+1].key,i),
                                    keys[i][t+1].inTreeWeight);
                    }
            }
            if(i)
            {
                if(strcmp(keys[i][0].opName, "MPI_Barrier")==0)
                  addEdge(Graph, 0, 1, keys[i][0].inTreeWeight);  
                else addEdge(Graph, 0, idInGraph(Graph, keys[i][0].key, i), keys[i][0].inTreeWeight);
                if(strcmp(keys[i][t].opName, "MPI_Barrier")==0)
                  addEdge(Graph, numVertices[0]-2, numVertices[0] - 1, 
                          keys[0][numVertices[0]-1].inTreeWeight);  
                else addEdge(Graph, rankOffsetInMatrix[i+1] - 1, 
                            numVertices[0] -1, keys[0][numVertices[0]-1].inTreeWeight);
                addEdge(Graph, rankOffsetInMatrix[i+1] - 1, numVertices[0] - 1, 1);
            }
        }        

    }   
    _wrap_py_return_val = PMPI_Finalize();
    
    return _wrap_py_return_val;
}