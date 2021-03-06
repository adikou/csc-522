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
#define _NUM_COLLECTIVES_ 6
#define _NUM_UNCOLLECTIVES_ 6

char* collectives[_NUM_COLLECTIVES_] = {"MPI_Barrier", "MPI_Scatter", 
                                        "MPI_Gather" , "MPI_Reduce", 
                                        "MPI_Alltoall","MPI_Allreduce"};

char* mpiOpNames[_NUM_MPI_OPS_] = {"MPI_Init", "MPI_Send", "MPI_Recv",
                                   "MPI_Isend", "MPI_Irecv", "MPI_Barrier",
                                   "MPI_Scatter", "MPI_Gather", "MPI_Reduce",
                                   "MPI_Allreduce", "MPI_Wait", "MPI_Waitall",
                                   "MPI_Alltoall", "MPI_Finalize"};

#define WHITE 0
#define GRAY  1
#define BLACK 2

typedef struct times
{
    int t;
    struct times *next;
}times;

times *mpiOpTimes[_NUM_MPI_OPS_];
int stats[_NUM_MPI_OPS_][5];

struct opCountNode
{
    int tag, count;
    struct opCountNode *next;
};

typedef struct opCount
{
    struct opCountNode *head;
}opCount;

opCount **opSeqCount;

typedef struct vertexMetaData
{
    char *key;
    char *interTreeChild;
} vertexMetaData;

typedef struct adjList
{
    int src;
    int dest;
    long int weight;
    long int bytes;
    int isCritical;
    struct adjList* next;
} adjList;
 
adjList *path = NULL, *criticalPath = NULL;

typedef struct graphVertex
{
    char *key;
    char *sendTarget;
    char *recvTarget;
    char *opName;
    int fromRank, toRank;
    int tag, opSeq;
    int inTreeWeight, interTreeWeight;
    long int bytes;
    int parent, Sv, color;
    int d, f;
    int id;
    adjList *root;
        
} graphVertex;

typedef struct graph
{
    int numGraphVertices;
    graphVertex *vertexListArray;
} graph; 

struct llist
{
    int vertex;
    struct llist *next;
}*root = NULL, *front = NULL, *rear = NULL;

struct requestList
{
    char *waitingOn;
    MPI_Request *request;
    int visited;
    struct requestList *next;
} *requestHead = NULL;

MPI_Request *curRequest = NULL; 

int myRank, numNodes;
FILE *fin, *fout;
char *baseFileName = "tmp", *fileName;

double startTime, endTime, discoveryTime, stime, etime;
double globalEndTime = 0;

int totalOps = 0, *numVertices, numCollectives = 0;
graphVertex **keys;
graph *Graph;
int numGraphVertices, *rankOffsetInMatrix, pathLength = 0;

long int *distTo;
adjList **edgeTo;

char *lastParentKey;

/* ================== End of Constants for MPI Operations ================== */

void insert(times **timeroot, int t)
{
    times *new = (times*)malloc(sizeof(times));
    new->t = t;
    new->next = *timeroot;
    *timeroot = new;
}

adjList* newAdjListNode(int dest, long int weight, long int bytes)
{
    adjList* new = (adjList*) malloc(sizeof(adjList));
    new->dest = dest;
    new->weight = weight;
    new->bytes = bytes;
    new->next = NULL;
    return new;
} 

graph* createGraph(int numGraphVertices)
{
    int i;

    graph* Graph = (graph*) malloc(sizeof(graph));
    Graph->numGraphVertices = numGraphVertices;
 
    Graph->vertexListArray = (graphVertex*) malloc(numGraphVertices * sizeof(graphVertex));
 
    for (i = 0; i < numGraphVertices; ++i)
    {
        Graph->vertexListArray[i].root = NULL;
        Graph->vertexListArray[i].key = Graph->vertexListArray[i].opName = NULL;
        Graph->vertexListArray[i].recvTarget = NULL;
        Graph->vertexListArray[i].id = i;
    }
 
    return Graph;
}
 
void addEdge(graph* Graph, int src, int dest, long int weight, long int bytes)
{
    adjList* new = newAdjListNode(dest, weight, bytes);
    new->src = src;
    new->isCritical = 0;
    new->next = Graph->vertexListArray[src].root;
    Graph->vertexListArray[src].root = new; 
}

void printGraph(graph *Graph)
{
    int i;
    for (i = 0; i < Graph->numGraphVertices; ++i)
    {
        adjList *cur = Graph->vertexListArray[i].root;
        printf("\nAdjacency list of vertex %d %s \nroot ", i, Graph->vertexListArray[i].key);
        while (cur)
        {
            printf("-> %d %s %ld", cur->dest, Graph->vertexListArray[cur->dest].key, 
                   cur->weight);
            cur = cur->next;
        }
        printf("\n");
    }
}
opCount **allocateOpSeqCount()
{
    int i;
    opCount *vals, **temp;

    // allocate values
    vals = (opCount *) malloc (numNodes * _NUM_MPI_OPS_ * sizeof(opCount));

    // allocate vector of pointers
    temp = (opCount **) malloc (_NUM_MPI_OPS_ * sizeof(opCount*));

    for(i=0; i < _NUM_MPI_OPS_; i++)
        temp[i] = &(vals[i * numNodes]);

    return temp;
}

void initOpSeqCount()
{
    int i, j;
    for(i = 0; i < _NUM_MPI_OPS_; i++)
        for(j = 0; j < numNodes; j++)
            opSeqCount[i][j].head = NULL;
}


char* getWaitingOnKey(MPI_Request *request)
{
    struct requestList *cur = requestHead;
    while(cur)
    {
        if(cur->request == request && cur->visited != 1)
        {
            cur->visited = 1;
            return cur->waitingOn;
        }
        cur = cur->next;
    }
}

void insertRequest(char *waitingOn, MPI_Request *request)
{
    struct requestList *new = (struct requestList*)malloc(sizeof(struct requestList));
    new->waitingOn = waitingOn;
    new->request = request;
    new->next = requestHead;
    requestHead = new;
}

int setAndGetCount(opCount *op, int tag)
{
    int retval;
    struct opCountNode *cur = op->head, *new;
    while(cur != NULL && cur->tag != tag)
        cur = cur->next;

    if(cur == NULL)
    {   
        new = (struct opCountNode*)malloc(sizeof(struct opCountNode));
        new->tag = tag;
        new->count = 1;
        retval = new->count;
        new->next = op->head;
        op->head = new;
    }
    else
    {
        cur->count = cur->count+1;
        retval = cur->count;
    }

    return retval;
}

int isCollective(char *key)
{
    int i, retval = 0;
    for(i = 0; i < _NUM_COLLECTIVES_; ++i)
    {
        if(strcmp(key, collectives[i])==0)
        {
            retval = 1;
            break;
        }
    } 
    return retval;
}


void countCollectives()
{
    int j;
    for (j = 0; j < numVertices[0]; ++j)
    {
        if(isCollective(keys[0][j].opName))
            numCollectives++;
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
        case _MPI_REDUCE_       : baseOp = "MPI_Reduce"; break;
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
                   long int bytes, double latency, int opSeq)
{
    vertexMetaData tmp;
    char *buf, *bufStr, *_buf;

    if(lastParentKey == NULL)
        lastParentKey = "NULL";

    switch(mpiOp)
    {
        case _MPI_INIT_ : 
                        fout = fopen(fileName, "w");
                        
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
        case _MPI_SEND_ : case _MPI_ISEND_ :
                        fout = fopen(fileName, "a");
                        

                        tmp.key = malloc(50);
                        tmp.interTreeChild = malloc(50);
                        generateKey(tmp.key, mpiOp, fromRank, 
                                    toRank, tag, opSeq);
                        lastParentKey = tmp.key;
                        fprintf(fout, "%s", tmp.key);
                        fprintf(fout, "%s", FOUR_SPACES);
                        buf = malloc(10);
                        sprintf(buf, "%d", (int)round((endTime - startTime)));
                        fprintf(fout, "%s", buf);
                        fprintf(fout, "%s", FOUR_SPACES);
                        generateKey(tmp.interTreeChild, _MPI_RECV_, fromRank, 
                                    toRank, tag, opSeq);
                        _buf = (char*)malloc(50);
                        strcpy(_buf, tmp.interTreeChild);
                        strcat(_buf, "||");
                        generateKey(tmp.interTreeChild, _MPI_IRECV_, fromRank, 
                                    toRank, tag, opSeq);
                        strcat(_buf, tmp.interTreeChild);
                        fprintf(fout, "%s", _buf);
                        fprintf(fout, "%s", FOUR_SPACES);
                        sprintf(buf, "%d", (int)round(latency));
                        fprintf(fout, "%s", buf);
                        fprintf(fout, "%s", FOUR_SPACES);
                        sprintf(buf, "%ld", bytes);
                        fprintf(fout, "%s", buf);
                        fprintf(fout, "%s", "\n");

                        fclose(fout);
                        break;
        case _MPI_IRECV_ : 
        case _MPI_RECV_ : 
                        fout = fopen(fileName, "a");

                        tmp.key = malloc(50);
                        generateKey(tmp.key, mpiOp, toRank, 
                                    fromRank, tag, opSeq);
                        fprintf(fout, "%s", tmp.key);
                        fprintf(fout, "%s", FOUR_SPACES);
                        
                        buf = malloc(10);
                        sprintf(buf, "%d", (int)round((endTime - startTime)));
                        fprintf(fout, "%s", buf);
                        fprintf(fout, "%s", "\n");

                        fclose(fout);
 
                        break;
                        
        case _MPI_BARRIER_ : case _MPI_SCATTER_ : 
        case _MPI_GATHER_  : case _MPI_ALLREDUCE_ :
        case _MPI_ALLTOALL_ : case _MPI_REDUCE_:
                        fout = fopen(fileName, "a");
                        tmp.key = malloc(50);
                        generateKey(tmp.key, mpiOp, fromRank, 
                                    toRank, tag, opSeq);
                        fprintf(fout, "%s", tmp.key);
                        fprintf(fout, "%s", FOUR_SPACES);   
                        buf = malloc(10);
                        sprintf(buf, "%d", (int)round((endTime - startTime)));
                        fprintf(fout, "%s", buf);
                        fprintf(fout, "%s", FOUR_SPACES);
                        fprintf(fout, "%s", "\n");
                        
                        fclose(fout);
                        break;
        
        case _MPI_WAIT_         : 
                        fout = fopen(fileName, "a");
                        tmp.key = malloc(50);
                        generateKey(tmp.key, mpiOp, fromRank, 
                                    toRank, tag, opSeq);
                        fprintf(fout, "%s", tmp.key);
                        fprintf(fout, "%s", FOUR_SPACES);   
                        buf = (char*)malloc(10);
                        sprintf(buf, "%d", (int)round((endTime - startTime)));
                        fprintf(fout, "%s", buf);
                        fprintf(fout, "%s", FOUR_SPACES);
                        _buf = getWaitingOnKey(curRequest);
                        fprintf(fout, "%s", _buf);
                        fprintf(fout, "%s", FOUR_SPACES);
                        fprintf(fout, "%s", "\n");
                        
                        fclose(fout);
                        break;
        case _MPI_WAITALL_      : break;
        
        case _MPI_FINALIZE_ : 
                        fout = fopen(fileName, "a");
                        //if(myRank == 0)
                        //{
                            tmp.key = malloc(30);
                            generateKey(tmp.key, mpiOp, fromRank, 
                                        toRank, tag, opSeq);

                            lastParentKey = tmp.key;
                            fprintf(fout, "%s", tmp.key);
                            fprintf(fout, "%s", FOUR_SPACES);
                            buf = malloc(10);
                            sprintf(buf, "%d", (int)round((endTime - startTime)));
                            fprintf(fout, "%s", buf);
                            fprintf(fout, "%s", "\n");
                        //}
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

void push(adjList *edge)
{
    adjList *new = (adjList*)malloc(sizeof(adjList));
    memcpy(new, edge, sizeof(adjList));
    new->next = path;
    path = new;
}

void insertVertex(int vertex)
{
    struct llist *new = (struct llist*)malloc(sizeof(struct llist));
    new->vertex = vertex;
    new->next = root;
    root = new;   
}

void DFS_visit(graph *Graph, graphVertex *u)
{
    int i = 0;
    graphVertex *v;
    discoveryTime++;
    u->d = discoveryTime;
    u->color = GRAY;

    adjList *cur = u->root;
    while(cur)
    {
        v = &(Graph->vertexListArray[cur->dest]);
        if(v->color == WHITE)
        {
            v->parent = u->id;
            DFS_visit(Graph,  v);
        }
        cur = cur->next;
    }

    u->color = BLACK;
    discoveryTime = discoveryTime + 1;
    u->f = discoveryTime;

    insertVertex(u->id);
}

void DFS(graph *Graph)
{
    int i;
    for(i = 0; i < Graph->numGraphVertices; ++i)
    {
        Graph->vertexListArray[i].color = WHITE;
        Graph->vertexListArray[i].parent = -1;
    }
    discoveryTime = 0;

    for(i = 0; i < Graph->numGraphVertices; ++i)
    {
        if(Graph->vertexListArray[i].color == WHITE)
            DFS_visit(Graph, &(Graph->vertexListArray[i]));
    }

}

void relax(graph* Graph, adjList *e)
{
    int v , w;
    v = e->src; w = e->dest;
    if(distTo[w] < distTo[v] + e->weight)
    {
        distTo[w] = distTo[v] + e->weight;
        edgeTo[w] = e;
    }
}

void longestPathInGraph(graph *Graph)
{
    int i;
    struct llist *cur = root;
    adjList *e; 
    
    distTo = (long int*)malloc(Graph->numGraphVertices * sizeof(long int));
    edgeTo = (struct adjList**)malloc(Graph->numGraphVertices * sizeof(struct adjList*));

    for(i = 0; i < Graph->numGraphVertices; ++i)
    {
        distTo[i] = (long int)-INFINITY;
        edgeTo[i] = NULL;
    }

    distTo[0] = 0;
    // cur is for topological ordering of G
    while(cur)
    {
        // e is for adjList of cur->vertex
        e = Graph->vertexListArray[cur->vertex].root;
        while(e)
        {
            relax(Graph, e);
            e = e->next;
        }
        //Done with this vertex
        cur = cur->next;

    }

    for(e = edgeTo[rankOffsetInMatrix[1]-1]; e != NULL; e = edgeTo[e->src])
        push(e);

    e  = path;
    while(e)
    {
        pathLength++; 
        e = e->next;
    }
    criticalPath = (adjList*)malloc(pathLength * sizeof(adjList));
    e = path; i = 0;
    while(e)
    {
        criticalPath[i].src = e->src;
        criticalPath[i].dest = e->dest;
        criticalPath[i].weight = e->weight;
        criticalPath[i].bytes = e->bytes;
        ++i;
        e = e->next;
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

void isEdgeCritical(adjList *e)
{
    int i, flag = 1;
    while(e && flag)
    {
        for(i = 0; i < pathLength; ++i)
            if(criticalPath[i].src == e->src &&
               criticalPath[i].dest == e->dest &&
               criticalPath[i].weight == e->weight)
            {
                if(flag)
                    e->isCritical = 1;
                flag = 0;
            }
        e = e->next;
    }
}

void fillDotGraph(graph *Graph)
{
    int i, j, boolExp;
    char *buf, *colorAttr = "color=\"RED\"";
    fout = fopen("dotGraph.txt", "w");

    buf = (char*)malloc(200);
    sprintf(buf, "digraph g{");
    fprintf(fout, "%s\n", buf);
    sprintf(buf, "overlap=scalexy");
    fprintf(fout, "%s\n", buf);
    sprintf(buf, "nodesep=0.6");
    fprintf(fout, "%s\n", buf);
    sprintf(buf, "node[fontsize=11]");
    fprintf(fout, "%s\n", buf);

    for(i = 0; i < numNodes; i++)
    {
        sprintf(buf, "subgraph cluster_%d{", i);
        fprintf(fout, "%s\n", buf);
        sprintf(buf, "label=\"Rank %d\"", i);
        fprintf(fout, "%s\n", buf);

        //Create rank-wise non-collective vertices
        for(j = rankOffsetInMatrix[i]; j < rankOffsetInMatrix[i+1]; ++j)
        {
            boolExp = (isCollective(Graph->vertexListArray[j].opName))
                || (strcmp(Graph->vertexListArray[j].opName, "MPI_Init")== 0)
                || (strcmp(Graph->vertexListArray[j].opName, "MPI_Finalize")== 0);
            if(!boolExp)
            {
                sprintf(buf, "%d [label=\"%s\"];",j, Graph->vertexListArray[j].opName);
                fprintf(fout, "%s\n", buf);
            }
        }

        sprintf(buf,"}");
        fprintf(fout, "%s\n", buf);
    }

    for(j = rankOffsetInMatrix[0]; j < rankOffsetInMatrix[1]; j++)
    {
        boolExp = (isCollective(Graph->vertexListArray[j].opName))
                || (strcmp(Graph->vertexListArray[j].opName, "MPI_Init")== 0)
                || (strcmp(Graph->vertexListArray[j].opName, "MPI_Finalize")== 0);
        if(boolExp)
        {
            sprintf(buf, "%d [label=%s];", j, Graph->vertexListArray[j].opName);
            fprintf(fout, "%s\n", buf);
        }
       
    }
    fprintf(fout, "\n");
    //Add edges
    for (i = 0; i < Graph->numGraphVertices; ++i)
    {
        adjList *cur = Graph->vertexListArray[i].root;
        while (cur)
        {
            isEdgeCritical(cur);
            if((strcmp(Graph->vertexListArray[cur->src].opName, "MPI_Send")== 0)
                || (strcmp(Graph->vertexListArray[cur->src].opName, "MPI_Isend")== 0))
            {   
                if(cur->next == NULL)
                {
                    if(cur->isCritical)
                        sprintf(buf, "%d->%d [label=\"%ld(%ld)\" %s]", cur->src, cur->dest, 
                            cur->weight, cur->bytes, colorAttr);
                    else
                        sprintf(buf, "%d->%d [label=\"%ld(%ld)\"]", cur->src, cur->dest, 
                            cur->weight, cur->bytes);
                    fprintf(fout, "%s\n", buf);
                }
                else
                {
                    if(cur->isCritical)
                        sprintf(buf, "%d->%d [label=\"%ld\" %s]", 
                                cur->src, cur->dest, cur->weight, colorAttr);
                    else
                        sprintf(buf, "%d->%d [label=\"%ld\"]", 
                                cur->src, cur->dest, cur->weight);
                    fprintf(fout, "%s\n", buf);
                }
            }
            else
            {
                if(cur->isCritical)
                    sprintf(buf, "%d->%d [label=\"%ld\" %s]", 
                            cur->src, cur->dest, cur->weight, colorAttr);
                else
                    sprintf(buf, "%d->%d [label=\"%ld\"]", 
                            cur->src, cur->dest, cur->weight);
                fprintf(fout, "%s\n", buf);
            }
            cur = cur->next;
        }
    }
    fprintf(fout, "}");
    fclose(fout);
}

void writeToCritPathOut(graph *Graph)
{
    int i;
    char *buf;
    FILE *fout;

    buf = (char*)malloc(100);
    fout = fopen("critPath.out", "w");
    for(i = 0; i < pathLength; ++i)
    {
        if(isCollective(Graph->vertexListArray[criticalPath[i].src].opName)
            || (strcmp(Graph->vertexListArray[criticalPath[i].src].opName, "MPI_Init")== 0))
        {
            sprintf(buf, "%s -1", Graph->vertexListArray[criticalPath[i].src].opName);
            fprintf(fout, "%s\n", buf);
            sprintf(buf, "%ld", criticalPath[i].weight);
            fprintf(fout, "%s\n", buf);
        }
        else if((strcmp(Graph->vertexListArray[criticalPath[i].src].opName, "MPI_Send")== 0)
                ||(strcmp(Graph->vertexListArray[criticalPath[i].src].opName, "MPI_Isend")== 0))
        {
            sprintf(buf, "%s %d", Graph->vertexListArray[criticalPath[i].src].opName,
                    Graph->vertexListArray[criticalPath[i].src].fromRank);
            fprintf(fout, "%s\n", buf);
            if((strcmp(Graph->vertexListArray[criticalPath[i].src].opName, "MPI_Recv")== 0)
                ||(strcmp(Graph->vertexListArray[criticalPath[i].src].opName, "MPI_Irecv")== 0)
                ||(strcmp(Graph->vertexListArray[criticalPath[i].src].opName, "MPI_Wait")== 0))
                sprintf(buf, "%ld", criticalPath[i].bytes);
            else sprintf(buf, "%ld", criticalPath[i].weight);
            fprintf(fout, "%s\n", buf);   
        }
        else if((strcmp(Graph->vertexListArray[criticalPath[i].src].opName, "MPI_Recv")== 0)
                ||(strcmp(Graph->vertexListArray[criticalPath[i].src].opName, "MPI_Irecv")== 0))
        {
            sprintf(buf, "%s %d", Graph->vertexListArray[criticalPath[i].src].opName,
                    Graph->vertexListArray[criticalPath[i].src].toRank);
            fprintf(fout, "%s\n", buf);
            sprintf(buf, "%ld", criticalPath[i].weight);
            fprintf(fout, "%s\n", buf);   
        }
        else
        {
            sprintf(buf, "%s %d", Graph->vertexListArray[criticalPath[i].src].opName,
                    Graph->vertexListArray[criticalPath[i].src].fromRank);
            fprintf(fout, "%s\n", buf);
            sprintf(buf, "%ld", criticalPath[i].weight);
            fprintf(fout, "%s\n", buf);       
        }
        if(strcmp(Graph->vertexListArray[criticalPath[i].dest].opName, "MPI_Finalize")== 0)
        {
            sprintf(buf, "%s -1", Graph->vertexListArray[criticalPath[i].dest].opName);
            fprintf(fout, "%s", buf);
        }
    }

    fclose(fout);
}

int partition(int* A, int p, int r)
{
    int i, x, j, tmp;
    x = A[r];
    i = p-1;
    for(j = p; j < r; j++)
    {
        if(A[j] <= x)
        {
            i = i + 1;
            tmp = A[i];
            A[i] = A[j];
            A[j] = tmp;
        }
    }
    tmp = A[i+1];
    A[i+1] = A[r];
    A[r] = tmp;

    return i+1;
}

void quicksort(int* A, int p, int r)
{
    int i,q;
    if(p < r)
    {
        q = partition(A, p, r);
        quicksort(A, p, q-1);
        quicksort(A, q+1, r);
    }
}
void writeToStatsDat()
{
    int i,j,k, count, sum;
    char *timeFile, *buf, *tmpFileName, *line = NULL, *bufLine = NULL;
    char *bufToken = NULL;
    ssize_t read; 
    size_t len = 0;
    int *timeVals[_NUM_MPI_OPS_], index[_NUM_MPI_OPS_];

    times *curt;
    buf = (char*)malloc(20);
    timeFile = (char*)malloc(strlen("time") + 5);
    strcpy(timeFile, "time");
    sprintf(buf, "%d", myRank); strncat(timeFile, buf, sizeof(buf));
    strncat(timeFile, ".txt", sizeof(".txt"));

    fout = fopen(timeFile, "w");

    for(k = 0; k < _NUM_MPI_OPS_; ++k)
    {
        curt = mpiOpTimes[k];
        sprintf(buf, "%s ", mpiOpNames[k]);
        fprintf(fout, "%s", buf);
        while(curt)
        {
            if(curt->next != NULL)
            {
                sprintf(buf, "%d ", curt->t);
                fprintf(fout, "%s", buf);
            }
            else
            {
                sprintf(buf, "%d", curt->t);
                fprintf(fout, "%s", buf);   
            }
            curt = curt->next;
        }
        
        if(k != _NUM_MPI_OPS_ - 1)
        {
            sprintf(buf, "%s", "\n");
            fprintf(fout, "%s", buf);
        }
    }

    fclose(fout);
    PMPI_Barrier(MPI_COMM_WORLD);
    if(myRank == 0)
    {
        fout = fopen("stats.dat", "w");
        for(i = 0; i < numNodes; i++)
        {
            tmpFileName = malloc(strlen("time") + 5);
            strcpy(tmpFileName, "time");
            sprintf(buf, "%d", i); strncat(tmpFileName, buf, sizeof(buf));
            strncat(tmpFileName, ".txt", sizeof(".txt"));      
            fin = fopen(tmpFileName, "r");
            for(k = 0; k < _NUM_MPI_OPS_; k++)
            {
                read = getline(&line, &len, fin);
                bufLine = (char*)malloc(strlen(line)+1);
                strcpy(bufLine, line);
                bufToken = strtok(bufLine, " \n");
                bufToken = strtok(NULL, " \n");
                while(bufToken)
                {
                    stats[k][0]++;
                    bufToken = strtok(NULL, " \n");
                }

                index[k] = 0;
                if(stats[k][0] > 0)
                    timeVals[k] = (int*)malloc(stats[k][0] * sizeof(int));
            }
            fclose(fin);
        }

        for(i = 0; i < numNodes; ++i)
        {
            tmpFileName = malloc(strlen("time") + 5);
            strcpy(tmpFileName, "time");
            sprintf(buf, "%d", i); strncat(tmpFileName, buf, sizeof(buf));
            strncat(tmpFileName, ".txt", sizeof(".txt"));      
            fin = fopen(tmpFileName, "r");
            for(k = 0; k < _NUM_MPI_OPS_; k++)
            {
                read = getline(&line, &len, fin);
                bufLine = (char*)malloc(strlen(line)+1);
                strcpy(bufLine, line);
                bufToken = strtok(bufLine, " \n");
                bufToken = strtok(NULL, " \n");
                while(bufToken)
                {
                    timeVals[k][index[k]] = atoi(bufToken);
                    bufToken = strtok(NULL, " \n");
                    ++index[k];
                }
            }
            fclose(fin);
        }

        for(i = 0; i < _NUM_MPI_OPS_; ++i)
            quicksort(timeVals[i], 0, stats[i][0]-1);
       
        for(i = 0 ; i < _NUM_MPI_OPS_; ++i)
            if(stats[i][0] > 0)
            {
                stats[i][2] = timeVals[i][0];
                stats[i][3] = timeVals[i][stats[i][0]/2];
                stats[i][4] = timeVals[i][stats[i][0]-1];
                sum = 0;
                for(j = 0; j < stats[i][0]; j++)
                    sum += timeVals[i][j];
                stats[i][1] = sum/stats[i][0];
            }
        
        sprintf(buf, "Function\tInvocations\tMean\tMin\tMedian\tMax\n");
        fprintf(fout, "%s", buf);

        for(i = 0; i < _NUM_MPI_OPS_-1; i++)
        {
            sprintf(buf, "%s\t", mpiOpNames[i]);
            fprintf(fout, "%s", buf);
            for(j = 0; j < 5; ++j)
            {
                sprintf(buf, "%d\t", stats[i][j]);
                fprintf(fout, "%s", buf);
            }
            sprintf(buf, "\n");
            fprintf(fout, "%s", buf);
        }
        fclose(fout);
    }

}

/* ================== C Wrappers for MPI_Barrier ================== */
_EXTERN_C_ int PMPI_Barrier(MPI_Comm arg_0);
_EXTERN_C_ int MPI_Barrier(MPI_Comm arg_0) 
{ 
    int _wrap_py_return_val = 0, count;
    totalOps++;
    
    
    endTime = MPI_Wtime();
    
  {
    stime = MPI_Wtime();
    _wrap_py_return_val = PMPI_Barrier(arg_0);
    etime = MPI_Wtime();
    insert(&mpiOpTimes[_MPI_BARRIER_], (int)round((etime - stime)));
  }
  count = setAndGetCount(&opSeqCount[_MPI_BARRIER_][myRank], -1);
    writeMetaData(_MPI_BARRIER_, 0, 0, -1, 0, 0, count); 
    startTime = MPI_Wtime();
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
    int _wrap_py_return_val = 0, count;
    totalOps++;
    
    endTime = MPI_Wtime();
    
  {
    stime = MPI_Wtime();
    _wrap_py_return_val = PMPI_Alltoall(arg_0, arg_1, arg_2, arg_3, arg_4, 
                                        arg_5, arg_6);
    etime = MPI_Wtime();
    insert(&mpiOpTimes[_MPI_ALLTOALL_], (int)round((etime - stime)));
  }
  count = setAndGetCount(&opSeqCount[_MPI_ALLTOALL_][myRank], -1);
    writeMetaData(_MPI_ALLTOALL_, 0, 0, -1, 0, 0, count); 
    startTime = MPI_Wtime();
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
    int _wrap_py_return_val = 0, count;
    totalOps++;
    endTime = MPI_Wtime();
    {
        stime = MPI_Wtime();
      _wrap_py_return_val = PMPI_Scatter(arg_0, arg_1, arg_2, arg_3, arg_4, 
                                         arg_5, arg_6, arg_7);
      etime = MPI_Wtime();
      insert(&mpiOpTimes[_MPI_SCATTER_], (int)round((etime - stime)));
    }
    count = setAndGetCount(&opSeqCount[_MPI_SCATTER_][myRank], -1);
    writeMetaData(_MPI_SCATTER_, 0, 0, -1, 0, 0, count);
    startTime = MPI_Wtime();
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
    int _wrap_py_return_val = 0, count;
    totalOps++;
    
    endTime = MPI_Wtime();
 
    {
        stime = MPI_Wtime();
      _wrap_py_return_val = PMPI_Gather(arg_0, arg_1, arg_2, arg_3, arg_4, 
                                        arg_5, arg_6, arg_7);
      etime = MPI_Wtime();
      insert(&mpiOpTimes[_MPI_GATHER_], (int)round((etime - stime)));
    }
    count = setAndGetCount(&opSeqCount[_MPI_GATHER_][myRank], -1);
    writeMetaData(_MPI_GATHER_, 0, 0, -1, 0, 0, count);
    startTime = MPI_Wtime();
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
    int _wrap_py_return_val = 0, count;
    totalOps++;
     endTime = MPI_Wtime();
    {
        stime = MPI_Wtime();
      _wrap_py_return_val = PMPI_Reduce(arg_0, arg_1, arg_2, arg_3, arg_4, 
                                        arg_5, arg_6);
      etime = MPI_Wtime();
      insert(&mpiOpTimes[_MPI_REDUCE_], (int)round((etime - stime)));
    }

    count = setAndGetCount(&opSeqCount[_MPI_REDUCE_][myRank], -1);
    writeMetaData(_MPI_REDUCE_, 0, 0, -1, 0, 0, count);
    startTime = MPI_Wtime();
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
    int _wrap_py_return_val = 0, count;
    totalOps++;
    endTime = MPI_Wtime();
    {
        stime = MPI_Wtime();
      _wrap_py_return_val = PMPI_Allreduce(arg_0, arg_1, arg_2, arg_3, 
                                           arg_4, arg_5);
      etime = MPI_Wtime();
      insert(&mpiOpTimes[_MPI_ALLREDUCE_], (int)round((etime - stime)));
    }

    count = setAndGetCount(&opSeqCount[_MPI_ALLREDUCE_][myRank], -1);
    writeMetaData(_MPI_ALLREDUCE_, 0, 0, -1, 0, 0, count);
    startTime = MPI_Wtime();
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Send ================== */
_EXTERN_C_ int PMPI_Send(void *buf, int cnt, MPI_Datatype datatype, int dest, 
                         int tag, MPI_Comm comm);
_EXTERN_C_ int MPI_Send(void *buf, int cnt, MPI_Datatype datatype, int dest, 
                         int tag, MPI_Comm comm) 
{ 
    int _wrap_py_return_val = 0, count;
    double latency;
    int dtypeSize;
    MPI_Type_size(datatype, &dtypeSize);
    latency = (0.000000291935 * cnt *sizeof(datatype)) + 0.000598493;
    totalOps++;
    endTime = MPI_Wtime();
    
    {
        stime = MPI_Wtime();
      _wrap_py_return_val = PMPI_Send(buf, cnt, datatype, dest, tag, comm);
      etime = MPI_Wtime();
      insert(&mpiOpTimes[_MPI_SEND_], (int)round((etime - stime)));
    }

    count = setAndGetCount(&opSeqCount[_MPI_SEND_][dest], tag);
    writeMetaData(_MPI_SEND_, myRank, dest, tag, cnt * dtypeSize, latency, count);
    startTime = MPI_Wtime();

    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Isend ================== */
_EXTERN_C_ int PMPI_Isend(void *buf, int cnt, MPI_Datatype datatype, int dest, 
                         int tag, MPI_Comm comm, MPI_Request *request);
_EXTERN_C_ int MPI_Isend(void *buf, int cnt, MPI_Datatype datatype, int dest, 
                         int tag, MPI_Comm comm, MPI_Request *request) 
{ 
    int _wrap_py_return_val = 0, count;
    double latency;
    int dtypeSize;
    char *tmp = NULL;
    MPI_Type_size(datatype, &dtypeSize);
    latency = (0.000000291935 * cnt *sizeof(datatype)) + 0.000598493;
    totalOps++;
    endTime = MPI_Wtime();
    {
        stime = MPI_Wtime();
      _wrap_py_return_val = PMPI_Isend(buf, cnt, datatype, dest, tag, comm, request);
      etime = MPI_Wtime();
      insert(&mpiOpTimes[_MPI_ISEND_], (int)round((etime - stime)));
    }

    tmp = (char*)malloc(50);
    count = setAndGetCount(&opSeqCount[_MPI_SEND_][dest], tag);
    writeMetaData(_MPI_ISEND_, myRank, dest, tag, cnt * dtypeSize, latency, count);
    generateKey(tmp, _MPI_ISEND_, myRank, dest, tag, count);
    insertRequest(tmp, request);
    startTime = MPI_Wtime();
    return _wrap_py_return_val;
}


/* ================== C Wrappers for MPI_Recv ================== */
_EXTERN_C_ int PMPI_Recv(void *buf, int cnt, MPI_Datatype datatype, int source, 
                         int tag, MPI_Comm comm, MPI_Status *status);
_EXTERN_C_ int MPI_Recv(void *buf, int cnt, MPI_Datatype datatype, int source, 
                         int tag, MPI_Comm comm, MPI_Status *status) 
{ 
    int _wrap_py_return_val = 0, count;
    totalOps++;
    endTime = MPI_Wtime();
     
    {
        stime = MPI_Wtime();
      _wrap_py_return_val = PMPI_Recv(buf, cnt, datatype, source, tag, 
                                      comm, status);
      etime = MPI_Wtime();
      insert(&mpiOpTimes[_MPI_RECV_], (int)round((etime - stime)));
    }

    count = setAndGetCount(&opSeqCount[_MPI_RECV_][source], tag);
    writeMetaData(_MPI_RECV_, myRank, source, tag, cnt, 0, count);
    startTime = MPI_Wtime();
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Irecv ================== */
_EXTERN_C_  int PMPI_Irecv(void *buf, int cnt, MPI_Datatype datatype, int source,
                             int tag, MPI_Comm comm, MPI_Request *request);

_EXTERN_C_  int MPI_Irecv(void *buf, int cnt, MPI_Datatype datatype, int source,
                             int tag, MPI_Comm comm, MPI_Request *request)

{ 
    int _wrap_py_return_val = 0, count;
    char *tmp;
    totalOps++;
    endTime = MPI_Wtime();
    {
        stime = MPI_Wtime();
      _wrap_py_return_val = PMPI_Irecv(buf, cnt, datatype, source, tag, 
                                      comm, request);
      etime = MPI_Wtime();
      insert(&mpiOpTimes[_MPI_IRECV_], (int)round((etime - stime)));
    }

    tmp = (char*)malloc(50);
    count = setAndGetCount(&opSeqCount[_MPI_RECV_][source], tag);
    writeMetaData(_MPI_IRECV_, myRank, source, tag, cnt, 0, count);
    generateKey(tmp, _MPI_IRECV_, source, myRank, tag, count);
    insertRequest(tmp, request);
    startTime = MPI_Wtime();
    return _wrap_py_return_val;
}
/* ================== C Wrappers for MPI_Wait ================== */
_EXTERN_C_ int PMPI_Wait(MPI_Request *request, MPI_Status *arg_1);
_EXTERN_C_ int MPI_Wait(MPI_Request *request, MPI_Status *arg_1) 
{ 
    int _wrap_py_return_val = 0, count;
    totalOps++;
    endTime = MPI_Wtime();
    curRequest = request;
    {
      stime = MPI_Wtime();
      _wrap_py_return_val = PMPI_Wait(request, arg_1);
      etime = MPI_Wtime();
      insert(&mpiOpTimes[_MPI_WAIT_], (int)round((etime - stime)));
    }

    count = setAndGetCount(&opSeqCount[_MPI_WAIT_][myRank], -1);
    writeMetaData(_MPI_WAIT_, myRank, 0, -1, 0, 0, count);
    startTime = MPI_Wtime();
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Waitall ================== */
_EXTERN_C_ int PMPI_Waitall(int reqCount, MPI_Request *request, MPI_Status *arg_2);
_EXTERN_C_ int MPI_Waitall(int reqCount, MPI_Request *request, MPI_Status *arg_2) 
{ 
    int _wrap_py_return_val = 0, i, count;

    for(i = 0; i < reqCount; ++i)
    {
        totalOps++;
        
        endTime = MPI_Wtime();
        curRequest = &request[i];
        stime = MPI_Wtime();
      _wrap_py_return_val = PMPI_Waitall(reqCount, request, arg_2);

      etime = MPI_Wtime();
      count = setAndGetCount(&opSeqCount[_MPI_WAIT_][myRank], -1);
        writeMetaData(_MPI_WAIT_, myRank, 0, -1, 0, 0, count);
        startTime = MPI_Wtime();
      insert(&mpiOpTimes[_MPI_WAIT_], (int)round((etime - stime)));
    }
    
   
    return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Init ================== */
_EXTERN_C_ int PMPI_Init(int *argc, char ***argv);
_EXTERN_C_ int MPI_Init(int *argc, char ***argv) 
{ 
    int _wrap_py_return_val = 0;
    char *buf;
    int i, j, count;

    fileName = malloc(strlen(baseFileName) + 4);
    buf = (char*)malloc(3);

    for(i = 0; i < _NUM_MPI_OPS_; ++i)
    {
        mpiOpTimes[i] = NULL;
        for(j = 0; j < 5; j++)
        {
                stats[i][j] = 0;
        }
    }

    strcpy(fileName, baseFileName);
    
      stime = MPI_Wtime();
      _wrap_py_return_val = PMPI_Init(argc, argv);
      etime = MPI_Wtime();
      insert(&mpiOpTimes[_MPI_INIT_], (int)round((etime - stime)));
    

    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); 
    MPI_Comm_size(MPI_COMM_WORLD, &numNodes);

    opSeqCount = allocateOpSeqCount();
    initOpSeqCount();

    sprintf(buf, "%d", myRank); strncat(fileName, buf, sizeof(buf));
    strncat(fileName, ".txt", 4);
    count = setAndGetCount(&opSeqCount[_MPI_INIT_][myRank], -1);
    startTime = MPI_Wtime();
    writeMetaData(_MPI_INIT_, 0, 0, -1, 0, 0, count);
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
    char *buf = malloc(20);
    char *tmpFileName;
    char ch, *bufKey, *line = NULL, *bufLine = NULL;
    char *starget, *tokens;
    ssize_t read; 
    size_t len = 0;
    int idInMatrix, u, count;
    adjList *v;
    int countTillEquals = 0, i, j, k, t, totalVerticesSoFar = 0;
    int firstNode, lastNode, ci, cj, flag = 10;
    FILE *fp; 
    graphVertex *vals = NULL;
    struct llist *cur;

    //if(myRank == 0)
        totalOps++;

    count = setAndGetCount(&opSeqCount[_MPI_FINALIZE_][myRank], -1);
    endTime = MPI_Wtime();
    //PMPI_Reduce(&endTime, &globalEndTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //endTime = globalEndTime;
    writeMetaData(_MPI_FINALIZE_, 0, 0, -1, 0, 0, count);    
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
                                        keys[i][j].recvTarget = 
                                    (char*)malloc(strlen(bufKey)+1);
                                        strcpy(keys[i][j].recvTarget, bufKey);
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
            //remove(tmpFileName);
        }

        countCollectives();


        numGraphVertices = numVertices[0];
        rankOffsetInMatrix[0] = 0;
        for (i = 1; i < numNodes; ++i)
        {
            numGraphVertices += numVertices[i] - numCollectives - 1; 
            rankOffsetInMatrix[i] = numGraphVertices - 
                                    (numVertices[i] - numCollectives) + 1;
        }
        rankOffsetInMatrix[numNodes] = numGraphVertices;
        
        Graph = createGraph(numGraphVertices);        
        j = 0;

        for(i = 0; i < numNodes; i++)
        {
            if(i == 0)
                lastNode = numVertices[i];
            else lastNode = numVertices[i] - 1;
            for(k = 0; k < lastNode; k++)
            {
                if(keys[i][k].opName != NULL)
                {
                    if(i == 0)
                    {
                        Graph->vertexListArray[j].key = (char*)malloc(strlen(keys[i][k].key)+1);
                        strcpy(Graph->vertexListArray[j].key, keys[i][k].key);
                        if(strcmp(keys[i][k].opName, "MPI_Wait")==0)
                        {
                            Graph->vertexListArray[j].recvTarget = (char*)malloc(strlen(keys[i][k].recvTarget)+1);
                            strcpy(Graph->vertexListArray[j].recvTarget, keys[i][k].recvTarget);
                        }
                        decomposeKey(&Graph->vertexListArray[j]);
                        j++;
                    }
                    else
                    {
                        if(!isCollective(keys[i][k].opName))
                        {
                            Graph->vertexListArray[j].key = (char*)malloc(strlen(keys[i][k].key)+1);
                            strcpy(Graph->vertexListArray[j].key, keys[i][k].key);
                            if(strcmp(keys[i][k].opName, "MPI_Wait")==0)
                            {
                                Graph->vertexListArray[j].recvTarget = (char*)malloc(strlen(keys[i][k].recvTarget)+1);
                                strcpy(Graph->vertexListArray[j].recvTarget, keys[i][k].recvTarget);
                            }
                            decomposeKey(&Graph->vertexListArray[j]);
                            j++;
                        }
                    }
                }
            }
            
        }

        for(i = 0; i< numNodes; ++i)
        {
            if(i == 0)
            {
                firstNode = 1;
                lastNode = numVertices[i] - 1;
            }
            else
            {
                firstNode = 0;
                lastNode = numVertices[i] - 2;
            }

            if(i)
            {
                if(firstNode>=0 && lastNode >=0 && 
                    (keys[i][firstNode].opName != NULL || keys[i][lastNode].opName != NULL))
                { 
                    if(isCollective(keys[i][firstNode].opName))
                        addEdge(Graph, 0, 1, keys[i][firstNode].inTreeWeight, 0);  
                    else addEdge(Graph, 0, idInGraph(Graph, keys[i][firstNode].key, i), 
                                 keys[i][firstNode].inTreeWeight, 0);
                    
                    if(isCollective(keys[i][lastNode].opName))
                        addEdge(Graph, idInGraph(Graph, keys[i][lastNode].key, 0), 
                                numVertices[0] - 1, 
                                keys[i][numVertices[i]-1].inTreeWeight, 0);  
                    else
                    { 
                        if((strcmp(keys[i][lastNode].opName,"MPI_Send")==0)
                            ||(strcmp(keys[i][lastNode].opName,"MPI_Isend")==0))
                        {
                            starget = (char*)malloc(50);
                            strcpy(starget, keys[i][lastNode].sendTarget);
                            tokens = strtok(starget, "|");
                            if(idInGraph(Graph, tokens,keys[i][lastNode].toRank) == -1)
                                tokens = strtok(NULL, "|");

                            else flag = 0;
                           
                            // Check for a wait to add to
                            for(ci = rankOffsetInMatrix[keys[i][lastNode].toRank]; 
                                ci < rankOffsetInMatrix[keys[i][lastNode].toRank + 1]; ++ci)
                            {
                                if(strcmp(Graph->vertexListArray[ci].opName, "MPI_Wait")==0)
                                {
                                    if(strcmp(Graph->vertexListArray[ci].recvTarget, tokens)==0)
                                    {
                                        flag = 0;
                                        addEdge(Graph, idInGraph(Graph,keys[i][lastNode].key,i), 
                                                    ci, keys[i][t].interTreeWeight, keys[i][lastNode].bytes);
                                    }
                                }
                            }
                            if(idInGraph(Graph, tokens,keys[i][t].toRank) == -1 && flag)
                                tokens = strtok(NULL, "|");
                            if(flag)
                                addEdge(Graph, idInGraph(Graph,keys[i][lastNode].key,i), 
                                        idInGraph(Graph, tokens, keys[i][lastNode].toRank), 
                                        keys[i][lastNode].interTreeWeight, keys[i][lastNode].bytes);
                        }
                        addEdge(Graph, rankOffsetInMatrix[i+1] - 1, 
                                numVertices[0] -1, keys[i][numVertices[i]-1].inTreeWeight, 0);
                    }
                }
            }
            for(t = 0; t < lastNode; t++)
            {
                if(keys[i][t].opName != NULL || keys[i][t+1].opName != NULL)
                {
                    // Current Single op Send/Isend - One Span to Recv also
                        if(isCollective(keys[i][t].opName))
                        { 
                            if(isCollective(keys[i][t+1].opName))
                            {
                                addEdge(Graph, idInGraph(Graph,keys[i][t].key,0),
                                    idInGraph(Graph,keys[i][t+1].key,0),
                                    keys[i][t+1].inTreeWeight, 0);
                            }
                            else addEdge(Graph, idInGraph(Graph,keys[i][t].key,0),
                                    idInGraph(Graph,keys[i][t+1].key,i),
                                    keys[i][t+1].inTreeWeight, 0);
                        }
                        else 
                        {
                            if((strcmp(keys[i][t].opName, "MPI_Send") == 0) 
                            || (strcmp(keys[i][t].opName, "MPI_Isend") == 0))
                            {
                                starget = (char*)malloc(50);
                                strcpy(starget, keys[i][t].sendTarget);
                                tokens = strtok(starget, "|");
                                if(idInGraph(Graph, tokens,keys[i][t].toRank) == -1)
                                    tokens = strtok(NULL, "|");

                                // Check for a wait to add to
                                for(ci = rankOffsetInMatrix[keys[i][t].toRank]; 
                                    ci < rankOffsetInMatrix[keys[i][t].toRank + 1]; ++ci)
                                {
                                    if(strcmp(Graph->vertexListArray[ci].opName, "MPI_Wait")==0)
                                    {
                                        if(strcmp(Graph->vertexListArray[ci].recvTarget, tokens)==0)
                                        {
                                            flag = 0;
                                            addEdge(Graph, idInGraph(Graph,keys[i][t].key,i), 
                                                    ci, keys[i][t].interTreeWeight, keys[i][t].bytes);
                                        }
                                    }
                                }
                                if(idInGraph(Graph, tokens,keys[i][t].toRank) == -1 && flag)
                                    tokens = strtok(NULL, "|");

                                if(flag)
                                    addEdge(Graph, idInGraph(Graph,keys[i][t].key,i), 
                                        idInGraph(Graph, tokens, keys[i][t].toRank), 
                                        keys[i][t].interTreeWeight, keys[i][t].bytes);
                            }

                            if(isCollective(keys[i][t+1].opName))
                                addEdge(Graph, idInGraph(Graph,keys[i][t].key,i),
                                        idInGraph(Graph,keys[i][t+1].key,0),
                                        keys[i][t+1].inTreeWeight, 0);
                            else addEdge(Graph, idInGraph(Graph,keys[i][t].key,i),
                                        idInGraph(Graph,keys[i][t+1].key,i),
                                        keys[i][t+1].inTreeWeight, 0);

                        }
                }
            }
        }  
        //printGraph(Graph);

        DFS(Graph);
        cur = root;
        longestPathInGraph(Graph);
        writeToCritPathOut(Graph);

        u = 0;
        cur = root;
        fillDotGraph(Graph);    
        k = system("mv dotGraph.txt dotGraph.dot");

    }   
    writeToStatsDat();
    _wrap_py_return_val = PMPI_Finalize();

    if(myRank == 0)
        k = system("dot -Tpng -oCritPathGraph.png dotGraph.dot");  

    return _wrap_py_return_val;
}