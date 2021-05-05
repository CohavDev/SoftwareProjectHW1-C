#include <stdio.h>
#include <stdlib.h>


 struct Cluster{
    int d;
    int count;
    //pointers to arrays
    double *sum;
    double *mean;

};
typedef struct Cluster Cluster;

struct Vector{
    double *coordinates;
//    struct cluster *clust;

};
typedef struct Vector Vector;

//input: Array of clusters and its length
//output: None, calcs means of all clusters
//TODO: didnt check correctness
int reCalcMeans(Cluster *clusterArray,int length){
    int calcMean(Cluster *clust);
    int changed = 0;
    int i;
    Cluster *t;
    t = clusterArray;
    for(i = 0;i<length;i++){
        if(calcMean(t) == 1){
            changed = 1;
        }
        t++;
    }
    return changed;
}
//calcs mean of specific cluster
//returns 1 if changed, 0 otherwise
int calcMean(Cluster *clust){
    int i;
    int changed = 0;
    double *sum, *mean; //pointers to sum and mean of clust
    sum = (*clust).sum;
    mean = (*clust).mean;
    for(i=0;i<(*clust).d;i++){
        if(*mean == (*sum / (*clust).count)){
            changed = 1;
        }
        *mean = *sum / (*clust).count;
        mean++;
        sum++;
    }
    return changed;
}
//adds specific vector to specific cluster
void addVector(Cluster *clust, Vector *vec){
    int i;
    double *sum, *coordinates;
    sum = (*clust).sum;
    coordinates = (*vec).coordinates;
    for(i = 0;i<(*clust).d;i++){
        *sum += *coordinates;
        sum++;
        coordinates++;
    }
}

//input: Array of clusters, it's length and a vector.
//output: None. the function find the closest cluster to the vector
// and insert vector to cluster.
void findCluster(Cluster *clusterArray, Vector *vec, int arrayLen){
    //declarations
    int i;
    Cluster *minCluster, *tempCluster;
    double minDistance,tempDistance;
    double distance(double *x,double *y, int length);
    void addVector(Cluster *clust, Vector *vec);

    //default minimum
    minCluster = clusterArray;
    tempCluster = minCluster;
    minDistance = distance((*vec).coordinates, (*minCluster).mean,(*minCluster).d);

    //loop through all clusters
    for(i = 1;i < arrayLen;i++){
        tempCluster++;
        tempDistance = distance((*vec).coordinates, (*tempCluster).mean, (*tempCluster).d);
        if(tempDistance < minDistance){
            minDistance = tempDistance;
            minCluster = tempCluster;
        }
    }
    addVector(minCluster, vec);
}
int main() {
    //test
//    double *p, *g, *t;
//    int i;
//    p = malloc(3*sizeof(double));
//    t = p;
//    for (i = 0;i < 3;i++){
//        *t = i;
//        t++;
//
//    }
//    g = malloc(3*sizeof(double));
//    t = g;
//    for (i = 0;i < 3;i++){
//        *t = i*2;
//        t++;
//
//    }
//    double distance(double *x,double *y, int length);
//    double d = distance(p, g, 3);
//    printf("Distance is: %.2f",d);
//    free(p);
//    free(g);
    return 0;
}

//kMeans function
void kMeans(int k, int maxIter ){
    void printMeans(Cluster *clusterArray);
    void initFromFile(int k, Cluster *clusterArray, Vector *vecArray);
    void findCluster(Cluster *clusterArray, Vector *vec, int arrayLen);
    void refreshClusters(Cluster *clusterArray, int clusterLength, int d);

    int reCalcMeans(Cluster *clusterArray,int length);

    Cluster *clusterArray = NULL;
    Vector *vecArray = NULL;
    initFromFile(k,clusterArray, vecArray);
    int lengthVectors; // number of vectors
    int lengthClusters; //number of clusters
    int i, iterCount;
    Vector *currentVec;
    currentVec = vecArray;
    int changed = 1;
    iterCount = 0;
    while(iterCount < maxIter && changed == 1){
        //loop through vectors and insert to clusters
        for(i=0;i<lengthVectors;i++){
            findCluster(clusterArray, currentVec, lengthClusters);
        }
        //recalcmeans returns 1 if mean changed, 0 otherwise
        changed = reCalcMeans(clusterArray, lengthClusters);
        //init sum and count to zero's:
        refreshClusters(clusterArray, lengthClusters, (*clusterArray).d);
        iterCount++;
    }
    printMeans(clusterArray);
    //TODO: continue from here
}
//functions

//input: 2 arrays(vectors coordinates)
//output: squered distance between them
//TODO: this works! :)
double distance(double *x, double *y, int length){
    double *a,*b;
    double sum = 0;
    int i;
    a = x;
    b = y;
    for(i=0; i<length;i++){
        sum += (*a-*b)*(*a-*b);
        a++;
        b++;
    }
    return sum;
}

//input: Array of clusters and it's length.number of coordinates.
//output: None. All clusters sum and count is initalized to zero
void refreshClusters(Cluster *clusterArray, int clusterLength, int d){
    int i,j;
    Cluster *temp;
    double *indexer;
    temp = clusterArray;
    for(i=0;i<clusterLength;i++){
        (*temp).count = 0;
        indexer = ((*temp).sum);
        //loop through array of sum and init all to zero (0)
        for(j = 0;j < d;j++){
            *indexer = 0;
            indexer++;
        }
        temp++;
    }
}
void initFromFile(int k, Cluster *clusterArray, Vector *vecArray){
    //TODO:this should initialize vectors and clusters arrays
}

void printMeans(Cluster *clusterArray){
    //TODO:this should print all cluster's means to cmd
}

