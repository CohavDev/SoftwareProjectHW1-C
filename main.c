#include <stdio.h>
#include <stdlib.h>


struct Cluster {
    int d;
    int count;
    //pointers to arrays
    double *sum;
    double *mean;

};
typedef struct Cluster Cluster;

struct Vector {
    double *coordinates;
//    struct cluster *clust;

};
typedef struct Vector Vector;

//input: Array of clusters and its length
//output: None, calcs means of all clusters
//TODO: didnt check correctness
int reCalcMeans(Cluster *clusterArray, int length) {
    int calcMean(Cluster *clust);
    int changed = 0; // "False"
    int i;
    for (i = 0; i < length; i++) {
        if (calcMean(&(clusterArray[i])) == 1) {
            changed = 1;
        }
    }
    return changed;
}

//calcs mean of specific cluster
//returns 1 if changed, 0 otherwise
int calcMean(Cluster *clust) {
    int i;
    int changed = 0; //"False"
    double *sum, *mean; //pointers to sum and mean of clust
    double calcVal;
    sum = (*clust).sum;
    mean = (*clust).mean;
    for (i = 0; i < (*clust).d; i++) {
        calcVal = (sum[i] / (*clust).count);
        if (mean[i] != calcVal) {
            changed = 1; //"True"
        }
        mean[i] = calcVal;
    }
    return changed;
}

//adds specific vector to specific cluster
void addVector(Cluster *clust, Vector *vec) {
    int i;
    double *sum, *coordinates;
    sum = (*clust).sum;
    coordinates = (*vec).coordinates;
    for (i = 0; i < (*clust).d; i++) {
        sum[i] += *coordinates;
        coordinates++;
    }
    (*clust).count++;
}

//input: Array of clusters, it's length and a vector.
//output: None. the function find the closest cluster to the vector
// and insert vector to cluster.
void findCluster(Cluster *clusterArray, Vector *vec, int arrayLen) {
    //declarations
    int i;
    Cluster *minCluster, *tempCluster;
    double minDistance, tempDistance;
    double distance(double *x, double *y, int length);
    void addVector(Cluster *clust, Vector *vec);

    //default minimum
    minCluster = clusterArray;
    tempCluster = minCluster;
    minDistance = distance((*vec).coordinates, (*minCluster).mean, (*minCluster).d);

    //loop through all clusters, except for 1st
    for (i = 1; i < arrayLen; i++) {
        tempCluster++;
        tempDistance = distance((*vec).coordinates, (*tempCluster).mean, (*tempCluster).d);
        if (tempDistance < minDistance) {
            minDistance = tempDistance;
            minCluster = tempCluster;
        }
    }
    addVector(minCluster, vec);
}

int main() {
    //test
    char ch;
    while (ch != '\n') {
        scanf("%c", &ch);
        printf("%c", ch);
    }
//    int *p = malloc(3*sizeof(int));
//    int i;
//    for(i=0;i<3;i++){
//        p[i] = i+1;
//    }
//
//    printf("%d", p[1]);
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
void kMeans(int k, int maxIter) {
    //TODO: Google when to declare functions (block below probably should be deleted) and where
    void printMeans(Cluster *clusterArray, int length);
    void initFromFile(int k, Cluster *clusterArray, Vector *vecArray);
    void findCluster(Cluster *clusterArray, Vector *vec, int arrayLen);
    void refreshClusters(Cluster *clusterArray, int clusterLength, int d);
    void freeMem(Cluster *clusterArray, Vector *vecArray);
    int reCalcMeans(Cluster *clusterArray, int length);

    Cluster *clusterArray = NULL;
    Vector *vecArray = NULL;
    initFromFile(k, clusterArray, vecArray);
    int lengthVectors; // number of vectors
    int lengthClusters; //number of clusters
    int i, iterCount;
    Vector *currentVec;
    currentVec = vecArray;
    int changed = 1;
    iterCount = 0;
    while (iterCount < maxIter && changed == 1) {
        //loop through vectors and insert to clusters
        for (i = 0; i < lengthVectors; i++) {
            findCluster(clusterArray, currentVec, lengthClusters);
        }
        //recalcmeans returns 1 if mean changed, 0 otherwise
        changed = reCalcMeans(clusterArray, lengthClusters);
        //init sum and count to zero's:
        refreshClusters(clusterArray, lengthClusters, (*clusterArray).d);
        iterCount++;
    }
    printMeans(clusterArray, lengthClusters);
    //TODO: continue from here
    freeMem(clusterArray, vecArray);
}
//functions

//input: 2 arrays(vectors coordinates)
//output: squered distance between them
//TODO: this works! :)
double distance(double *x, double *y, int length) {
    double sum = 0;
    int i;
    for (i = 0; i < length; i++) {
        sum += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sum;
}

//input: Array of clusters and it's length.number of coordinates.
//output: None. All clusters sum and count is initalized to zero
void refreshClusters(Cluster *clusterArray, int clusterLength, int d) {
    int i, j;
    Cluster *temp;
    double *indexer;
    temp = clusterArray;
    for (i = 0; i < clusterLength; i++) {
        (*temp).count = 0;
        indexer = ((*temp).sum);
        //loop through array of sum and init all to zero (0)
        for (j = 0; j < d; j++) {
            *indexer = 0;
            indexer++;
        }
        temp++;
    }
}

void initFromFile(int k, Cluster *clusterArray, Vector *vecArray) {
    //TODO:this should initialize vectors and clusters arrays, NOT FINISHED

    //find paramater d with first line
    int d = 0;
    int n = 0;
    char str;
    double number;
    const int LINE_MAX_LENGTH = 1000; // according to forum
    int i,j;
    for (i = 0; i < LINE_MAX_LENGTH; i++) { ;
        if (scanf("%c", &str) == EOF || str == '\n') {
            break; //reached end of line TODO:check if correct
        }
        if (str == ',') {
            d++;
        }
    }
    //found d
    rewind(stdin);//TODO:check if correct this way
    vecArray = malloc(sizeof(Vector) * 50); //init to size of 50 vectors
    double *arr = malloc(sizeof(double) * d); // array of coordinates
    int reached_end = 0; //"False"

    //Read vectors from file
    while (reached_end == 0) {
        for (i = 0; i < d; i++) {
            if (scanf("%lf%c", &number, &str) == EOF) {
                reached_end = 1;
                break;
            }
            arr[i] = number;
        }
        n++;
        Vector vec;
        vec.coordinates = arr;
        vecArray[n] = vec;//TODO:should be array of pointers
    }

    //init clusters
    for(i=0; i < k;i++){
        Cluster cl;
        //deep copy
        for(j=0;j<d;j++){
             cl.mean[j] = vecArray[i].coordinates[j];
        }
        clusterArray[i] = cl;
    }

}

void printMeans(Cluster *clusterArray, int length) {
    //TODO:didnt check
    int i, j, d;
    double *temp;
    d = clusterArray->d;
    //loop through clusters
    for (i = 0; i < length - 1; i++) {
        temp = (clusterArray[i]).mean;
        //loop through mean array
        for (j = 0; j < d; j++) {
            printf("%.4f,", *temp);
            temp++;
        }
        printf("%.4f", *temp);
    }
}

//deallocate memory
void freeMem(Cluster *clusterArray, Vector *vecArray) {
    free(clusterArray);
    free(vecArray);
}

