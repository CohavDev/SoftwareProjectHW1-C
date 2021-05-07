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

/** global variables **/
double **vecArray;
Cluster **clusterArray;
int k,d,n;
int maxIter;

/**function declarations**/
int reCalcMeans();
int calcMean(Cluster *clust);
void addVector(Cluster *clust, double *vec);
void findCluster(double *vec);
void kMeans(int maxIter);
double distance(double *x, double *y, int length);
void refreshClusters();
void initFromFile();
void printMeans();
void freeMem();

/** main **/
int main(int argc, char *argv[]) {
//    printf("%s", argv[1]);
    k = atoi(argv[1]);
    if(argc == 3){
        maxIter = atoi(argv[2]);
    }
    else{
        maxIter = 200;//TODO:check default value instructions
    }
//    printf("%d , %d",k,maxIter);
    kMeans(maxIter);
    return 0;
}

/** functions **/

//input: Array of pointers of clusters and its length
//output: None, calcs means of all clusters
//TODO: didnt check correctness
int reCalcMeans() {
    int changed = 0; // "False"
    int i;
    for (i = 0; i < k; i++) {
        if (calcMean(clusterArray[i]) == 1) {
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
    for (i = 0; i < d; i++) {
        calcVal = (sum[i] / (*clust).count);
        if (mean[i] != calcVal) {
            changed = 1; //"True"
        }
        mean[i] = calcVal;
    }
    return changed;
}

//adds specific vector to specific cluster
void addVector(Cluster *clust, double *vec) {
    int i;
    double *sum;
    sum = (*clust).sum;
    for (i = 0; i < (*clust).d; i++) {
        sum[i] += vec[i];
    }
    (*clust).count++;
}

//input: Array of clusters, it's length and a vector.
//output: None. the function find the closest cluster to the vector
// and insert vector to cluster.
void findCluster(double *vec) {
    //declarations
    int i;
    Cluster *minCluster, *tempCluster;
    double minDistance, tempDistance;

    //default minimum
    minCluster = *clusterArray;
    tempCluster = minCluster;
    minDistance = distance(vec, (*minCluster).mean, (*minCluster).d);

    //loop through all clusters, except for 1st
    for (i = 1; i < k; i++) {
        tempCluster++;
        tempDistance = distance(vec, (*tempCluster).mean, (*tempCluster).d);
        if (tempDistance < minDistance) {
            minDistance = tempDistance;
            minCluster = tempCluster;
        }
    }
    addVector(minCluster, vec);
}

//kMeans function
void kMeans(int maxIter) {
    Cluster **clusterArray = NULL;
    double **vecArray = NULL;
    initFromFile(k, clusterArray, vecArray);
    int i, iterCount;
    double *currentVec;
    currentVec = *vecArray;
    int changed = 1;
    iterCount = 0;
    while (iterCount < maxIter && changed == 1) {
        //loop through vectors and insert to clusters
        for (i = 0; i < n; i++) {
            findCluster(currentVec);
        }
        //recalcmeans returns 1 if mean changed, 0 otherwise
        changed = reCalcMeans(clusterArray, k);
        //init sum and count to zero's:
        refreshClusters(clusterArray, k, (*clusterArray)->d);
        iterCount++;
    }
    printMeans(clusterArray, k);
    //TODO: continue from here
    freeMem(clusterArray, vecArray);
}
//functions

//input: 2 arrays(vectors coordinates)
//output: squered distance between them
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
void refreshClusters() {
    int i, j;
    Cluster *temp;
    double *indexer;
    temp = *clusterArray;
    for (i = 0; i < k; i++) {
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

void initFromFile() {
    //TODO:this should initialize vectors and clusters arrays, NOT FINISHED

    //find paramater d with first line
    d = 0;
    n = 0;
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
    int size = 50;
    vecArray = malloc(sizeof(double*) * size); //init to size of size = 50 pointers
    double *arr = malloc(sizeof(double) * d); // array of coordinates
    int reached_end = 0; //"False"

    //Read vectors from file
    while (reached_end == 0) {
        if(n == size){
            vecArray = realloc(vecArray, size*2*sizeof(double*));
            size += size;
        }
        for (i = 0; i < d; i++) {
            if (scanf("%lf%c", &number, &str) == EOF) {
                reached_end = 1;
                break;
            }
            arr[i] = number;
        }
        n++;
        vecArray[n] = arr;
    }
    vecArray = realloc(vecArray, sizeof(double*) * n);

    //init clusters
    clusterArray = (Cluster**)malloc(k * sizeof(Cluster*));
    for(i=0; i < k;i++){
        Cluster cl;
        //deep copy
        for(j=0;j<d;j++){
             cl.mean[j] = vecArray[i][j];
        }
        clusterArray[i] = &cl;
    }

}

void printMeans() {
    //TODO:didnt check
    int i, j;
    double *temp;
    //loop through clusters
    for (i = 0; i < k; i++) {
        temp = (clusterArray[i])->mean;
        //loop through mean array
        for (j = 0; j < d-1; j++) {
            printf("%.4f,", temp[j]);
        }
        printf("%.4f", temp[d]);
    }
}

//deallocate memory
//TODO:need fix
void freeMem() {
    free(clusterArray);
    free(vecArray);
}

