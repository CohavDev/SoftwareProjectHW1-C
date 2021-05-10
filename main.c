#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


struct Cluster {
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
void kMeans();
double distance(double *x, double *y);
void refreshClusters();
void initFromFile();
void printMeans();
void freeMem();

/** main **/
int main(int argc, char *argv[]) {
    //TODO:validate cmd arguments
    assert(argc == 2 || argc == 3);
    k = atoi(argv[1]);
    if(argc == 3){
        maxIter = atoi(argv[2]);
    }
    else{
        maxIter = 200; //default value
    }
    kMeans();
    return 0;
}

/** functions **/

//input: Array of pointers of clusters and its length
//output: None, calcs means of all clusters
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
    sum = clust->sum;
    mean = clust->mean;
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
    for (i = 0; i < d; i++) {
        sum[i] = sum[i] + vec[i];
    }
    (*clust).count++;
}

//input: Array of clusters, it's length and a vector.
//output: None. the function find the closest cluster to the vector
// and insert vector to cluster.
void findCluster(double *vec) {
    //declarations
    int i;
    Cluster *minCluster;
    double minDistance, tempDistance;

    //default minimum
    minCluster = clusterArray[0];
    minDistance = distance(vec, minCluster->mean);

    //loop through all clusters, except for 1st
    for (i = 1; i < k; i++) {
        tempDistance = distance(vec, (clusterArray[i])->mean);
        if (tempDistance < minDistance) {
            minDistance = tempDistance;
            minCluster = clusterArray[i];
        }
    }
    addVector(minCluster, vec);
}

//kMeans function
void kMeans() {
    initFromFile(k, clusterArray, vecArray);
    int i, iterCount;
    int changed = 1;
    iterCount = 0;
    while (iterCount < maxIter && changed == 1) {
        //loop through vectors and insert to clusters
        for (i = 0; i < n; i++) {
            findCluster(vecArray[i]);
        }
        changed = reCalcMeans(); //re-calculate means. returns 1 if mean changed, 0 otherwise
        refreshClusters(); //init sum and count to zero's:
        iterCount++;
    }
    printMeans();
    //TODO: continue from here
    freeMem(clusterArray, vecArray);
}
//functions

//input: 2 arrays(vectors coordinates)
//output: squered distance between them
double distance(double *x, double *y) {
    double sum = 0;
    int i;
    for (i = 0; i < d; i++) {
        sum += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sum;
}

//input: Array of clusters and it's length.number of coordinates.
//output: None. All clusters sum and count is initalized to zero
void refreshClusters() {
    int i, j;
    for (i = 0; i < k; i++) {
        clusterArray[i]->count = 0;
        //loop through array of sum and init all to zero (0)
        for (j = 0; j < d; j++) {
            (clusterArray[i]->sum)[j] = 0;
        }
    }
}

void initFromFile() {
    //find paramater d with first line
    d = 0;
    n = 0;
    char str;
    double number;
    const int LINE_MAX_LENGTH = 1000; // according to forum
    int i,j;
    for (i = 0; i < LINE_MAX_LENGTH; i++) { ;
        if (scanf("%lf%c", &number, &str) == EOF || str == '\n') {
            d++;
            break;
        }
        if (str == ',') {
            d++;
        }
    }

    rewind(stdin);
    int size = 50;
    vecArray = malloc(sizeof(double*) * size); //init to size of size
    assert(vecArray !=NULL);
    double *arr;
    int reached_end = 0; //"False"

    //Read vectors from file
    while (reached_end == 0) {
        arr  = malloc(sizeof(double) * d); // array of coordinates
        if(n == size){
            vecArray = realloc(vecArray, (2*size*sizeof(double*)));
            assert(vecArray !=NULL);
            size += size;
        }
        for (i = 0; i < d; i++) {
            if (scanf("%lf%c", &number, &str) == EOF) {
                reached_end = 1;
                break;
            }
            arr[i] = number;

        }
        if(reached_end != 1){
            n++;
            vecArray[n-1] = arr;
        }

    }
    vecArray = realloc(vecArray, sizeof(double*) * n);
    assert(vecArray != NULL);

    //init clusters
    clusterArray = malloc(k * sizeof(Cluster*));
    assert(clusterArray != NULL);
    for(i=0; i < k;i++){
        Cluster *cl = malloc(sizeof(Cluster));
        cl->mean = malloc(d*sizeof(double));
        cl->sum = calloc(d,sizeof(double ));
        assert(cl->sum !=NULL && cl->mean != NULL);
        //deep copy
        for(j=0;j<d;j++){
            (cl->mean)[j] = vecArray[i][j];
        }
        cl->count = 0;
        clusterArray[i] = cl;
    }
}

void printMeans() {
    int i, j;
    double *temp;
    //loop through clusters
    for (i = 0; i < k; i++) {
        temp = (clusterArray[i])->mean;
        //loop through mean array
        for (j = 0; j < d-1; j++) {
            printf("%0.4lf,", temp[j]);
        }
        printf("%0.4lf\n", temp[d-1]);
    }
}

//deallocate memory
//TODO:need fix
void freeMem() {
    free(clusterArray);
    free(vecArray);
}

