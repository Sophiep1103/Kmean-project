#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "spkmeans.h"

#define EPSILON 0.00001 /*EPSILON for convergence in the process of building A*/

/*define global variables*/
int k=0;
int iter = 0;
double epsilon = 0.0;
double* old_mean = NULL;
double **means = NULL;
double **points = NULL;
int *len_clusters = NULL;
double ***clusters =NULL;
int numPoints;
int pointLength;
double ** W = NULL;
double ** D = NULL;
double ** L = NULL;
double ** V = NULL;
double ** U = NULL;
double ** JACOBI = NULL;
double ** SPK = NULL;
double **dataPoints = NULL;
double ** A = NULL;
double ** P = NULL;
double ** resV = NULL;
double* eigenvalues= NULL;
double* eigenvalues_columns=NULL;
double ** A_prime = NULL;
int* i_j = NULL;
double c=0.0;
double s=0.0;

/*prototypes*/
double** wam_calc(double** dataPoints2, int numPoints2, int pointLenght2);
double** ddg_calc(double** dataPoints2, int numPoints2, int pointLength2);
double** gl_calc(double** dataPoints2, int numPoints2, int pointLength2);
double** jacobi_calc(double** dataPoints2, int numPoints2, int pointLength2);
double** spk_calc(double** dataPoints2, int arg_k, int numPoints2, int pointLength2);
double sqrtDistance(double* point, double* mean);
void calc_point_length(FILE *file);
void calc_num_points(FILE *file);
void print_matrix(double** matrix, int num_row, int num_col);
void make_data_points(char* file_name);
void buildP();
void calc_c_s();
void calc_pivot();
void calc_A_prime();
void copyMatBtoMatA(double** a, double** b);
double calc_off(double** matrix);
void update_eigenvalues();
void sort_eigenvalues_in_place();
void find_k();
int isConverged();
void updateV();
void updateU();
double** initial2Darray(int dim1, int dim2);
double distance(double* point, double* mean);
int find_closest_mean(double* point, double** means);
void calc_mean(double **cluster, int cluster_len, int index);
double** k_means(int k,int iter,double epsilon, int numPoints3, int pointLength3,double** means, double**points);
int update_centroids();
void assign_to_centroids();
void init_len_clusters();
void freeMemoryAndExit_kmeans();
void freeMemoryAndExit();
void freeMemory_kmeans();
void freeMemory();
void update_JACOBI();
void update_SPK();
void mult_matrix(double **mat1, double **mat2, double **res);
int sign(double val);
int MinusZero(double num);

/*calculates W, return w*/
double** wam_calc(double** dataPoints2, int numPoints2, int pointLength2){
    int i;
    int j;
    double distance;
    dataPoints = dataPoints2;
    numPoints = numPoints2;
    pointLength = pointLength2;
    W=initial2Darray(numPoints, numPoints);
    for (i=0; i<numPoints; i++){
        /*compute half of the matrix (symetric matrix)*/ 
        for (j=i; j<numPoints; j++){
            if (j!=i){
                distance= sqrtDistance(dataPoints[i], dataPoints[j]);
                W[i][j]= exp(-0.5*distance);
                W[j][i]=W[i][j]; /*W is symertic*/
            }
        }
    }
    return W;
}



/*calculates W and then D, returnes D*/
double** ddg_calc(double** dataPoints2, int numPoints2, int pointLength2){
    int i;
    int j;
    double sum=0;
    wam_calc(dataPoints2, numPoints2,pointLength2);
    D=initial2Darray(numPoints, numPoints);
    for (i=0; i<numPoints; i++){
        sum=0;
        for (j=0; j<numPoints; j++){
            sum+=W[i][j];
        }
        D[i][i]=sum;
    }
    return D;
}

/*calculates W, then D and then L, returnes L*/
double** gl_calc(double** dataPoints2, int numPoints2, int pointLength2){
    int i;
    int j;
    ddg_calc(dataPoints2, numPoints2,pointLength2);
    L=initial2Darray(numPoints, numPoints);
    for (i=0; i<numPoints; i++){
        for(j=0; j<numPoints; j++){
            L[i][j]=D[i][j]-W[i][j];
        }
    }
    return L;
}

/*calculates eigthenvalues and eightenvectors, return them in JACOBI matrix*/
double** jacobi_calc(double** dataPoints2, int numPoints2, int pointLength2){
    int num_iterations=0;
    A = dataPoints2;
    numPoints = numPoints2;
    pointLength = pointLength2;
    /*allocationg space for an array which will store the i and j indeces for calculating the jacobi matrix*/
    i_j=calloc(2, sizeof(int));
    if(i_j==NULL){
        freeMemoryAndExit();
    }
    /*allocationg space for P, A_prime , V , JACOBI and resV matrixes.*/
    P=initial2Darray(numPoints, numPoints);
    A_prime=initial2Darray(numPoints, numPoints);
    resV=initial2Darray(numPoints, numPoints); /*help matrix for calculate V*/
    V=initial2Darray(numPoints, numPoints);
    JACOBI= initial2Darray(numPoints+1, numPoints); /*stores the eigthenvalues and eightenvectors */

    while(num_iterations<100){
        calc_pivot();   
        calc_c_s(); /*calculate c and s values*/
        buildP(); /*build the P matrix*/
        if(num_iterations==0){ /*initialize V*/
            copyMatBtoMatA(V,P);
        }
        else{ /* V=V*P */
            mult_matrix(V, P, resV); /*resV=V*P*/
            copyMatBtoMatA(V,resV); /*V=resV*/
        }
        calc_A_prime();
        if(isConverged()==1){ /*if we reached the convergence condotion*/
            copyMatBtoMatA(A,A_prime);
            break;
        }
        else{
            copyMatBtoMatA(A,A_prime);
        }
        num_iterations+=1;
    }

    update_eigenvalues(); /*insert the diagnoal of P into eigenvalues array*/
    update_JACOBI(); /*insert the eigthenvalues and eightenvectors into JACOBI matrix */
    return JACOBI;
}

/*insert in to JACOBI the eigenvalues in the first row and then the eigenvectors as columns*/
void update_JACOBI(){
    int i;
    int j;
    for(i=0; i<numPoints; i++){
        JACOBI[0][i]=eigenvalues[i];
    }
    for(i=1; i<numPoints+1; i++){
        for(j=0; j<numPoints; j++){
        	JACOBI[i][j]=V[i-1][j];
        }

    }
}

/*calculates L from data points, calculates JACOBI from L ,returnes SPK matrix that contains the U matrix and also the k value*/
double** spk_calc(double** dataPoints2, int arg_k, int numPoints2, int pointLength2){
    gl_calc(dataPoints2, numPoints2, pointLength2); /*calculate L*/
    jacobi_calc(L, numPoints, numPoints); /*calculates JACOBI from L*/
    sort_eigenvalues_in_place(); /*sort eigenvalues and keep their new places in eigenvalues_columns*/
    if (arg_k==-1){ /*k wasn't given as an argument*/
        find_k();
    }
    else{
        k=arg_k;
    }
    U= initial2Darray(numPoints,k);
    SPK= initial2Darray(numPoints+1,k);
    updateU(); /*contains the first k eigenvectors*/
    update_SPK(); /*contains the U matrix and also the k value*/
    return SPK;
}

/*insert U and k in to SPK*/
void update_SPK(){
    int i;
    int j;

    for(i=0; i<numPoints; i++){
        for(j=0; j<k; j++){
            SPK[i][j]=U[i][j];
        }
    }
    SPK[i][0]=k;
}

/*initializes a 2 dimentional matrix*/
double** initial2Darray(int dim1, int dim2){
    int i;
    double** arr=NULL;
    arr= calloc(dim1,sizeof(double *));
    if (arr == NULL){
        freeMemoryAndExit();
    }
    for (i=0; i<dim1; i++){
        arr[i] = calloc(dim2,sizeof(double));
        if (arr[i]==NULL){
            freeMemoryAndExit();
        }
    }
    return arr;
}

/*insert first k eigtenvectors to U, treat the -0.000 case*/
void updateU(){
    int i;
    int j;
    int curr_col;
    double curr_value;
    int* indicator_minus_zero = calloc(k,sizeof(int)); /*indicator_minus_zero[i]==1 iff the i'th eigenvalue is -0.0*/
    char* minusZeroStr= calloc(50, sizeof(char));
    for (i = 0; i < k; i++){ /*first k eigenvectors*/
        curr_col = eigenvalues_columns[i];
        curr_value= eigenvalues[i];
        sprintf(minusZeroStr, "%.4f", curr_value); /*minusZeroStr=str(round(curr_value,4))*/
        if(strcmp(minusZeroStr, "-0.0000")==0){ /*if the i'th eigenvalue is approximate -0.0*/
            indicator_minus_zero[i]=1;
        }
        for (j = 0; j < numPoints; j++){ /*copy the k'th eigenvector into U*/
            if(indicator_minus_zero[i]==1){
                U[j][i] = (-1)*V[j][curr_col]; /*multiply the vector by -1*/
            }
            else{
                U[j][i] = V[j][curr_col];
            }
        }
    }
    free(indicator_minus_zero);
    free(minusZeroStr);
}


/*calcs A_prime from A,c,s following the rules*/
void calc_A_prime(){
    int r;
    int i= i_j[0];
    int j=i_j[1];
    copyMatBtoMatA(A_prime,A);
    for (r = 0; r < numPoints; r++){
        if(r!=i && r!=j){
            A_prime[r][i]=c*A[r][i]-s*A[r][j];
            A_prime[r][j]=c*A[r][j]+s*A[r][i];
            A_prime[i][r]=A_prime[r][i];
            A_prime[j][r]=A_prime[r][j];
        }
    }
    A_prime[i][i]=c*c*A[i][i]+s*s*A[j][j]-2*s*c*A[i][j];
    A_prime[j][j]=s*s*A[i][i]+c*c*A[j][j]+2*s*c*A[i][j];
    A_prime[i][j]=0;
    A_prime[j][i]=0;
}

/*find k by eigengap heuristic*/
void find_k(){
    double max_delta=0;
    int i;
    double delta;
    k=0;
    for(i=0; i<numPoints/2; i++){
        delta= fabs(eigenvalues[i]-eigenvalues[i+1]);
        if (delta>max_delta){
            max_delta= delta;
            k=i;
        }    
    }
    k=k+1; /*eigtenvalues are numbered 1,...,n*/
}

/*sort eigenvalues and eigenvalues_columns in place*/
void sort_eigenvalues_in_place(){
    int i;
    int j;
    double tmp1;
    double tmp2;

    eigenvalues_columns=calloc(numPoints, sizeof(double));
    if (eigenvalues_columns == NULL){
        freeMemoryAndExit();
    }

    for(i=0; i<numPoints; i++){ /*initialize eigenvalues_columns to 0,....,numPoints-1*/
        eigenvalues_columns[i]=i;
    }

    for (i = 0; i < numPoints-1; i++){ /*bubble sort of eigenvalues and their columns*/
        for (j = 0; j < numPoints-1-i; j++){
            if (eigenvalues[j] > eigenvalues[j+1]){
                tmp1 = eigenvalues[j+1];
                tmp2 = eigenvalues_columns[j+1];
                eigenvalues[j+1] = eigenvalues[j];
                eigenvalues_columns[j+1] = eigenvalues_columns[j];
                eigenvalues[j] = tmp1;
                eigenvalues_columns[j] = tmp2;
            }
        }
    }
}

/*inserts the diagnol of A to eigenvalues array*/
void update_eigenvalues(){
    int i;
    eigenvalues=calloc(numPoints, sizeof(double));
    if (eigenvalues == NULL){
        freeMemoryAndExit();
    }
    for (i=0; i<numPoints; i++){
        eigenvalues[i]=A[i][i];
    }
}

/*update P following the rules*/
void buildP(){
    int k;
    int z;
    int i;
    int j;
    i = i_j[0]; /*the indice of the pivot*/
    j = i_j[1];
    for (k = 0; k < numPoints; k++) {
        for (z = 0; z < numPoints; z++) {
            if (k == z) {
                if (k == i || k == j) {
                    P[k][k] = c;
                } 
                else {
                    P[k][k] = 1;
                }
            }
            else {
                P[k][z] = 0.0;
            }
        }
    }
    if (s != 0.0){
        P[i][j] = s;
        P[j][i] = -1 * s;
    }
}

/*multiply mat1 and mat2 and insert the result to res*/
void mult_matrix(double **mat1, double **mat2, double **res){
    int i;
    int j;
    int z;
    double sum = 0.0;
    for (i = 0; i < numPoints; i++){
        for (j = 0; j < numPoints; j++){
            for (z = 0; z < numPoints; z++){
                sum += mat1[i][z] * mat2[z][j];
            }
            res[i][j] = sum;
            sum = 0.0;
        }
    }
}

/*multiply V by P*/
void updateV(){
    int i;
    int j;
    int z;
    double sum;
    copyMatBtoMatA(resV, V);
    sum = 0.0;
    for (i = 0; i < numPoints; i++){
        for (j = 0; j < numPoints; j++){
            for (z = 0; z < numPoints; z++){
                sum += resV[i][z] * P[z][j];
            }
            V[i][j] = sum;
            sum = 0.0;
        }
    }
}

/*calculate s and c as defined in the assignment*/
void calc_c_s(){
    double theta;
    double t;
    int i= i_j[0];
    int j= i_j[1];
    if (A[i][j] == 0) {
        s = 0;
        c = 1;
    }
    else{
        theta = (double)((A[j][j]-A[i][i]) / (2*A[i][j]));
        t = (double)(sign(theta) / (fabs(theta) + sqrt(pow(theta, 2) + 1)));
        c= (double)(1.0 / sqrt(pow(t, 2) + 1));
        s=t*c;
    }
}

/*if val>=0: return 1. else: return -1*/
int sign(double val) {
    if (val >= 0.0) {
        return 1;
    }
    else {
        return -1;
    }
}

/*create a copy of matrix b in matrix a*/
void copyMatBtoMatA(double** a, double** b){
    int i;
    int j;
    for (i = 0; i < numPoints; i++){
        for (j = 0; j < numPoints; j++){
            a[i][j] = b[i][j];
        }
    }
}

/*calcs the pivot's indexes following the rules*/
void calc_pivot(){
    int i;
    int j;
    double max_num=-1.0;
    for (i=0; i<numPoints; i++){
        for(j=i+1; j<numPoints; j++){ /* going throw the upper triangle of the matrix*/
            if(fabs(A[i][j])>max_num){ /*> and not >= to make sure we get the upper left element in case of equality*/
                max_num=fabs(A[i][j]);
                i_j[0]=i;
                i_j[1]=j;
            }
        }
    }
}

/*returnes 1 if is converged, 0 otherwise*/
int isConverged(){
    if(calc_off(A)-calc_off(A_prime)<=EPSILON){
        return 1;
    }
    return 0;
}

/*calcs off(matrix)^2*/
double calc_off(double** matrix){
    int i;
    int j;
    double sum=0;
    for(i=0; i<numPoints; i++){
        for(j=0; j<numPoints; j++){
            if (i!=j){
                sum+=matrix[i][j]*matrix[i][j];
            }
        }
    }
    return sum;
}

/*calculate the Euclidean sqrt Distance d(point, mean)*/ 
double sqrtDistance(double* point, double* mean){
    double dist = 0 ;
    int i;
    for (i=0; i<pointLength; i++){
        double res = pow((double)(point[i] - mean[i]), 2.0);
        dist = dist+ res;
        
    }
    return dist;
}


/*free memory from k_means algorithem and exit with error*/
void freeMemoryAndExit_kmeans() {
    printf("An Error Has Occurred");
    freeMemory_kmeans();
    exit(1);
}

/*free memory from spk algoritem and exit with error*/
void freeMemoryAndExit() {
    printf("An Error Has Occurred");
    freeMemory();
    exit(1);
}

/*calculates the length of a point in the data*/
void calc_point_length(FILE *file){
    int d, c;
    d = 1;
    while ((c = fgetc(file)) != EOF) { 
        if (c == ',') {
            d++;
        }
        if (c == '\n') {
            break;
        }
    }
    pointLength= d;
}


/*calculate the number of points in the data*/
void calc_num_points(FILE *file){
    int n, c;
    n = 1;
    while ((c = fgetc(file)) != EOF) {
        if (c == '\n') {
            n++;
        }
    }
    numPoints = n-1;
}

/*prints matrix*/
void print_matrix(double** matrix, int num_row, int num_col){
    int i,j;
    for (i=0; i<num_row; i++){
        for (j=0; j<num_col; j++){
            if (j!=num_col-1)
                printf("%.4f,", matrix[i][j]); 
            else
                printf("%.4f", matrix[i][j]); 

        }
        printf("\n");
    }
}

/*create a matrix of the points in the data*/
void make_data_points(char* file_name){
    int i;
    int j;
    int num=0;
    FILE *file = fopen(file_name, "r");
    calc_num_points(file);
    rewind(file);
    calc_point_length(file);
    rewind(file);

    dataPoints=initial2Darray(numPoints, pointLength);

    for (i=0; i<numPoints; i++){
        for (j=0; j<pointLength; j++){
            num=fscanf(file, "%lf,", &(dataPoints[i][j]));
            if (num<0){
                freeMemoryAndExit();
            }
        }
    }
    fclose(file);
}

/*k_means algorithem, return the means of the final clusters*/
double** k_means(int k1,int iter1,double epsilon1, int numPoints1, int pointLength1,double** means1, double**points1){
    int i;
    int is_not_converge; /*indicates if the distance is smaller than epsilon.*/
    int count_loop; /*counts the num of iterations*/
    k = k1;
    iter=iter1; /*will be 300*/
    epsilon= epsilon1; /*will be 0*/
    pointLength= pointLength1;
    means= means1; /*first means calculated by k_means++*/
    points =points1;
    numPoints= numPoints1;
    is_not_converge = 1;
    count_loop =0;

    clusters = calloc(k, sizeof(double **));
    if(clusters==NULL){
        freeMemoryAndExit();
    }
    len_clusters = calloc(k, sizeof(int));
    if(len_clusters==NULL){
        freeMemoryAndExit();
    }
    for (i = 0; i < k; i++) {
        clusters[i] = malloc(sizeof(double*));
        if (clusters[i] == NULL) {
            freeMemoryAndExit();
        }
    }
    old_mean = calloc(pointLength, sizeof(double));
    if (old_mean == NULL) {
        freeMemoryAndExit();
    }

    while((is_not_converge==1) && (count_loop<iter)){
        count_loop++;
        assign_to_centroids(); /*assign points to the closest centroids*/
        is_not_converge = update_centroids(); /*update the new means of the centroids and checks convergence*/
        init_len_clusters(); /*initialize the array that keeps the length of the clusters*/
    }
    return means; /*returns final means*/
}

/*calculate the Euclidean Distance d(point, mean)*/ 
double distance(double* point, double* mean){
    double dist = 0 ;
    int i;
    for (i=0; i<pointLength; i++){
        double res = pow((double)(point[i] - mean[i]), 2.0);
        dist = dist+ res;
        
    }
    return sqrt(dist);
}

/*finds the closest mean to a specific point from all the means that are currently in **means.
returns the (int) index (in **mean array) of this closest mean. */ 
int find_closest_mean( double* point, double** means){
    double min_dist = 1.79769e+308;
    int index_min_point = 0 ;
    int i;
    double dist;
    for (i = 0; i<k; i++){
        /*calculates the distance of the point from every mean in order to find the closest mean*/
        dist = distance(point, means[i]);
        if (dist< min_dist){
            min_dist = dist;
            index_min_point=i;
        }
    }
    return index_min_point;
}

/* calculates the mean of the points in cluster, index is the index of the cluster in clusters array*/
void calc_mean(double **cluster, int cluster_len, int index){
    int i;
    int j;
    int len;

    for(i=0; i<pointLength; i++){
        means[index][i]=0;
    }

    len= cluster_len;
    for (i=0; i<len; i++){
        for(j=0; j<pointLength; j++){
            means[index][j] += cluster[i][j];
        }
    }
    for (i=0; i<pointLength; i++){
         means[index][i] = (double)(means[index][i]/len);
    }
}

/*assign all the points to the correct cluster. Here, we update the clusters*/ 
void assign_to_centroids(){
    int i;
    int closest_mean_index;
    for (i=0; i<numPoints; i++){ /*iterates over all the points.*/
        closest_mean_index = find_closest_mean(points[i], means); /*find the closest mean to a point*/
        /*adding (one) more space to the correct cluster in order to add this point to the cluster*/
        clusters[closest_mean_index] = realloc(clusters[closest_mean_index],len_clusters[closest_mean_index]*sizeof(double*)+sizeof(double*));
        if(clusters[closest_mean_index]==NULL){
            freeMemoryAndExit();
        }
        len_clusters[closest_mean_index]++; /*updating the array of the clusters lengths.*/
        clusters[closest_mean_index][len_clusters[closest_mean_index]-1]=points[i]; /*adding the point to a cluster of the closest mean.*/
    }
}

/*initialize len_clusters to be zeroes*/
void init_len_clusters(){
    int i;
    for (i=0; i<k; i++){
        len_clusters[i]=0;
    }
}


/*updated the centroids (means of the clusters) after assignning every point to a cluster. 
check convergence and returns 0 if we reached convergence, 1 otherwise*/
int update_centroids(){
    int i;
    int j;
    double cord;
    int is_not_converge = 0;
    for(i=0; i<k; i++){
        for(j=0; j<pointLength; j++){
            cord=means[i][j];
            old_mean[j]=cord; /*old mean saves a specific mean temporary*/
        }
        calc_mean(clusters[i], len_clusters[i], i); /*calcs the new mean of the cluster*/
        if(distance(old_mean, means[i])> epsilon) /*check convergence*/
            is_not_converge = 1; /*if one mean at least didnt converged- return 1*/
    }
    return is_not_converge; /*if all the means did converged- return 0*/
}

/*free memory of all arrays in k_means algorithem*/
void freeMemory_kmeans() {
    int i;
    if (points != NULL) {
        for (i = 0; i < numPoints; i++) {
            if(points[i]!=NULL){
                free(points[i]);
            }
        }
        free(points);
    }
    if (clusters != NULL){
        for (i = 0; i < k; i++) {
            if(clusters[i]!=NULL){
                free(clusters[i]);
            }
        }
        free(clusters);
    }
    if (means != NULL){
        for (i = 0; i < k; i++) {
            free(means[i]);
        }
        free(means);
    }
    if (len_clusters != NULL) {
        free(len_clusters);
    }
    if (old_mean != NULL){
        free(old_mean);
    }
}


/*free memory of all arrays in k_means algorithem*/
void freeMemory() {
    int i;
    if (dataPoints != NULL) {
        for (i = 0; i < numPoints; i++) {
            if(dataPoints[i]!=NULL){
                free(dataPoints[i]);
            }
        }
        free(dataPoints);
    }
    if (W!= NULL){
        for (i = 0; i < numPoints; i++) {
            free(W[i]);
        }
        free(W);
    }
    if (L!= NULL){
        for (i = 0; i < numPoints; i++) {
            free(L[i]);
        }
        free(L);
    }
    if (D!= NULL){
        for (i = 0; i < numPoints; i++) {
            free(D[i]);
        }
        free(D);
    }
    if (V!= NULL){
        for (i = 0; i < numPoints; i++) {
            free(V[i]);
        }
        free(V);
    }
    if (U!= NULL){
        for (i = 0; i < numPoints; i++) {
            free(U[i]);
        }
        free(U);
    }
    if (P!= NULL){
        for (i = 0; i < numPoints; i++) {
            free(P[i]);
        }
        free(P);
    }
    if (resV!= NULL){
        for (i = 0; i < numPoints; i++) {
            free(resV[i]);
        }
        free(resV);
    }
    if (A_prime!= NULL){
        for (i = 0; i < numPoints; i++) {
            free(A_prime[i]);
        }
        free(A_prime);
    }
    if (JACOBI!= NULL){
        for (i = 0; i < numPoints+1; i++) {
            free(JACOBI[i]);
        }
        free(JACOBI);
    }
    if (SPK!= NULL){
        for (i = 0; i < numPoints+1; i++) {
            free(SPK[i]);
        }
        free(SPK);
    }
    if(eigenvalues != NULL){
        free(eigenvalues);
    }
    if(eigenvalues_columns != NULL){
        free(eigenvalues_columns);
    }
    if(i_j != NULL){
        free(i_j);
    }
        if (points != NULL) {
        for (i = 0; i < numPoints; i++) {
            if(points[i]!=NULL){
                free(points[i]);
            }
        }
        free(points);
    }
    if (clusters != NULL){
        for (i = 0; i < k; i++) {
            if(clusters[i]!=NULL){
                free(clusters[i]);
            }
        }
        free(clusters);
    }
    if (means != NULL){
        for (i = 0; i < k; i++) {
            free(means[i]);
        }
        free(means);
    }
    if (len_clusters != NULL) {
        free(len_clusters);
    }
    if (old_mean != NULL){
        free(old_mean);
    }
}


/*get goal, filename from Command line arguments- returns the required matrix*/
int main(int argc, char** argv){
    int i;
    char* goal;
    char* file_name;
    if(argc!=3){
        printf("An Error Has Occurred");
        return 1;
    }
    goal=malloc(strlen(argv[1])+1);
    file_name=malloc(strlen(argv[2])+1);
    strcpy(goal,argv[1]);
    strcpy(file_name,argv[2]);
    make_data_points(file_name); /*build the array of points*/
    if (strcmp(goal, "wam")==0){ /*goal==wam*/
        wam_calc(dataPoints , numPoints , pointLength);
        print_matrix(W, numPoints, numPoints);
    }
    else if (strcmp(goal, "ddg")==0){ /*goal==ddg*/
        ddg_calc(dataPoints , numPoints , pointLength);
        print_matrix(D, numPoints, numPoints);
    }
    else if (strcmp(goal, "gl")==0){ /*goal==gl*/
        gl_calc(dataPoints , numPoints , pointLength);
        print_matrix(L, numPoints, numPoints);
    }
    else if (strcmp(goal, "jacobi")==0){ /*goal==jacobi*/
        jacobi_calc(dataPoints , numPoints , pointLength);
        for (i=0; i<numPoints-1; i++){ /*print eigenvalues*/
            printf("%0.4f,", eigenvalues[i]);
        }
        printf("%0.4f\n", eigenvalues[i]);
        print_matrix(V, numPoints, numPoints); /*print eigenvectors*/
    }
    free(goal);
    free(file_name);
    freeMemory();
    return 0;
}

