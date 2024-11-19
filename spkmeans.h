
#ifndef SPKMEANS_H
#define SPKMEANS_H
/*prototypes of functions that are used in spkmeansmodule.c and are written in spkmeans.c*/
double** wam_calc(double** dataPoints2, int numPoints2, int pointLenght2);
double** ddg_calc(double** dataPoints2, int numPoints2, int pointLength2);
double** gl_calc(double** dataPoints2, int numPoints2, int pointLength2);
double** jacobi_calc(double** dataPoints2, int numPoints2, int pointLength2);
double** spk_calc(double** dataPoints2, int arg_k, int numPoints2, int pointLength2);
double** initial2Darray(int dim1, int dim2);
double** k_means(int k,int iter,double epsilon, int numPoints3, int pointLength3,double** means, double**points);
void freeMemoryAndExit_kmeans();
void freeMemory();
void freeMemoryAndExit();
void freeMemory_kmeans();

#endif 
