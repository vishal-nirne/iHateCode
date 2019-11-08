/* 
   for compilation of code gcc -fopenmp mat.c
   for execution of the code ./a.out 

   the minimum time taken is   time_taken 94.762835	
   the maximum time taken is  time_taken 121.966892
*/


#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
                                                      
//size of the matrix
#define N 512

//it all begins from here

//memory allocation to the matrix
double **alloc_matrix(int size){
    double** matrix;
    matrix = (double**) malloc(size * sizeof(double *));
    if(matrix == NULL){
    	printf("\n Error in mem alloc\n");
    	return NULL;
    }

    for (int i = 0; i < size; i++){
    	matrix[i] = (double*) malloc(size*sizeof(double));
        //matrix[i] = matrix[0] + i*size;
        //printf(" %lf\n",**matrix);
    }
    return matrix;
}

//set to NULL and free the pointers
double **free_matrix(double **mat, int n) {

    int i;

    if (mat == NULL) {
        return (NULL);
    }

    for (i = 0; i < n; i++) { 
       //printf("\nFree row %d",i);
        free(mat[i]); 
            
    } 

    free(mat);         
    mat = NULL;

    return (NULL);   
}

void sum(double **a, double **b, double **result, int sz) {
    int i, j;

    	for (i = 0; i < sz; i++) {
        	for (j = 0; j < sz; j++) {
            	result[i][j] = a[i][j] + b[i][j];
        	}
    	}
}

void subtract(double **a, double **b, double **result, int sz) {
    int i, j;

    for (i = 0; i < sz; i++) {
        for (j = 0; j < sz; j++) {
            result[i][j] = a[i][j] - b[i][j];
        }
    }	
}


double** strassen(double **a, double **b, double **c, int tam) {

    // trivial case: when the matrix is 1 X 1:
    if (tam == 1) {
        c[0][0] = a[0][0] * b[0][0];
        return c;
    }

    // other cases are treated here:
    else {
        int newTam = tam/2;
        double **a11, **a12, **a21, **a22;
        double **b11, **b12, **b21, **b22;
        double **c11, **c12, **c21, **c22;
        double **p1, **p2, **p3, **p4, **p5, **p6, **p7;

        // memory allocation:
        a11 = alloc_matrix(newTam);
        a12 = alloc_matrix(newTam);
        a21 = alloc_matrix(newTam);
        a22 = alloc_matrix(newTam);

        b11 = alloc_matrix(newTam);
        b12 = alloc_matrix(newTam);
        b21 = alloc_matrix(newTam);
        b22 = alloc_matrix(newTam);

        c11 = alloc_matrix(newTam);
        c12 = alloc_matrix(newTam);
        c21 = alloc_matrix(newTam);
        c22 = alloc_matrix(newTam);

        p1 = alloc_matrix(newTam);
        p2 = alloc_matrix(newTam);
        p3 = alloc_matrix(newTam);
        p4 = alloc_matrix(newTam);
        p5 = alloc_matrix(newTam);
        p6 = alloc_matrix(newTam);
        p7 = alloc_matrix(newTam);

        double **aResult = alloc_matrix(newTam);
        double **bResult = alloc_matrix(newTam);

        int i, j;

        //dividing the matrices in 4 sub-matrices:
        for (i = 0; i < newTam; i++) {
            for (j = 0; j < newTam; j++) {
                a11[i][j] = a[i][j];
                a12[i][j] = a[i][j + newTam];
                a21[i][j] = a[i + newTam][j];
                a22[i][j] = a[i + newTam][j + newTam];

                b11[i][j] = b[i][j];
                b12[i][j] = b[i][j + newTam];
                b21[i][j] = b[i + newTam][j];
                b22[i][j] = b[i + newTam][j + newTam];
            }
        }

        // Calculating p1 to p7:

        sum(a11, a22, aResult, newTam); // a11 + a22
        sum(b11, b22, bResult, newTam); // b11 + b22
        strassen(aResult, bResult, p1, newTam); // p1 = (a11+a22) * (b11+b22)

        sum(a21, a22, aResult, newTam); // a21 + a22
        strassen(aResult, b11, p2, newTam); // p2 = (a21+a22) * (b11)

        subtract(b12, b22, bResult, newTam); // b12 - b22
        strassen(a11, bResult, p3, newTam); // p3 = (a11) * (b12 - b22)

        subtract(b21, b11, bResult, newTam); // b21 - b11
        strassen(a22, bResult, p4, newTam); // p4 = (a22) * (b21 - b11)

        sum(a11, a12, aResult, newTam); // a11 + a12
        strassen(aResult, b22, p5, newTam); // p5 = (a11+a12) * (b22)	

        subtract(a21, a11, aResult, newTam); // a21 - a11
        sum(b11, b12, bResult, newTam); // b11 + b12
        strassen(aResult, bResult, p6, newTam); // p6 = (a21-a11) * (b11+b12)

        subtract(a12, a22, aResult, newTam); // a12 - a22
        sum(b21, b22, bResult, newTam); // b21 + b22
        strassen(aResult, bResult, p7, newTam); // p7 = (a12-a22) * (b21+b22)

        // calculating c21, c21, c11 e c22:

        sum(p3, p5, c12, newTam); // c12 = p3 + p5
        sum(p2, p4, c21, newTam); // c21 = p2 + p4

        sum(p1, p4, aResult, newTam); // p1 + p4
        sum(aResult, p7, bResult, newTam); // p1 + p4 + p7
        subtract(bResult, p5, c11, newTam); // c11 = p1 + p4 - p5 + p7

        sum(p1, p3, aResult, newTam); // p1 + p3
        sum(aResult, p6, bResult, newTam); // p1 + p3 + p6
        subtract(bResult, p2, c22, newTam); // c22 = p1 + p3 - p2 + p6

        // Grouping the results obtained in a single matrix:
        for (i = 0; i < newTam ; i++) {
            for (j = 0 ; j < newTam ; j++) {
                c[i][j] = c11[i][j];
                c[i][j + newTam] = c12[i][j];
                c[i + newTam][j] = c21[i][j];
                c[i + newTam][j + newTam] = c22[i][j];
            }
        }

        // deallocating memory (free):
        a11 = free_matrix(a11, newTam);
        a12 = free_matrix(a12, newTam);
        a21 = free_matrix(a21, newTam);
        a22 = free_matrix(a22, newTam);

        b11 = free_matrix(b11, newTam);
        b12 = free_matrix(b12, newTam);
        b21 = free_matrix(b21, newTam);
        b22 = free_matrix(b22, newTam);

        c11 = free_matrix(c11, newTam);
        c12 = free_matrix(c12, newTam);
        c21 = free_matrix(c21, newTam);
        c22 = free_matrix(c22, newTam);

        p1 = free_matrix(p1, newTam);
        p2 = free_matrix(p2, newTam);
        p3 = free_matrix(p3, newTam);
        p4 = free_matrix(p4, newTam);
        p5 = free_matrix(p5, newTam);
        p6 = free_matrix(p6, newTam);
        p7 = free_matrix(p7, newTam);
        aResult = free_matrix(aResult, newTam);
        bResult = free_matrix(bResult, newTam);
    } // end of else
    return c;
} // end of Strassen function


void display_matrix(double **mat){

	printf("printing only first 4 results \n");
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){

			
			printf("%lf \t",mat[i][j]);
		}
		printf("\n");
	}
}

int main(int argc, char* argv[]){
	
	double **a,**b,**c;
	int tid,nthreads;
	bool correct=false;
	clock_t t;	
	srand(time(NULL));
	
	printf("\n MATRIX MULTIPLICATION USING STRASSEN'S METHOD with N=%d\n",N);

		tid = omp_get_thread_num();
	
	if(tid == 0){
    	//printf("Initializing matrices...\n");
    a = alloc_matrix(N);
	b = alloc_matrix(N);
	c = alloc_matrix(N);

	}

	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			a[i][j]= 1;
			//fscanf(file,"lf",&a[i][j]);
		}
	}
	
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			b[i][j]= 2;
			//fscanf(file,"lf",&b[i][j]);
		}
	}
	
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			c[i][j] = 0.0;
		}
	}
		
	printf("\n The matrix A is\n");
	display_matrix(a);
	printf("\n The matrix B is\n");
	display_matrix(b);

	//multiply the matrices
	t = clock();
	//matrix_multiply(a,b,c,N);
	strassen(a,b,c,N);
	t = clock()-t;

	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
	printf("**************strassen matrix multiplication****************");
	printf("\n time_taken %f seconds",time_taken);
	printf("\n");
	printf("\n The matrix C is\n");
	display_matrix(c);
	
	free_matrix(a,N);
	free_matrix(b,N);
	free_matrix(c,N);
	
	return 0;

}