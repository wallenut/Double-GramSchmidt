/*

	CS440 Homework 1

	Allen Wang
	2/27/2016

	Objective: Write a code that implements the double Gram Schmidt procedure to calculate 
	the inverse of a given square matrix.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void inv_double_gs(double* a, int n, double* u, double* b);	

int main(){
/*Setting A to be a matrix of size n = 4*/
	int n = 3;
	int i, j;

//Sample matrix a, matrix u, and matrix b
	// A is multiplied to test eigenvector magnitude differences and lost digits
	double a[3][3] = {
		{1.0 , 1.0 , 5.0},
		{0.0 , 1.0 , 7.0},
		{9.0 , 9.0, 0.0}
	};
	double u[3][3] = {
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0}
	};

	double b[3][3] = {
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0}
	};


	inv_double_gs(&a[0][0], 3, &u[0][0], &b[0][0]);

	printf("Matrices:\n");

	//Print Matrix A
	printf("\nMatrix A:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("\t%e\t", a[i][j]);
		}
		printf("\n");
	}

	//Print Matrix B
	printf("\nMatrix B:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("\t%e\t", b[i][j]);
		}
		printf("\n");
	}

	//Print Matrix C
	printf("\nMatrix U:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("\t%e\t", u[i][j]);
		}
		printf("\n");
	}
}

void inv_double_gs(double* a, int n, double* u, double* b){
	/*need U and U*
	Given A.
	A * G = U
	Obtain G by doing gs to Identity as you would to a copy of A.
	U = A multiplied by G
	Transpose U = U* 
	Multiply G by U* = A inverse. */
	double (*matrixa)[n] = a;
	double (*matrixu)[n] = u;
	double (*matrixb)[n] = b;
	double g[n][n]; //Initial G matrix set to Identity
	double btranspose[n][n];
	double utranspose[n][n];
	double identitymatrix[n][n];
	int i, j, k, element_n, v1, v2;


	for (i = 0; i < n ; i++){//Initialize G as the identity matrix
		g[i][i] = 1;
	}
	
	for (i = 0; i < n ; i++){//Make copy of A transpose set as U, the results of Gram schmidt on A
		for (j = 0; j < n ; j++)
		{
			matrixu[i][j] = matrixa[j][i];
		}
	}


	
	for (v1 = 0; v1 < n; v1++){//Double GramSchmidt on U and identity to get G
		for (v2 = 0; v2 < v1; v2++){//For loop for left orthogonalization
			double dotproduct = 0.0;
			for (element_n = 0; element_n < n; element_n++){
				dotproduct = dotproduct + (matrixu[v2][element_n] * matrixu[v1][element_n]);
			}
			for (element_n = 0; element_n < n; element_n++){
				matrixu[v1][element_n] = matrixu[v1][element_n] - (dotproduct*matrixu[v2][element_n]);
				g[v1][element_n] = g[v1][element_n] - (dotproduct*g[v2][element_n]);
			}
		}
		for (v2 = v1 + 1; v2 < n; v2++){//For loop for right orthogonalization
			double dotproduct = 0.0;
			for (element_n = 0; element_n < n; element_n++){
				dotproduct = dotproduct + (matrixu[v1][element_n] * matrixu[v2][element_n]);
			}
			for (element_n = 0; element_n < n; element_n++){
				matrixu[v2][element_n] = matrixu[v2][element_n] - (dotproduct*matrixu[v1][element_n]);
				g[v2][element_n] = g[v2][element_n] - dotproduct*g[v1][element_n];
			}
		}

		double magnitude = 0.0;
		for (i = 0; i < n; i++){//normalization to magnitude of U
			magnitude += matrixu[v1][i]*matrixu[v1][i];
		}
		magnitude = sqrt(magnitude);
		for (i = 0; i < n; i++){//Divide by magnitude
			matrixu[v1][i] = matrixu[v1][i]/magnitude;
			g[v1][i] = g[v1][i]/magnitude;
		}

	}
	//G is now complete
	
		
	//Transpose U
	for (i = 0; i < n ; i++){
		for (j = 0; j < n ;j++){
			utranspose[i][j] = matrixu[j][i];
		}
	}
	for (i = 0; i < n ; i++){
		for (j = 0; j < n ;j++){
			matrixu[i][j] = utranspose[i][j];
		}
	}
	

	//Multiply new G by U transpose to equal B (A inverse).
	for (i = 0; i < n ; i++){	
		for (j = 0; j < n ; j++){
			double multiplication = 0.0;
			for (k = 0 ; k < n; k++){
				multiplication = multiplication + (g[k][j]*utranspose[i][k]);
			}
			matrixb[i][j] = multiplication;
		}

	}
	//Transpose B back to get A inverse
	for (i = 0; i < n ; i++){
		for (j = 0; j < n ;j++){
			btranspose[i][j] = matrixb[j][i];
		}
	}
	for (i = 0; i < n ; i++){
		for (j = 0; j < n ;j++){
			matrixb[i][j] = btranspose[i][j];
		}
	}

	//TEST identity

	//Get U transpose again

	for (i = 0; i < n ; i++){
		for (j = 0; j < n ;j++){
			utranspose[i][j] = matrixu[j][i];
		}
	}

	for (i = 0; i < n ; i++){	
		for (j = 0; j < n ; j++){
			double multiplication = 0.0;
			for (k = 0 ; k < n; k++){
				multiplication = multiplication + (matrixu[k][j]*utranspose[i][k]);
			}
			identitymatrix[i][j] = multiplication;
		}

	}

	
}