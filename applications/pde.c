/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: primjena GS,SOR i CG metode na računanje parcijalne diferencijalne jednadžbe
 *
 *
 * BUILD: gcc pde.c -o pde_exe -lblas -llapack -lm
 * RUN  : ./pde_exe
 *
 ====================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../lib/f2c.h"
#include "../lib/fblaswr.h"
#include "../lib/clapack.h"
#include <math.h>
#include "sor.h"
#include "cg.h"
#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define CYN   "\x1B[36m"
#define RESET "\x1B[0m"

void ispis_matrice( doublereal *A, integer n ){
	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++)
			printf(" %f ",A[i+j*n]);
		printf("\n");
	}
	printf("\n\n\n");	
}



void pretty_printer(doublereal* A, integer n) {
    int i,j, max=A[0];
    for(i=1;i<n*n;i++)
		if(A[i]>max) 
			max=A[i];
	doublereal step=max/4;
	 
	for(int k=0;k<n;k++){
		for(j=0;j<n;j++){
			i=k+j*n;
			if(A[i]>=(max-step)) printf("\x1B[31m¤ \x1B[0m");
			if(A[i]<(max-step) && A[i]>=(max-2*step)) printf("\x1B[33m¤ \x1B[0m");
			if(A[i]<(max-2*step) && A[i]>=(max-3*step)) printf("\x1B[32m¤ \x1B[0m");
			if(A[i]<(max-3*step) && A[i]>=(max-4*step)) printf("\x1B[34m¤ \x1B[0m");
			if(A[i]<(max-4*step)) printf("\x1B[36m¤ \x1B[0m");
			 			
		}

			//printf(" %f ",A[i+j*n]);



		printf("\n");
	}
	printf("\n\n\n");
}

int main() {

    integer n=81,N=n*n,i;
    doublereal* A=(doublereal*)malloc(N*sizeof(doublereal));
    doublereal* b=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal* x=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal tol = 1e-8,omega;

    //SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
    char uplo;
    doublereal alpha=0.0,beta=400.0;
    dlaset_(&uplo,&n,&n,&alpha,&beta,A,&n);

    for(i=0;i<n-1;i++)
	if((i+1)%9!=0)
	    A[i+1+n*i]=-100.0;

    for(i=1;i<n;i++)
	if(i%9!=0)
	    A[i-1+n*i]=-100.0;
	
    for(i=0;i<(n-9);i++) A[i+9+n*i]=-100.0;
    for(i=9;i<n;i++) A[i-9+n*i]=-100.0;

    for(i=0;i<n;i++){
	b[i]=0.0;
	if(i==40) b[i]=10000;
	x[i]=0.0;
    }


    //Gauss-Seidel
    omega=1.0;
    integer iterGS=sor(A,b,x,n,omega,tol);
    printf("\n\nGauss-Seidel\n");
    printf("\n%ld\n",iterGS);
    //ispis_matrice(x, sqrt(n));
    pretty_printer(x, sqrt(n));

    //SOR(omega)
    omega=1.53;
    for(i=0;i<n;i++) x[i]=0.0;
    integer iterSOR=sor(A,b,x,n,omega,tol);
    printf("\n\nSOR\n");
    printf("\n%ld\n",iterSOR);
	//ispis_matrice(x, sqrt(n));
	pretty_printer(x, sqrt(n));

    //CG
    for(i=0;i<n;i++) x[i]=0.0;
    integer iterCG=cg(A,x,b,tol,n);
    printf("\n\nCG\n");
    printf("\n%ld\n",iterCG);
	//ispis_matrice(x, sqrt(n));
	pretty_printer(x, sqrt(n));

    free(A);free(b);free(x);
    return 0;
}
