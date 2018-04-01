/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: metoda konjugiranih gradijenata - primijena na Stieltjesovu matricu
 *
 *
 * BUILD: gcc pcg2.c -o pcg2_exe -lblas -llapack -lm
 * RUN  : ./pcg2_exe
 *
 ====================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>


integer cg(doublereal* A, doublereal* x0, doublereal* b, doublereal tol, integer n){
    doublereal* r=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal* d=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal* pom=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal norma_b, broj_alpha, broj_beta, alpha=1.0, beta=0.0,a;
    char trans='N';
    integer inc=1, k=0;
    dcopy_(&n, b, &inc, r, &inc);
    dcopy_(&n, b, &inc, d, &inc);
    norma_b=dnrm2_(&n, b, &inc);
    doublereal norma_r=dnrm2_(&n,r,&inc);
    while(tol<norma_r/norma_b){

	a=ddot_(&n,r,&inc,r,&inc);
	trans='N'; alpha=1.0; beta=0.0;
	dgemv_(&trans,&n,&n,&alpha,A,&n,d,&inc,&beta,pom,&inc); // alpha*A*d = pom

	broj_alpha=a/ddot_(&n,d,&inc,pom,&inc);        

        daxpy_(&n, &broj_alpha, d, &inc, x0, &inc); //x0 = x0 + broj_alpha * d

        broj_alpha*=-1; 
	daxpy_(&n, &broj_alpha, pom, &inc, r, &inc); //r = r - broj_alpha * A * d= r - broj_alpha * pom        

        broj_beta = ddot_(&n, r, &inc, r, &inc)/a; //beta = r^T * r / a
        
	alpha=1.0;
	dscal_(&n, &broj_beta, d, &inc); //d = broj_beta * d
        daxpy_(&n, &alpha, r, &inc, d, &inc); //d = alpha * r + d;
	norma_r=dnrm2_(&n,r,&inc);
        k++;
    }
    return k;
}

int main() {
    FILE* f;
    doublereal* A;
    integer i, n=100;
    A=(doublereal*)malloc(n*n*sizeof(doublereal));
    f=fopen("stieltjes_matr.txt","r");
    for (i=0;i<n*n;i++){
        fscanf(f,"%lf",A+i);
    }
    fclose(f);

    doublereal* x0=(doublereal*)malloc(n*sizeof(doublereal));
doublereal* x1=(doublereal*)malloc(n*sizeof(doublereal));
    for(i=0;i<n;i++) x0[i]=1.0;
    doublereal* b=(doublereal*)malloc(n*sizeof(doublereal));

    //trazimo b tako da rjese bude 1
    // DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    // y := alpha*A'*x + beta*y
    char trans='N';
    doublereal alpha=1.0, beta=0.0;
    integer inc = 1;
    dgemv_(&trans, &n, &n, &alpha, A, &n, x0, &inc, &beta, b, &inc);

	dcopy_(&n,x0,&inc,x1,&inc);


    //postavimo x0 nazad na nul vektor
    //DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
    char uplo='p';
    alpha=0.0;
    beta=0.0;
    integer ldx;
    dlaset_(&uplo, &n, &inc, &alpha, &beta, x0, &ldx);

    //potrebne vrijednosti
    doublereal tol=1.e-8;

    integer broj_iteracija=cg(A,x0,b,tol,n);

    printf("\n potreban broj iteracija je: %ld \n",broj_iteracija);
    printf("\n");
    for(i=0;i<n;i++) printf(" %.8lf ",x0[i]);
    printf("\n");

    free(A); free(b); free(x0);
    return 0;
}
