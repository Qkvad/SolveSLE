/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: prekondicionirana metoda konjugiranih gradijenata - dijagonalna matrica prekondicioniranja
 *
 *
 * BUILD: gcc pcg.c -o pcg_exe -lblas -llapack -lm
 * RUN  : ./pcg_exe
 *
 ====================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

integer pcg (integer n, doublereal* A, doublereal* b, doublereal* x0, doublereal tol){
    integer inc=1, k=0;
    char side='L', uplo='L',trans,diag='N';
    doublereal alphaf,betaf, norm_r, norm_b,kriterij,a,alpha,beta,a_novi;
    doublereal* r=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal* d=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal* p=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal* pom=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal* D=(doublereal*)malloc(n*n*sizeof(doublereal));

    dcopy_(&n, b, &inc, r, &inc);
    norm_b=dnrm2_(&n,b,&inc); //racunanje norme od b
    norm_r=norm_b; //racunanje norme od r
    kriterij=norm_r/norm_b;

    for(integer i=0;i<n;i++) D[i+i*n]=A[i+i*n]; 

    dcopy_(&n,r,&inc,p,&inc);
    for(integer i=0;i<n;i++) p[i]=1/D[i+i*n]*p[i];

    printf("\n");
    for(integer i=0;i<n;i++) printf(" %lf ",p[i]);
    printf("\n");

    dcopy_(&n,p,&inc,d,&inc);

    while(kriterij>tol){

        a=ddot_(&n,r,&inc,p,&inc); //a=rk^t *pk

        alphaf=1.0; betaf=0.0;trans='N';
        dgemv_(&trans,&n,&n,&alphaf,A,&n,d,&inc,&betaf,pom,&inc); //pom=A*dk

        alpha=a/ddot_(&n,d,&inc,pom,&inc);

        daxpy_(&n,&alpha,d,&inc,x0,&inc); //xk+1=xk+alpha*dk

        alpha=-alpha;
        daxpy_(&n,&alpha,pom,&inc,r,&inc); //rk+1=rk-alpha*A*dk=rk-alpha*pom

	dcopy_(&n,r,&inc,p,&inc);
        for(integer i=0;i<n;i++) p[i]=1/D[i+i*n]*p[i];

        a_novi=ddot_(&n,r,&inc,p,&inc); //a=rk+1^t *pk+1
        beta=a_novi/a; //beta = rk+1^t *pk+1 / rk^t *pk

        dscal_(&n,&beta,d,&inc); // beta* d = d
        betaf=1.0;
        daxpy_(&n,&betaf,p,&inc,d,&inc); //beta * d + p = d

        norm_r=dnrm2_(&n,r,&inc);

        kriterij=norm_r/norm_b;
        printf("\n %ld. kriterij = %.9lf ",k,kriterij);

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
    for(i=0;i<n;i++) x0[i]=1.0;
    doublereal* b=(doublereal*)malloc(n*sizeof(doublereal));

    //trazimo b tako da rjese bude 1
    // DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    // y := alpha*A'*x + beta*y
    char trans='N';
    doublereal alpha=1.0, beta=0.0;
    integer inc = 1;
    dgemv_(&trans, &n, &n, &alpha, A, &n, x0, &inc, &beta, b, &inc);

    //postavimo x0 nazad na nul vektor
    //DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
    char uplo='p';
    alpha=0.0;
    beta=0.0;
    integer ldx;
    dlaset_(&uplo, &n, &inc, &alpha, &beta, x0, &ldx);

    //potrebne vrijednosti
    doublereal tol=1.e-8;

    integer broj_iteracija=pcg(n,A,b,x0,tol);

    printf("\n potreban broj iteracija je: %ld \n",broj_iteracija);
    printf("\n");
    for(i=0;i<n;i++) printf(" %lf ",x0[i]);
    printf("\n");

    free(A); free(b); free(x0);
    return 0;
}
