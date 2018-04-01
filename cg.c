/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: implementacija metode konjugiranih gradijenata - primjena na matricu 100x100 simetricna, pozitivno definitna, sa svojstvenim vrijednostima {1,4,9,...,10000} dobivena kao produkt A=QDQ^t 
 *
 *
 * BUILD: gcc sor2.c -o sor2_exe -lblas -llapack -lm
 * RUN  : ./sor2_exe
 *
 ====================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "clapack.h"
#include "fblaswr.h"
#include <math.h>

doublereal kriterij_zaustavljanja(doublereal* A, doublereal* x, doublereal* b, doublereal norma_b, integer n, doublereal* pom){
    char trans='N';
    integer inc=1;
    doublereal alpha=-1.0, beta=1.0;
    dcopy_(&n, b, &inc, pom, &inc);
    dgemv_(&trans, &n, &n, &alpha, A, &n, x, &inc, &beta, pom, &inc); // pom = alpha * A * x + beta * pom
    return dnrm2_(&n, pom, &inc)/norma_b;
}

integer cg(doublereal* A, doublereal* x, doublereal* b, doublereal tol, integer n){
    doublereal* r=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal* d=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal* pom=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal norma_b, broj_alpha, broj_beta, alpha=1.0, beta=0.0;
    char trans='N';
    integer inc=1, k=0, i;
    dcopy_(&n, b, &inc, r, &inc);
    dcopy_(&n, b, &inc, d, &inc);
    norma_b=dnrm2_(&n, b, &inc);
    while(tol<kriterij_zaustavljanja(A,x,b,norma_b,n,pom)){
	alpha=1.0; beta=0.0;
        dgemv_(&trans, &n, &n, &alpha, A, &n, d, &inc, &beta, pom, &inc); //broj_alpha=r^T * r/ d^T * A * d
        broj_alpha= ddot_(&n, r, &inc, r, &inc)/ddot_(&n, d, &inc, pom, &inc); //broj_alpha=r^T * r/ d^T * A * d
        daxpy_(&n, &broj_alpha, d, &inc, x, &inc); //x = x + broj_alpha * d
        dcopy_(&n, r, &inc, pom, &inc); //spremi r u pom za racunanje broj_beta
        broj_alpha*=-1; beta=1.0;
        dgemv_(&trans, &n, &n, &broj_alpha, A, &n, d, &inc, &beta, r, &inc); //r = r - broj_alpha * A * d
        broj_beta = ddot_(&n, r, &inc, r, &inc)/ddot_(&n, pom, &inc, pom, &inc); //beta = r^T * r / pom^T * pom
        dscal_(&n, &broj_beta, d, &inc); //d = broj_beta * d
        daxpy_(&n, &alpha, r, &inc, d, &inc); //d = alpha * r + d;
        k++;
    }
    k++;
    return k;
}

int main() {
    integer n=100, k=1,i, info, N=n*n,inc=1, idist = 2, iseed[4]={2000, 1000, 2000, 3999};
    char trans = 'N', side='L';
    doublereal tol = 1e-8, alpha=1.0, beta=0.0;

    doublereal* A = (doublereal*)malloc(n*n*sizeof(doublereal));
    doublereal* b = (doublereal*)malloc(n*sizeof(doublereal));
    doublereal* x = (doublereal*)malloc(n*sizeof(doublereal));
    doublereal* D = (doublereal*)malloc(n*n*sizeof(doublereal));
    doublereal* work = (doublereal*)malloc(n*sizeof(doublereal));
    doublereal* tau = (doublereal*)malloc(n*sizeof(doublereal));
    for(i=0;i<n*n;i++) { //generiranje dijagonale svojstvenih vrijednosti
        if(i/n==i%n) {
            D[i]=k*k;
            k++; 
        }        
        else D[i]=0.0;
    }
    dlarnv_(&idist, iseed, &N, A); //generiranje matrice Q=A - uniformna distribucija
    dgeqrf_(&n, &n, A, &n, tau, work, &n, &info); //QRF
    dormqr_(&side, &trans, &n, &n, &n, A, &n, tau, D, &n, work, &n, &info); //D=A*D
    trans='T'; side='R';
    dormqr_(&side, &trans, &n, &n, &n, A, &n, tau, D, &n, work, &n, &info); //A=D*A^t 
    for(i=0;i<n;i++) x[i]=1.0;
    trans='N'; 
    dgemv_(&trans, &n, &n, &alpha, D, &n, x, &inc, &beta, b, &inc); //b=alpha*A*x
    for(i=0;i<n;i++) x[i]=0.0;    
    integer kraj=cg(D,x,b,tol,n);
    printf("broj iteracija je %ld \nx:\n",kraj);
    //for(i=0;i<n;i++) printf("%lf\n",x[i]);
    
    free(A);
    free(b);
    free(x);
    free(D);
    free(tau);
    return 0;
}
