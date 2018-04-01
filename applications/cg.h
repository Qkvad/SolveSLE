/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: metoda konjugiranih gradijenata
 *
 ====================================================================================================================*/

#include <stdlib.h>
#include "../lib/f2c.h"
#include "../lib/fblaswr.h"
#include "../lib/clapack.h"
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


