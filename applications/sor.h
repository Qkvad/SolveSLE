#include <stdlib.h>
#include "../lib/f2c.h"
#include "../lib/fblaswr.h"
#include "../lib/clapack.h"

int sor(doublereal* A, doublereal* b, doublereal* x0, integer n, doublereal omega, doublereal epsilon){

    integer k=0, i, j, inc=1;
    char trans='N';
    doublereal pom,kriterij,alpha=-1.0,beta=1.0;
    doublereal normab=dnrm2_(&n,b,&inc);
    doublereal* b2=(doublereal*)malloc(n*sizeof(doublereal));
    kriterij=1.0;
    while(kriterij>epsilon){
        k++;
        for(i=0;i<n;i++){
            x0[i]=(1-omega)*x0[i];
            pom=b[i];
            for(j=0;j<i;j++)
                pom=pom-A[i+n*j]*x0[j];
            for(j=i+1;j<n;j++)
                pom=pom-A[i+n*j]*x0[j];
            x0[i]=x0[i]+pom*omega/A[i*n+i];
        }
	dcopy_(&n,b,&inc,b2,&inc);
	dgemv_(&trans,&n,&n,&alpha,A,&n,x0,&inc,&beta,b2,&inc);
	kriterij=dnrm2_(&n,b2,&inc)/normab;
    }
    return k;
}


