/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: primjena GS,SOR i CG metode na računanje obične diferencijalne jednadžbe
 *
 *
 * BUILD: gcc ode.c -o ode_exe -lblas -llapack -lm
 * RUN  : ./ode_exe
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

int main() {

    integer n=99,N=n*n,i;
    doublereal* A=(doublereal*)malloc(N*sizeof(doublereal));
    doublereal* b=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal* x=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal tol = 1e-8,omega;

    for(i=0;i<N;i++) A[i]=0.0;

    for(i=0;i<n;i++){
        A[i+n*i]=1.9999;
	A[i-1+n*i]=-1;
	A[i+1+n*i]=-1;
    }

    doublereal a=0.01;
    for(i=0;i<n;i++){
	b[i]=0.0002*sin(a);
	a+=0.01;
	if(i=98) b[i]+=cos(1);
	x[i]=0;
    }

    //SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
    char uplo='l';
    integer info,inc=1;
    doublereal* pom=(doublereal*)malloc(N*sizeof(doublereal));
    dcopy_(&N, A, &inc, pom, &inc);
    dpotrf_(&uplo,&n,pom,&n,&info);
    if(info==0) printf("Matrica je pozitivno definitna.");
    if(info>0) printf("Matrica nije pozitivno definitna.");


    //Gauss-Seidel
    omega=1.0;
    integer iterGS=sor(A,b,x,n,omega,tol);
    printf("\n\nGauss-Seidel\n");
    printf("\n%ld\n",iterGS);

    //egzaktno rjesenje
    a=0.01;
    for(i=0;i<n;i++) {
        printf("\n%lf     %lf",a*cos(a),x[i]);
        a+=0.01;
    }

    //SOR(omega)
    omega=1.95;
    for(i=0;i<n;i++) x[i]=0.0;
    integer iterSOR=sor(A,b,x,n,omega,tol);
    printf("\n\nSOR\n");
    printf("\n%ld\n",iterSOR);

    //egzaktno rjesenje
    a=0.01;
    for(i=0;i<n;i++) {
        printf("\n%lf     %lf",a*cos(a),x[i]);
        a+=0.01;
    }

    //CG
    for(i=0;i<n;i++) x[i]=0.0;
    integer iterCG=cg(A,x,b,tol,n);
    printf("\n\nCG\n");
    printf("\n%ld\n",iterCG);

    //egzaktno rjesenje
    a=0.01;
    for(i=0;i<n;i++) {
        printf("\n%lf     %lf",a*cos(a),x[i]);
        a+=0.01;
    }
    return 0;
}
