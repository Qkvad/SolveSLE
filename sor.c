/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (Applied Mathematics)
 * Project: SOR method for solving system of linear equations
 *
 *
 * BUILD: gcc sor.c -o sor_exe -lblas -llapack -lm
 * RUN  : ./sor_exe
 *
 ====================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include <math.h>

doublereal A[] = {101, -4, 8, 12, -4, 20, -7, 3, 8, -7, 78, 32, 12, 3, 32, 113};
doublereal b[] = {117, 12, 111, 160};
doublereal x0[]={0, 0, 0, 0};
doublereal omega=1.05;
doublereal epsilon=1e-5;
integer n=4;
integer N=16;

doublereal sor_norma(){
    doublereal* R=(doublereal*)malloc(n*n*sizeof(doublereal));
    doublereal* L=(doublereal*)malloc(n*n*sizeof(doublereal));
    doublereal* D=(doublereal*)malloc(n*n*sizeof(doublereal));
    doublereal* work=(doublereal*)malloc(n*n*sizeof(doublereal));
    char uploL='L',uploR='U',norm='i',trans='N';
    doublereal alpha=1.0;
    dlacpy_(&uploL, &n, &n, A, &n, L, &n);
    dlacpy_(&uploR, &n, &n, A, &n, R, &n);
    for(integer i=0; i<n*n; i++) {
        if (i / n == i % n) D[i] = A[i] ;
        else D[i]=0;
        L[i]*=omega;
        R[i]=-omega*R[i];
    }
    dlacpy_(&uploR, &n, &n, D, &n, L, &n);
    for(integer i=0; i<n*n;i++) D[i]=(1-omega)*D[i];
    dlacpy_(&uploL, &n, &n, D, &n, R, &n);
    dtrsm_(&uploL, &uploL, &trans, &trans, &n, &n, &alpha, L, &n, R, &n);
    return dlange_(&norm, &n, &n, R, &n, work);
}

integer sor_konvergencija(){
    if (sor_norma() < 1.0) return 1;
    return 0;
}

int sor_kriterij(){
    doublereal* x=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal* work=(doublereal*)malloc(n*sizeof(doublereal));
    doublereal pom;
    char norm='m';
    integer inc = 1;
    for(integer i=0;i<n;i++) x[i]=x0[i];
    for(integer i=0;i<n;i++){
        x[i]=(1-omega)*x[i];
        pom=b[i];
        for(integer j=0;j<i;j++) pom=pom-A[i+n*j]*x0[j];
        for(integer j=i+1;j<n;j++) pom=pom-A[i+n*j]*x0[j];
        x[i]=x[i]+pom*omega/A[i*n+i];
    }
    for(integer i=0;i<n;i++) x[i]=x[i]-x0[i];
    doublereal norma_razlike = dlange_(&norm,&inc, &n, x, &n, work);
    doublereal norma = sor_norma();
    return (int)ceil( log(epsilon*(1-norma)/norma_razlike)/log(norma));
}

int sor_rjesavac(doublereal* x){
    int k=1, i, j;
    doublereal pom;
    int kriterij_zaustavljanja = sor_kriterij();
    while(k<=kriterij_zaustavljanja){
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
    }
    for(i=0;i<n;i++) x[i]=x0[i];
    return kriterij_zaustavljanja;
}

int main(int argc, char *argv[]) {
    if(sor_konvergencija() == 1) {
        doublereal *x = (doublereal*) malloc(n * sizeof(doublereal));
        int broj_iteracija = sor_rjesavac(x);
        int i;
        printf("\nx = ( ");
        for (i = 0; i < n; i++) { i != n-1 ? printf("%.16f, ", x[i]) : printf("%.16f )", x[i]); }
        printf("\nbroj iteracija = %d\n", broj_iteracija);
    }
    else
        printf("Nije zadovoljen uvjet konvergencije.");
    return 0;
}
