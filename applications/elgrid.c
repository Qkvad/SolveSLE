/*=====================================================================================================================
 *
 * Author : Petra Brčić
 * College: PMF-MO (DS - Primijenjena Matematika)
 * Project: primjena GS,SOR i GMRES metode na računanje sustava električne mreže
 *
 *
 * BUILD: gcc elgrid.c -o elgrid_exe -lblas -llapack -lm
 * RUN  : ./elgrid_exe
 *
 ====================================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../lib/f2c.h"
#include "../lib/fblaswr.h"
#include "../lib/clapack.h"
#include <math.h>
#include "GMRES.c"
#include "sor.h"

integer n=6;
doublereal A[]={11,-20,0,0,0,-2,-5,41,-3,0,-3,0,0,-15,7,-1,0,0,0,0,-4,2,-10,0,0,-6,0,-1,28,-15,-1,0,0,0,-15,47};

void matvec(doublereal *alpha,doublereal *x,doublereal *beta,doublereal *y){
    char trans='n';
    integer inc=1;
    dgemv_(&trans,&n,&n,alpha,A,&n,x,&inc,beta,y,&inc);
}

void psolve(doublereal *x,doublereal *b){
    integer inc=1;
    dcopy_(&n,b,&inc,x,&inc);
}

int main() {

    integer N=n*n,i;
    doublereal b[] = {500.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    doublereal x[]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    doublereal tol = 1e-8,omega;

    //Gauss-Seidel
    omega=1.0;
    integer iterGS=sor(A,b,x,n,omega,tol);
    printf("\n\nGauss-Seidel\n");
    for(i=0;i<n;i++) printf("%lf ",x[i]);
    printf("\n%ld\n",iterGS);

    //SOR(omega)
    omega=1.35;
    for(i=0;i<n;i++) x[i]=0.0;
    integer iterSOR=sor(A,b,x,n,omega,tol);
    printf("\n\nSOR\n");
    for(i=0;i<n;i++) printf("%lf ",x[i]);
    printf("\n%ld\n",iterSOR);


    //GMRES
    for(i=0;i<n;i++) x[i]=0.0;
    doublereal *work2=malloc(n*(n+4)*sizeof(doublereal));
    doublereal *H=malloc((n+1)*(n+2)*sizeof(doublereal));
    integer restrt=n, ldh=n+1, iter=n, info;
    gmres_( &n, b, x, &restrt, work2, &n, H, &ldh, &iter, &tol, matvec, psolve, &info );
    printf("\n\nGMRES\n");
    for(i=0;i<n;i++) printf("%lf ",x[i]);
    printf("\n%ld\n",iter);
    
    return 0;
}
