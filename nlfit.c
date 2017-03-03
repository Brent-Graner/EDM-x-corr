#include <ansi_c.h>

#include "nrutil.h" 


#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

extern int  aph; 
extern int fit;
extern double chi[4];

void funcs(double x, double a[], double *y, double dyda[], int na);


void NLF(double data_x[], double data_y[], int AInumScans, double coefs[],int vary[], int ncoef, 
         double errors[]);

// Numerical Recipes Levenberg-Marquardt Routine function declarations

void mrqmin(double x[], double y[], double sig[], int ndata, double a[], 
    int ia[],
    int ma, double **covar, double **alpha, double *chisq, 
    void (*funcs)(double, double [], double *, double [], int), 
    double *alamda);

void mrqcof(double x[], double y[], double sig[], int ndata, double a[], 
     int ia[],
     int ma, double **alpha, double beta[], double *chisq, 
     void (*funcs)(double, double [], double *, double [], int));

void covsrt(double **covar, int ma, int ia[], int mfit);

void gaussj(double **a, int n, double **b, int m);
 



void NLF(double data_x[], double data_y[], double data_s[], int npoints, double coefs[],int vary[], int ncoef, 
         double errors[])
{
    int i,n;
    int maxiter, * ia;    
    double chilim,chan;  

    double ** covar, ** alpha;
    double  chisq,lamda, prevchi;
  
    
  	chilim= 0.01;      // limit for chi
  	maxiter=20;       // max number of iterations

	ia=ivector(1,ncoef);
	for (i=1;i<=ncoef;i++) ia[i]=vary[i-1];

    alpha=dmatrix(1,ncoef,1,ncoef);
    covar=dmatrix(1,ncoef,1,ncoef); 
    lamda=-1.0;								  // initialize fit 
    i=0;
    n=0;
    
    
    mrqmin(data_x,data_y, data_s, npoints, coefs, ia,
           ncoef, covar, alpha, &chisq, &funcs, &lamda);
      
    prevchi=chisq;  
             
    do 	
    {   	i++;
       	mrqmin(data_x,data_y, data_s, npoints, coefs, ia,
       	ncoef, covar, alpha, &chisq, &funcs, &lamda);
       	
       	chan=fabs(prevchi/chisq-1.0);
       	if (chan < chilim) n++;
       	else n=0;
       	prevchi=chisq; 

    }
    while (n<3 && i<maxiter) ;  
		

    lamda=0.0;
    mrqmin(data_x,data_y, data_s, npoints, coefs, ia,
      ncoef, covar, alpha, &chisq, &funcs, &lamda);

    for (i=1;i<=ncoef;i++) 		
       {
       errors[i-1]=sqrt(chisq*covar[i][i]/npoints);
       }
         	     
     	
     free_dmatrix(alpha,1,ncoef,1,ncoef); 
     free_dmatrix(covar,1,ncoef,1,ncoef);  
	 free_dvector(ia,1,ncoef); 
}   
       



// The actual Numerical Recipes Levenberg-Marquardt Routine 

void mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
    int ma, double **covar, double **alpha, double *chisq, 
    void (*funcs)(double, double [],double *, double [], int), double *alamda) 

{
      int j,k,l; 
      static int mfit; 
      static double ochisq,*atry,*beta,*da,**oneda;

      if (*alamda < 0.0) {   /* Initialization.  */
	  atry=dvector(0,ma-1); 
	  beta=dvector(1,ma); 
	  da=dvector(1,ma); 
	  for (mfit=0,j=1;j<=ma;j++) 
		    if (ia[j]) mfit++; 
	  oneda=dmatrix(1,mfit,1,1); 
	  *alamda=0.001; 
	  mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs); 
	  ochisq=(*chisq); 
	  for (j=1;j<=ma;j++) atry[j-1]=a[j-1];
      } 
      for (j=1;j<=mfit;j++) {
     /* Alter linearized fitting matrix, by augmenting diagonal elements. */
	   for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k]; 
	   covar[j][j]=alpha[j][j]*(1.0+(*alamda)); 
	   oneda[j][1]=beta[j]; 
      } 
      gaussj(covar,mfit,oneda,1); /* Matrix solution.   */
      for (j=1;j<=mfit;j++) da[j]=oneda[j][1]; 
      if (*alamda == 0.0) { /* Once converged, evaluate covariance matrix. */
	     covsrt(covar,ma,ia,mfit); 
	     free_dmatrix(oneda,1,mfit,1,1); 
	     free_dvector(da,1,ma); 
	     free_dvector(beta,1,ma); 
	     free_dvector(atry,1,ma); 
	     return;
      }
      for (j=0,l=1;l<=ma;l++) /* Did the trial succeed? */
	    if (ia[l]) atry[l-1]=a[l-1]+da[++j]; 
      mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs); 
      if (*chisq < ochisq) {  /* Success, accept the new solution.  */
	     *alamda *= 0.3; 
	      ochisq=(*chisq); 
	      for (j=1;j<=mfit;j++) { 
		   for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k]; 
		   beta[j]=da[j];
	      }
	      for (l=1;l<=ma;l++) a[l-1]=atry[l-1]; 
       } else { /* Failure, increase alamda and return.   */
       *alamda *= 3.0; 
       *chisq=ochisq;
       }
}


void mrqcof(double x[], double y[], double sig[], int ndata, double a[], int ia[],
     int ma, double **alpha, double beta[], double *chisq, 
     void (*funcs)(double, double [], double *, double [], int)) 
/* Used by mrqmin to evaluate the linearized fitting matrix alpha, and
vector beta as in (15.5.8), and calculate chisq. */
{
    int i,j,k,l,m,mfit=0; 
    double ymod,wt,sig2i,dy,*dyda;

    dyda=dvector(1,ma); 
    for (j=1;j<=ma;j++) 
	if (ia[j]) mfit++; 
    for (j=1;j<=mfit;j++) { /* Initialize (symmetric) alpha, beta.*/
	 	for (k=1;k<=j;k++) alpha[j][k]=0.0; 
	 	beta[j]=0.0; 
    } 
    *chisq=0.0; 
    for (i=1;i<=ndata;i++) { /* Summation loop over all data. */
       	(*funcs)(x[i-1],a,&ymod,dyda,ma); 
       	sig2i=1.0/(sig[i]*sig[i]); 
       	dy=y[i-1]-ymod; 
    	for (j=0,l=1;l<=ma;l++) { 
	   		if (ia[l]) { 
				wt=dyda[l]*sig2i; 
		  		for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) alpha[j][++k] += wt*dyda[m]; 
		  		beta[j] += dy*wt;
	      	}
		}
		//*chisq += dy*dy*sig2i; /* And find chisq */
		*chisq += dy*dy*sig2i/(a[1]*a[1]+1); /* And find chisq, but divide by (slope^2+1) to give transverse deviation */
    }
    for (j=2;j<=mfit;j++) /* Fill in the symmetric side. */
	for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
    free_dvector(dyda,1,ma); 
}



void covsrt(double **covar, int ma, int ia[], int mfit) 

/* Expand in storage the covariance matrix covar, so as to take into
account parameters that are being held fixed. (For the latter, return
zero covariances.)        */
{
     int i,j,k; 
     double swap;

	for (i=mfit+1;i<=ma;i++) 
	     for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0; 
	k=mfit; 
	for (j=ma;j>=1;j--) { 
	     if (ia[j]) { 
		  for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j]) 
		  for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i]) 
		   k--; 
	     }
	}
 }

void gaussj(double **a, int n, double **b, int m)

 /* Linear equation solution by Gauss-Jordan elimination, equation (2.1.1)
 above. a[1..n][1..n] is the input matrix. b[1..n][1..m] is input containing
 the m right-hand side vectors. On output, a is replaced by its matrix
 inverse, and b is replaced by the corresponding set of solution vectors. */

{
    int *indxc,*indxr,*ipiv;
    int i,icol,irow,j,k,l,ll; 
    double big,dum,pivinv,swap;

    indxc=ivector(1,n); 
 /* The integer arrays ipiv, indxr, andindxc are used for bookkeeping on
    the pivoting.   */
    indxr=ivector(1,n);
    ipiv=ivector(1,n); 
    for (j=1;j<=n;j++) ipiv[j]=0;  
    for (i=1;i<=n;i++) {
	/* This is the main loop over the columns to be reduced.  */
       big=0.0; 
       for (j=1;j<=n;j++)
      /* This is the outer loop of the search for a pivot element. */
	   if (ipiv[j] != 1) 
	       for (k=1;k<=n;k++) { 
		  if (ipiv[k] == 0) { 
			if (fabs(a[j][k]) >= big) {
			    big=fabs(a[j][k]); 
			    irow=j; 
			    icol=k; 
			}
		  } else if (ipiv[k] > 1) ResetTextBox (aph, AP_Error, "gaussj: Singular Matrix-1\n");
	       }
       ++(ipiv[icol]); 

/* We now have the pivot element, so we interchange rows, if needed, to put
the pivot element on the diagonal. The columns are not physically interchanged,
only relabeled: indxc[i], the column of the ith pivot element, is the ith
column that is reduced, while indxr[i] is the row in which that pivot element
was originally located. If indxr[i] 6 = indxc[i] there is an implied column
interchange. With this form of bookkeeping, the solution b's will end up
in the correct order, and the inverse matrix will be scrambled by columns. */

	if (irow != icol) { 
	     for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l]) 
	     for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l]) 
	}
	indxr[i]=irow;

/* We are now ready to divide the pivot row by the pivot element, located at
	irow and icol. */

	indxc[i]=icol;
	if (a[icol][icol] == 0.0) ResetTextBox (aph, AP_Error, "gaussj: Singular Matrix--2\n"); 
	pivinv=1.0/a[icol][icol]; 
	a[icol][icol]=1.0; 
	for (l=1;l<=n;l++) a[icol][l] *= pivinv; 
	for (l=1;l<=m;l++) b[icol][l] *= pivinv;
	for (ll=1;ll<=n;ll++)  /* Next, we reduce the rows...    */
	      if (ll != icol) { /* ...except for the pivot one, of course. */
		   dum=a[ll][icol]; 
		   a[ll][icol]=0.0; 
		   for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum; 
		   for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
	      }
    }

 /* This is the end of the main loop over columns of the reduction. It only
    remains to unscram- ble the solution in view of the column interchanges.
    We do this by interchanging pairs of columns in the reverse order that
    the permutation was built up.   */

	for (l=n;l>=1;l--) {
	   if (indxr[l] != indxc[l]) 
		 for (k=1;k<=n;k++) SWAP(a[k][indxr[l]],a[k][indxc[l]])
    } /* And we are done.  */
    free_ivector(ipiv,1,n);
    free_ivector(indxr,1,n); 
    free_ivector(indxc,1,n); 
}


// Calculates the function and derivatives 

void funcs(double x, double a[], double *y, double dyda[], int na) 

{
/*   a[0]=A, a[1]=B; f=A+Bx   */
    
    *y=a[0]*x+a[1];
    
    dyda[1]=ace*x;
    dyda[2]=ase*earg/a[1];
    dyda[3]=ace;
    dyda[4]=se;
        
    
}



