/******************
Numerical Recipes linear fit functions: sect 15.3
	changed to double precision,
	data arrays indexed from 0..(ndat-1)
	nrerror is string variable, not function
*******************/
#include <ansi_c.h>
#include <math.h>
#include "nrutil.h"
#define ITMAX 100 		// Maximum allowed number of iterations.
#define EPS 3.0e-8 		// Relative accuracy.
#define FPMIN 1.0e-30 	// Number near the smallest representable
						// floating-point number.

#define CGOLD 0.3819660  // brent def.
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

// mnbrak def.
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
// Here GOLD is the default ratio by which successive intervals are magnified; GLIMIT is the
// maximum magnification allowed for a parabolic-fit step.

// fitexy
#define POTN 1.571000
#define BIG 1.0e30
#define PI 3.14159265
#define ACC 1.0e-3
static int nn;
static double *xx,*yy,*sx,*sy,*ww,aa,offs;

extern int standardfit;

char nrerror[100];

double gammq(double a, double x);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);
void fit(double x[], double y[], int ndata, double sig[], double *a,
			double *b, double *siga, double *sigb, double *chi2, double *q);
void avevar(double data[], unsigned long n, double *ave, double *var);
double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
			double (*func)(double));
double zbrent(double (*func)(double), double x1, double x2, double tol);
double chixy(double bang);
void fitexy(double x[], double y[], int ndat, double sigx[], double sigy[],
			double *a, double *b, double *siga, double *sigb, double *chi2, double *q);


double gammq(double a, double x)
/* Returns the incomplete gamma function Q(a; x) = 1 - P(a; x). */
{
	double gamser,gammcf,gln;
	
	if (x < 0.0 || a <= 0.0) strcpy(nrerror,"Invalid arguments in routine gammq");
	if (x < (a+1.0)) { 		// Use the series representation
		gser(&gamser,a,x,&gln);
		return 1.0-gamser; // and take its complement.
		}
	else { 		// Use the continued fraction representation.
		gcf(&gammcf,a,x,&gln);
		return gammcf;
		}
}


void gcf(double *gammcf, double a, double x, double *gln)
/* Returns the incomplete gamma function Q(a; x) evaluated by its continued fraction represen-
tation as gammcf. Also returns ln Gamma(a) as gln. */
{
	double gammln(double xx);
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;			// Set up for evaluating continued fraction
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {	// Iterate to convergence.
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
		}
	if (i > ITMAX) strcpy(nrerror,"a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;		// Put factors in front.
}


void gser(double *gamser, double a, double x, double *gln)
/* Returns the incomplete gamma function P(a; x) evaluated by its series representation as gamser.
Also returns ln Gamma(a) as gln */
{
	double gammln(double xx);
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) strcpy(nrerror,"x less than 0 in routine gser");
		*gamser=0.0;
		return;
		}
	else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
				}
			}
		strcpy(nrerror,"a too large, ITMAX too small in routine gser");
		return;
		}
}


double gammln(double xx)
/* Returns the value ln[Gamma(xx)] for xx>0 */
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
							24.01409824083091,-1.231739572450155,
							0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


void fit(double x[], double y[], int ndata, double sig[], double *a,
			double *b, double *siga, double *sigb, double *chi2, double *q)
/* Given a set of data points x[1..ndata],y[1..ndata] with individual standard deviations
sig[1..ndata], fit them to a straight line y = a + bx by minimizing chi-squre. Returned are
a,b and their respective probable uncertainties siga and sigb, the chi-square chi2, andthe
goodness-of-fit probability q (that the fit would have chi-squre this large or larger). If mwt=0 on
input, then the standard deviations are assumed to be unavailable: q is returned as 1.0 and
the normalization of chi2 is to unit standard deviation on all points. */
{
	int i,mwt=0;
	double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

	if(sig[0]>0) mwt=1;
	*b=0.0;
	if (mwt) { 	// Accumulate sums ...
		ss=0.0;
		for (i=0;i<ndata;i++) { 	// ...with weights
			if(sig[i]==0)
				printf("%d\n",i);
			wt=1.0/SQR(sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
			}
		} 
	else {
		for (i=0;i<ndata;i++) { 	// ...or without weights.
			sx += x[i];
			sy += y[i];
			}
		ss=ndata;
		}
	sxoss=sx/ss;
	if (mwt) {
		for (i=0;i<ndata;i++) {
			t=(x[i]-sxoss)/sig[i];
			st2 += t*t;
			*b += t*y[i]/sig[i];
			}
		}
	else {
		for (i=0;i<ndata;i++) {
			t=x[i]-sxoss;
			st2 += t*t;
			*b += t*y[i];
			}
		}
	*b /= st2; 		// Solve for a, b, .a, and.b.
	*a=(sy-sx*(*b))/ss;
	*siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
	*sigb=sqrt(1.0/st2);
	*chi2=0.0; 		// Calculate chi-squre.
	*q=1.0;
	if (mwt == 0) {
		for (i=0;i<ndata;i++)
			*chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
		sigdat=sqrt((*chi2)/(ndata-2)); 	// For unweighted data evaluate typ-
		*siga *= sigdat;					//		ical sig using chi2, and ad-
		*sigb *= sigdat;					//		just the standard deviations.
		}
	else {
		for (i=0;i<ndata;i++)
			*chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
		if (ndata>2) *q=gammq(0.5*(ndata-2),0.5*(*chi2)); // Equation (15.2.12).
	}
}


void avevar(double data[], unsigned long n, double *ave, double *var)
/* Given array data[1..n], returns its mean as ave and its variance as var. */
{
	unsigned long j;
	double s,ep;
	
	for (*ave=0.0,j=0;j<n;j++) *ave += data[j];
	*ave /= n;
	*var=ep=0.0;
	for (j=0;j<n;j++) {
		s=data[j]-(*ave);
		ep += s;
		*var += s*s;
		}
	*var=(*var-ep*ep/n)/(n-1); 		// Corrected two-pass formula (14.1.8).
}


double brent(double ax, double bx, double cx, double (*f)(double), double tol,
			double *xmin)
/* Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
between ax and cx, andf(bx) is less than both f(ax) and f(cx)), this routine isolates
the minimum to a fractional precision of about tol using Brent's method. The abscissa of
the minimum is returned as xmin, and the minimum function value is returned as brent, the
returned function value. */
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0; 	// This will be the distance moved on
					// the step before last.
	a=(ax < cx ? ax : cx); 	// a and b must be in ascending order,
	b=(ax > cx ? ax : cx);	// but input abscissas need not be. 
	x=w=v=bx; 				// Initializations...
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++) { 		// Main program loop.
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		
		printf("%d  %g   %g\n",iter,x,fx);
		
		
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) { 	// Test for done here.
			*xmin=x;
			return fx;
			}
		if (fabs(e) > tol1) { 		// Construct a trial parabolic fit.
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			/* The above conditions determine the acceptability of the parabolic fit. Here we
			   take the golden section step into the larger of the two segments. */
			else {
				d=p/q; 		// Take the parabolic step.
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
				}
			}
		else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
			}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		// This is the one function evaluation per iteration.
		if (fu <= fx) { 	// Now decide what to do with our function evaluation.
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u) 	// Housekeeping follows:
			SHFT(fv,fw,fx,fu)
			}
		else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
				}
			else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
				}
			} 		// Done with housekeeping.. 
		}			// Back for another iteration
	strcpy(nrerror,"Too many iterations in brent");
	*xmin=x; 		// Never get here.
	return fx;
}


void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
			double (*func)(double))
/* Given a function func, and given distinct initial points ax and bx, this routine searches in
the downhill direction (defined by the function as evaluated at the initial points) and returns
new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
values at the three points, fa, fb, and fc. */
{
	double ulim,u,r,q,fu,dum;
	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) { 					// Switch roles of a and b so that we can go
		SHFT(dum,*ax,*bx,dum)			// downhill in the direction from a to b. 
		SHFT(dum,*fb,*fa,dum)
		}
	*cx=(*bx)+GOLD*(*bx-*ax); 			// First guess for c.
	*fc=(*func)(*cx);
	while (*fb > *fc) { 					// Keep returning here until we bracket.
		r=(*bx-*ax)*(*fb-*fc);				// Compute u by parabolic extrapolation from
		q=(*bx-*cx)*(*fb-*fa);				// a; b; c. TINY is used to prevent any pos-
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/	// 		sible division by zero.
	 		(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		// We won't go farther than this. Test various possibilities:
		if ((*bx-u)*(u-*cx) > 0.0) { 		// Parabolic u is between b and c: try it.
			fu=(*func)(u);
			if (fu < *fc) { 				// Got a minimum between b and c.
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
				}
			else if (fu > *fb) {			// Got a minimum between between a and u.
				*cx=u;
				*fc=fu;
				return;
				}
			u=(*cx)+GOLD*(*cx-*bx); 		// Parabolic fit was no use. Use default mag-nification. fu=(*func)(u);
			}
		else if ((*cx-u)*(u-ulim) > 0.0) { 		// Parabolic fit is between c and its
			fu=(*func)(u);									// allowed limit. 
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
				}
			}
		else if ((u-ulim)*(ulim-*cx) >= 0.0) { 	// Limit parabolic u to maximum allowed value.
			u=ulim;
			fu=(*func)(u);
			}
		else { 					// Reject parabolic u, use default magnification.
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
			}
		SHFT(*ax,*bx,*cx,u) 	// Eliminate oldest point and continue.
		SHFT(*fa,*fb,*fc,fu)
		}
}


double zbrent(double (*func)(double), double x1, double x2, double tol)
/* Using Brent's method, find the root of a function func known to lie between x1 and x2. The
root, returned as zbrent, will be refined until its accuracy is tol. */
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		strcpy(nrerror,"Root must be bracketed in zbrent");
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a; 		// Rename a, b, c and adjust bounding interval d.
			fc=fa;
			e=d=b-a;
			}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
			}
		tol1=2.0*EPS*fabs(b)+0.5*tol; 	// Convergence check.
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa; 		// Attempt inverse quadratic interpolation.
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
				}
			else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
				}
			if (p > 0.0) q = -q; 		// Check whether in bounds.
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;		// Accept interpolation.
				d=p/q;
				}
			else {
				d=xm; 		// Interpolation failed, use bisection.
				e=d;
				}
			}
		else { 		// Bounds decreasing too slowly, use bisection.
			d=xm;
			e=d;
  	 		}
		a=b; 		// Move last best guess to a.
		fa=fb;
		if (fabs(d) > tol1) 		// Evaluate new trial root.
			b +=d;
		else
			b += SIGN(tol1,xm);
		fb=(*func)(b);
		}
	strcpy(nrerror,"Maximum number of iterations exceeded in zbrent");
	return 0.0; 		// Never get here.
}


double chixy(double bang)
/* Captive function of fitexy, returns the value of (chisquare - offs) for the slope b=tan(bang).
Scaled data and offs are communicated via the global variables. */
{
	int j;
	double ans,avex=0.0,avey=0.0,sumw=0.0,b;

	b=tan(bang);
	for (j=0;j<nn;j++) {
		if(standardfit) ww[j] = SQR(b*sx[j])+SQR(sy[j]);
		else ww[j] = SQR(sx[j])+SQR(sy[j]);
		sumw += (ww[j] = (ww[j] < 1.0/BIG ? BIG : 1.0/ww[j]));  // ww is inverted here
		avex += ww[j]*xx[j];		  // compute weighted averages of x and y
		avey += ww[j]*yy[j];
		}
	avex /= sumw;
	avey /= sumw;
	aa=avey-b*avex;		// intercept for this slope
	for (ans = -offs,j=0;j<nn;j++)
		if (standardfit) ans += ww[j]*SQR(yy[j]-aa-b*xx[j]);
		else ans += ww[j]/(SQR(b)+1)*SQR(yy[j]-aa-b*xx[j]); // trial chi2 for perp. deviation
															  // seems to work but not with x error bars (always a minimum for b=inf. unless sx=0
		
	
	return ans;										   
}


void fitexy(double x[], double y[], int ndat, double sigx[], double sigy[],
			double *a, double *b, double *siga, double *sigb, double *chi2, double *q)
/* Straight-line fit to input data x[1..ndat] and y[1..ndat] with errors in both x and y, the re-
spective standard deviations being the input quantities sigx[1..ndat] and sigy[1..ndat].
Output quantities are a and b such that y = a + bx minimizes chisquare, whose value is returned
as chi2. The chisquare probability is returned as q, a small value indicating a poor fit (sometimes
indicating underestimated errors). Standard errors on a and b are returned as siga and sigb.
These are not meaningful if either (i) the fit is poor, or (ii) b is so large that the data are
consistent with a vertical (innite b) line. If siga and sigb are returned as BIG, then the data
are consistent with all values of b. */
{
	int j, noerrbars=0;
	double swap,amx,amn,varx,vary,ang[7],ch[7],scale,bmn,bmx,d1,d2,r2,
		dum1,dum2,dum3,dum4,dum5;
	double *zero;
	
	xx=malloc(sizeof(double)*ndat);
	yy=malloc(sizeof(double)*ndat);
	sx=malloc(sizeof(double)*ndat);
	sy=malloc(sizeof(double)*ndat);
	ww=malloc(sizeof(double)*ndat);
	avevar(x,ndat,&dum1,&varx);		// Find the x and y variances, and scale 
	avevar(y,ndat,&dum1,&vary); 	// 	the data into the global variables
	scale=sqrt(varx/vary);			// 	for communication with the function chixy. 
	nn=ndat;
	
	if (sigy[0]==1.0) noerrbars=1;
	
	for (j=0;j<ndat;j++) {
		xx[j]=x[j];
		yy[j]=y[j]*scale;
		sx[j]=sigx[j];
		sy[j]=sigy[j]*scale;
		ww[j]=sqrt(SQR(sx[j])+SQR(sy[j])); 		// Use both x and y weights in first trial fit
		}
	
	if (sigx[0] > 0 && !standardfit) {
		zero=calloc(ndat,sizeof(double));
		fitexy(xx,yy,nn,zero,sy,&dum1,b,&dum2,&dum3,&dum4,&dum5);
		for (j=0;j<ndat;j++) {
			xx[j]=x[j];
			yy[j]=y[j]*scale;
			sx[j]=sigx[j];
			sy[j]=sigy[j]*scale;
			ww[j]=sqrt(SQR(sx[j])+SQR(sy[j]));
			}
		}
		
	else fit(xx,yy,nn,ww,&dum1,b,&dum2,&dum3,&dum4,&dum5); 	// Trial fit for b.
	offs=ang[1]=0.0; 		// Construct several angles for reference points, and make b an angle.
	ang[2]=atan(*b);
	ang[4]=0.0;
	ang[5]=ang[2];
	ang[6]=POTN;
	for (j=4;j<=6;j++) ch[j]=chixy(ang[j]);
	mnbrak(&ang[1],&ang[2],&ang[3],&ch[1],&ch[2],&ch[3],chixy);
	// Bracket the chisquare minimum and then locate it with brent.
	*chi2=brent(ang[1],ang[2],ang[3],chixy,ACC,b);
	*chi2=chixy(*b);
	*a=aa;
	*q=gammq(0.5*(nn-2),*chi2*0.5); 		// Compute chisquare probability.
	for (r2=0.0,j=0;j<nn;j++) r2 += ww[j]; // Save the inverse sum of weights at the minimum.
	r2=1.0/r2;
	bmx=BIG; 				// Now, find standard errors for b as
	bmn=BIG;		   		// 		points where delta chisquare = 1. 
	offs=(*chi2)+1.0;
	for (j=1;j<=6;j++) { 		// Go through saved values to bracket
		if (ch[j] > offs) {		// 		the desired roots. Note periodicity in slope angles.
			d1=fabs(ang[j]-(*b));
			while (d1 >= PI) d1 -= PI;
			d2=PI-d1;
			if (ang[j] < *b) {
				swap=d1;
				d1=d2;
				d2=swap;
				}
			if (d1 < bmx) bmx=d1;
			if (d2 < bmn) bmn=d2;
			}
		}
	if (bmx < BIG) { 		// Call zbrent to find the roots.
		bmx=zbrent(chixy,*b,*b+bmx,ACC)-(*b);
		amx=aa-(*a);
		bmn=zbrent(chixy,*b,*b-bmn,ACC)-(*b);
		amn=aa-(*a);
		*sigb=sqrt(0.5*(bmx*bmx+bmn*bmn))/(scale*SQR(cos(*b)));
		*siga=sqrt(0.5*(amx*amx+amn*amn)+r2)/scale; 		// Error in a has additional piece r2.
		}
	else {
		if (noerrbars) {
			ang[1]=sy[0]; ang[2]=0.5*sy[0];
			//mnbrak(&ang[1],&ang[2],&ang[3],&ch[1],&ch[2],&ch[3],fabs(chixy-1.));
			}
	
		(*sigb)=(*siga)=BIG;
		}
	
	*a /= scale; 				// Unscale the answers.
	*b=tan(*b)/scale;
//	free(ww);
//	free(sy);
//	free(sx);
//	free(yy);
//	free(xx);
}
