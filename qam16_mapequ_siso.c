#include <math.h>
#include "mex.h"
#include "utilities.h"

// interface
#define y_if      prhs[0]
#define nvar_in   prhs[1]
#define ch_if  	  prhs[2]
#define llrin_if  prhs[3]
#define llrout_if plhs[0]
#define inputs    4
#define outputs   1

#define K10 0.316227766016838

double getProb (double* pd, int i, int symi)
	{
	double d;
	int i0=4*i, i1=4*i+1, i2=4*i+2, i3=4*i+3;
	switch (symi)
		{
		case  0: d= pd[i0]+pd[i1]+pd[i2]+pd[i3]; break;
		case  1: d=-pd[i0]+pd[i1]+pd[i2]+pd[i3]; break;
		case  2: d= pd[i0]-pd[i1]+pd[i2]+pd[i3]; break;
		case  3: d=-pd[i0]-pd[i1]+pd[i2]+pd[i3]; break;
		case  4: d= pd[i0]+pd[i1]-pd[i2]+pd[i3]; break;
		case  5: d=-pd[i0]+pd[i1]-pd[i2]+pd[i3]; break;
		case  6: d= pd[i0]-pd[i1]-pd[i2]+pd[i3]; break;
		case  7: d=-pd[i0]-pd[i1]-pd[i2]+pd[i3]; break;
		case  8: d= pd[i0]+pd[i1]+pd[i2]-pd[i3]; break;
		case  9: d=-pd[i0]+pd[i1]+pd[i2]-pd[i3]; break;
		case 10: d= pd[i0]-pd[i1]+pd[i2]-pd[i3]; break;
		case 11: d=-pd[i0]-pd[i1]+pd[i2]-pd[i3]; break;
		case 12: d= pd[i0]+pd[i1]-pd[i2]-pd[i3]; break;
		case 13: d=-pd[i0]+pd[i1]-pd[i2]-pd[i3]; break;
		case 14: d= pd[i0]-pd[i1]-pd[i2]-pd[i3]; break;
		case 15: d=-pd[i0]-pd[i1]-pd[i2]-pd[i3]; break;
		}
	return 0.5*d;
	}

double getProb3 (double* pd, int i, int symi)
	{
	double d;
	int i0=4*i, i1=4*i+1, i2=4*i+2;
	switch (symi)
		{
		case 0: case  8: d= pd[i0]+pd[i1]+pd[i2]; break;
		case 1: case  9: d=-pd[i0]+pd[i1]+pd[i2]; break;
		case 2: case 10: d= pd[i0]-pd[i1]+pd[i2]; break;
		case 3: case 11: d=-pd[i0]-pd[i1]+pd[i2]; break;
		case 4: case 12: d= pd[i0]+pd[i1]-pd[i2]; break;
		case 5: case 13: d=-pd[i0]+pd[i1]-pd[i2]; break;
		case 6: case 14: d= pd[i0]-pd[i1]-pd[i2]; break;
		case 7: case 15: d=-pd[i0]-pd[i1]-pd[i2]; break;
		}
	return 0.5*d;
	}

double getProb2 (double* pd, int i, int symi)
	{
	double d;
	int i0=4*i, i1=4*i+1, i3=4*i+3;
	switch (symi)
		{
		case  0: case  4: d= pd[i0]+pd[i1]+pd[i3]; break;
		case  1: case  5: d=-pd[i0]+pd[i1]+pd[i3]; break;
		case  2: case  6: d= pd[i0]-pd[i1]+pd[i3]; break;
		case  3: case  7: d=-pd[i0]-pd[i1]+pd[i3]; break;
		case  8: case 12: d= pd[i0]+pd[i1]-pd[i3]; break;
		case  9: case 13: d=-pd[i0]+pd[i1]-pd[i3]; break;
		case 10: case 14: d= pd[i0]-pd[i1]-pd[i3]; break;
		case 11: case 15: d=-pd[i0]-pd[i1]-pd[i3]; break;
		}
	return 0.5*d;
	}

double getProb1 (double* pd, int i, int symi)
	{
	double d;
	int i0=4*i, i2=4*i+2, i3=4*i+3;
	switch (symi)
		{
		case  0: case  2: d= pd[i0]+pd[i2]+pd[i3]; break;
		case  1: case  3: d=-pd[i0]+pd[i2]+pd[i3]; break;
		case  4: case  6: d= pd[i0]-pd[i2]+pd[i3]; break;
		case  5: case  7: d=-pd[i0]-pd[i2]+pd[i3]; break;
		case  8: case 10: d= pd[i0]+pd[i2]-pd[i3]; break;
		case  9: case 11: d=-pd[i0]+pd[i2]-pd[i3]; break;
		case 12: case 14: d= pd[i0]-pd[i2]-pd[i3]; break;
		case 13: case 15: d=-pd[i0]-pd[i2]-pd[i3]; break;
		}
	return 0.5*d;
	}

double getProb0 (double* pd, int i, int symi)
	{
	double d;
	int i1=4*i+1, i2=4*i+2, i3=4*i+3;
	switch (symi)
		{
		case  0: case  1: d=+pd[i1]+pd[i2]+pd[i3]; break;
		case  2: case  3: d=-pd[i1]+pd[i2]+pd[i3]; break;
		case  4: case  5: d=+pd[i1]-pd[i2]+pd[i3]; break;
		case  6: case  7: d=-pd[i1]-pd[i2]+pd[i3]; break;
		case  8: case  9: d=+pd[i1]+pd[i2]-pd[i3]; break;
		case 10: case 11: d=-pd[i1]+pd[i2]-pd[i3]; break;
		case 12: case 13: d=+pd[i1]-pd[i2]-pd[i3]; break;
		case 14: case 15: d=-pd[i1]-pd[i2]-pd[i3]; break;
		}
	return 0.5*d;
	}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{
	int Lc,Ly,M,Ls,i,j,k,l,st2lim,stmax,st2;
	double *yR,*yI,*chR,*chI,*llrout,*llrin,*alpha,*beta,nvar,vv,d,dd,*sig0,*sig1,*sig2,*sig3;
	double symR[16]={3*K10,-3*K10,1*K10,-1*K10,3*K10,-3*K10,1*K10,-1*K10,3*K10,-3*K10,1*K10,-1*K10,3*K10,-3*K10,1*K10,-1*K10};
  	double symI[16]={3*K10,3*K10,3*K10,3*K10,-3*K10,-3*K10,-3*K10,-3*K10,1*K10,1*K10,1*K10,1*K10,-1*K10,-1*K10,-1*K10,-1*K10};
	CV tv;
	CN d0,d1;

	if (!((nrhs==inputs)||(nrhs==inputs-1))||(nlhs!=outputs))
		mexErrMsgTxt("usage: llr_out=qam16_mapequ_siso(y,nvar,h) or llr_out=qam16_mapequ_siso(y,nvar,h,llr_in)\n");

	Lc=mxGetM(ch_if);
	Ly=mxGetM(y_if);
	nvar=mxGetScalar(nvar_in);
	vv=-1.0/nvar;
	M=Lc-1;
	Ls=1<<(4*M);
	stmax=Ls-1;

	yR=mxGetPr(y_if);
	yI=mxGetPi(y_if);
	chR=mxGetPr(ch_if);
	chI=mxGetPi(ch_if);
	if (nrhs==inputs)
		{
		llrin=(double*)mxGetPr(llrin_if);
		for(i=0;i<(Ly*4);i++) if (llrin[i]>LMAX) llrin[i]=LMAX; else if (llrin[i]<(-LMAX)) llrin[i]=-LMAX;
		}
	else
		{
		llrin=(double*)calloc(Ly*4,sizeof(double));
		for(i=0;i<(Ly*4);i++) llrin[i]=0;
		}
	llrout_if=mxCreateDoubleMatrix(Ly*4,1,mxREAL);
	llrout=(double*)mxGetPr(llrout_if);

	sig0=(double*)calloc(16,sizeof(double));
  	sig1=(double*)calloc(16,sizeof(double));
  	sig2=(double*)calloc(16,sizeof(double));
  	sig3=(double*)calloc(16,sizeof(double));
	alpha=(double*)calloc(Ls*Ly,sizeof(double));
	beta=(double*)calloc(Ls*Ly,sizeof(double));
	for(i=0;i<(Ly*Ls);i++) {alpha[i]=-INFINITY; beta[i]=-INFINITY;}
	for(i=0;i<Ls;i++) {beta[i*Ly+Ly-1]=0.0; alpha[i*Ly]=0.0;} // termination

 	CValloc(&tv,Ls<<4);
	for(i=0;i<(Ls<<4);i++)
		{
		d=0; dd=0;
		if (chI!=0) for(l=0;l<=M;l++)
			{
			d +=chR[l]*symR[(i>>(l*4))&15]-chI[l]*symI[(i>>(l*4))&15];
			dd+=chR[l]*symI[(i>>(l*4))&15]+chI[l]*symR[(i>>(l*4))&15];
			}
		else for(l=0;l<=M;l++)
			{
			d +=chR[l]*symR[(i>>(l*4))&15];
			dd+=chR[l]*symI[(i>>(l*4))&15];
			}
		tv.c[i].r=d; tv.c[i].i=dd;
		}

	if (Lc>1)
		{
		// obtain alphas
		for(i=0;i<Ly-1;i++)
			{
			for(j=0;j<Ls;j++)
				{
				st2=j<<4;
				d0.r=yR[i];
				d0.i=yI[i];
				st2lim=st2&stmax;
				for(k=0;k<16;k++)
					{
					CNsub(&tv.c[st2+k],&d0,&d1);
	        		addplogs(&alpha[i+1+(st2lim+k)*Ly],alpha[i+j*Ly]+vv*CNabs(&d1)+getProb(llrin,i,k));
					}
				}
			d=-INFINITY;
			for(l=0;l<Ls;l++) if (alpha[i+1+l*Ly]>d) d=alpha[i+1+l*Ly];
			for(l=0;l<Ls;l++) alpha[i+1+l*Ly]-=d;
			}

		// obtain betas
		for(i=Ly-2;i>=0;i--)
			{
			for(j=0;j<Ls;j++)
				{
				st2=j<<4;
				d0.r=yR[i+1];
				d0.i=yI[i+1];
				st2lim=st2&stmax;
				for(k=0;k<16;k++)
	        		{
	        		CNsub(&tv.c[st2+k],&d0,&d1);
	        		addplogs(&beta[i+j*Ly],beta[i+1+(st2lim+k)*Ly]+vv*CNabs(&d1)+getProb(llrin,i+1,k));
	        		}
				}
			d=-INFINITY;
			for(l=0;l<Ls;l++) if (beta[i+l*Ly]>d) d=beta[i+l*Ly];
			for(l=0;l<Ls;l++) beta[i+l*Ly]-=d;
			}

		// obtain output
		for(i=0;i<Ly;i++)
			{
			for(k=0;k<16;k++) {sig0[k]=-INFINITY; sig1[k]=-INFINITY; sig2[k]=-INFINITY; sig3[k]=-INFINITY;}
			for(j=0;j<Ls;j++)
				{
				st2=j<<4;
				d0.r=yR[i];
				d0.i=yI[i];
				st2lim=st2&stmax;
				for(k=0;k<16;k++)
			        {
					CNsub(&tv.c[st2+k],&d0,&d1); d=vv*CNabs(&d1);
			        d+=alpha[i+j*Ly]+beta[i+(st2lim+k)*Ly];
			        addplogs(&sig0[k],d+getProb0(llrin,i,k));
			        addplogs(&sig1[k],d+getProb1(llrin,i,k));
			        addplogs(&sig2[k],d+getProb2(llrin,i,k));
			        addplogs(&sig3[k],d+getProb3(llrin,i,k));
			        }
				}
		    k=4*i; d=-INFINITY; dd=-INFINITY;
			addplogs(&d ,sig0[ 0]); addplogs(&d ,sig0[ 2]); addplogs(&d ,sig0[ 4]); addplogs(&d ,sig0[ 6]);
		    addplogs(&d ,sig0[ 8]); addplogs(&d ,sig0[10]); addplogs(&d ,sig0[12]); addplogs(&d ,sig0[14]);
		    addplogs(&dd,sig0[ 1]); addplogs(&dd,sig0[ 3]); addplogs(&dd,sig0[ 5]); addplogs(&dd,sig0[ 7]);
		    addplogs(&dd,sig0[ 9]); addplogs(&dd,sig0[11]); addplogs(&dd,sig0[13]); addplogs(&dd,sig0[15]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
		    k++; d=-INFINITY; dd=-INFINITY;
		    addplogs(&d ,sig1[ 0]); addplogs(&d ,sig1[ 1]); addplogs(&d ,sig1[ 4]); addplogs(&d ,sig1[ 5]);
		    addplogs(&d ,sig1[ 8]); addplogs(&d ,sig1[ 9]); addplogs(&d ,sig1[12]); addplogs(&d ,sig1[13]);
		    addplogs(&dd,sig1[ 2]); addplogs(&dd,sig1[ 3]); addplogs(&dd,sig1[ 6]); addplogs(&dd,sig1[ 7]);
		    addplogs(&dd,sig1[10]); addplogs(&dd,sig1[11]); addplogs(&dd,sig1[14]); addplogs(&dd,sig1[15]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
		    k++; d=-INFINITY; dd=-INFINITY;
		    addplogs(&d ,sig2[ 0]); addplogs(&d ,sig2[ 1]); addplogs(&d ,sig2[ 2]); addplogs(&d ,sig2[ 3]);
		    addplogs(&d ,sig2[ 8]); addplogs(&d ,sig2[ 9]); addplogs(&d ,sig2[10]); addplogs(&d ,sig2[11]);
		    addplogs(&dd,sig2[ 4]); addplogs(&dd,sig2[ 5]); addplogs(&dd,sig2[ 6]); addplogs(&dd,sig2[ 7]);
		    addplogs(&dd,sig2[12]); addplogs(&dd,sig2[13]); addplogs(&dd,sig2[14]); addplogs(&dd,sig2[15]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
		    k++; d=-INFINITY; dd=-INFINITY;
		    addplogs(&d ,sig3[ 0]); addplogs(&d ,sig3[ 1]); addplogs(&d ,sig3[ 2]); addplogs(&d ,sig3[ 3]);
		    addplogs(&d ,sig3[ 4]); addplogs(&d ,sig3[ 5]); addplogs(&d ,sig3[ 6]); addplogs(&d ,sig3[ 7]);
		    addplogs(&dd,sig3[ 8]); addplogs(&dd,sig3[ 9]); addplogs(&dd,sig3[10]); addplogs(&dd,sig3[11]);
		    addplogs(&dd,sig3[12]); addplogs(&dd,sig3[13]); addplogs(&dd,sig3[14]); addplogs(&dd,sig3[15]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
			}
		}
	else
		for(i=0;i<Ly;i++)
			{
			for(k=0;k<16;k++) {sig0[k]=-INFINITY; sig1[k]=-INFINITY; sig2[k]=-INFINITY; sig3[k]=-INFINITY;}
			d0.r=yR[i];
			d0.i=yI[i];
			for(k=0;k<16;k++)
		        {
				CNsub(&tv.c[k],&d0,&d1); d=vv*CNabs(&d1);
		        addplogs(&sig0[k],d+getProb0(llrin,i,k));
		        addplogs(&sig1[k],d+getProb1(llrin,i,k));
		        addplogs(&sig2[k],d+getProb2(llrin,i,k));
		        addplogs(&sig3[k],d+getProb3(llrin,i,k));
		        }
		    k=4*i; d=-INFINITY; dd=-INFINITY;
			addplogs(&d ,sig0[ 0]); addplogs(&d ,sig0[ 2]); addplogs(&d ,sig0[ 4]); addplogs(&d ,sig0[ 6]);
		    addplogs(&d ,sig0[ 8]); addplogs(&d ,sig0[10]); addplogs(&d ,sig0[12]); addplogs(&d ,sig0[14]);
		    addplogs(&dd,sig0[ 1]); addplogs(&dd,sig0[ 3]); addplogs(&dd,sig0[ 5]); addplogs(&dd,sig0[ 7]);
		    addplogs(&dd,sig0[ 9]); addplogs(&dd,sig0[11]); addplogs(&dd,sig0[13]); addplogs(&dd,sig0[15]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
		    k++; d=-INFINITY; dd=-INFINITY;
		    addplogs(&d ,sig1[ 0]); addplogs(&d ,sig1[ 1]); addplogs(&d ,sig1[ 4]); addplogs(&d ,sig1[ 5]);
		    addplogs(&d ,sig1[ 8]); addplogs(&d ,sig1[ 9]); addplogs(&d ,sig1[12]); addplogs(&d ,sig1[13]);
		    addplogs(&dd,sig1[ 2]); addplogs(&dd,sig1[ 3]); addplogs(&dd,sig1[ 6]); addplogs(&dd,sig1[ 7]);
		    addplogs(&dd,sig1[10]); addplogs(&dd,sig1[11]); addplogs(&dd,sig1[14]); addplogs(&dd,sig1[15]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
		    k++; d=-INFINITY; dd=-INFINITY;
		    addplogs(&d ,sig2[ 0]); addplogs(&d ,sig2[ 1]); addplogs(&d ,sig2[ 2]); addplogs(&d ,sig2[ 3]);
		    addplogs(&d ,sig2[ 8]); addplogs(&d ,sig2[ 9]); addplogs(&d ,sig2[10]); addplogs(&d ,sig2[11]);
		    addplogs(&dd,sig2[ 4]); addplogs(&dd,sig2[ 5]); addplogs(&dd,sig2[ 6]); addplogs(&dd,sig2[ 7]);
		    addplogs(&dd,sig2[12]); addplogs(&dd,sig2[13]); addplogs(&dd,sig2[14]); addplogs(&dd,sig2[15]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
		    k++; d=-INFINITY; dd=-INFINITY;
		    addplogs(&d ,sig3[ 0]); addplogs(&d ,sig3[ 1]); addplogs(&d ,sig3[ 2]); addplogs(&d ,sig3[ 3]);
		    addplogs(&d ,sig3[ 4]); addplogs(&d ,sig3[ 5]); addplogs(&d ,sig3[ 6]); addplogs(&d ,sig3[ 7]);
		    addplogs(&dd,sig3[ 8]); addplogs(&dd,sig3[ 9]); addplogs(&dd,sig3[10]); addplogs(&dd,sig3[11]);
		    addplogs(&dd,sig3[12]); addplogs(&dd,sig3[13]); addplogs(&dd,sig3[14]); addplogs(&dd,sig3[15]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
			}

	free(alpha);
	free(beta);
	free(sig0);
  	free(sig1);
  	free(sig2);
  	free(sig3);
  	if (nrhs!=inputs) free(llrin);
	CVfree(&tv);
	}