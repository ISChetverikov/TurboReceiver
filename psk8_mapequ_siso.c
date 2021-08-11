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

double getProb (double* pd, int i, int symi)
	{
	double d;
	int i0=3*i, i1=3*i+1, i2=3*i+2;
	switch (symi)
		{
		case 0: d= pd[i0]+pd[i1]+pd[i2]; break;
		case 1: d=-pd[i0]+pd[i1]+pd[i2]; break;
		case 2: d= pd[i0]-pd[i1]+pd[i2]; break;
		case 3: d=-pd[i0]-pd[i1]+pd[i2]; break;
		case 4: d= pd[i0]+pd[i1]-pd[i2]; break;
		case 5: d=-pd[i0]+pd[i1]-pd[i2]; break;
		case 6: d= pd[i0]-pd[i1]-pd[i2]; break;
		case 7: d=-pd[i0]-pd[i1]-pd[i2]; break;
		}
	return 0.5*d;
	}

double getProb2 (double* pd, int i, int symi)
	{
	double d;
	int i0=3*i, i1=3*i+1;
	switch (symi)
		{
		case 0: case 4: d= pd[i0]+pd[i1]; break;
		case 1: case 5: d=-pd[i0]+pd[i1]; break;
		case 2: case 6: d= pd[i0]-pd[i1]; break;
		case 3: case 7: d=-pd[i0]-pd[i1]; break;
		}
	return 0.5*d;
	}

double getProb1 (double* pd, int i, int symi)
	{
	double d;
	int i0=3*i, i2=3*i+2;
	switch (symi)
		{
		case 0: case 2: d= pd[i0]+pd[i2]; break;
		case 1: case 3: d=-pd[i0]+pd[i2]; break;
		case 4: case 6: d= pd[i0]-pd[i2]; break;
		case 5: case 7: d=-pd[i0]-pd[i2]; break;
		}
	return 0.5*d;
	}

double getProb0 (double* pd, int i, int symi)
	{
	double d;
	int i1=3*i+1, i2=3*i+2;
	switch (symi)
		{
		case 0: case 1: d=+pd[i1]+pd[i2]; break;
		case 2: case 3: d=-pd[i1]+pd[i2]; break;
		case 4: case 5: d=+pd[i1]-pd[i2]; break;
		case 6: case 7: d=-pd[i1]-pd[i2]; break;
		}
	return 0.5*d;
	}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{
	int Lc,Ly,M,Ls,i,j,k,l,st2lim,stmax,st2;
	double *yR,*yI,*chR,*chI,*llrout,*llrin,*alpha,*beta,nvar,vv,d,dd,*sig0,*sig1,*sig2;
	double symR[8]={-SQRTH,0,0,SQRTH,-1,-SQRTH,SQRTH,1};
  	double symI[8]={SQRTH,-1,1,-SQRTH,0,-SQRTH,SQRTH,0};
	CV tv;
	CN d0,d1;

	if (!((nrhs==inputs)||(nrhs==inputs-1))||(nlhs!=outputs))
		mexErrMsgTxt("usage: llr_out=psk8_mapequ_siso(y,nvar,h) or llr_out=psk8_mapequ_siso(y,nvar,h,llr_in)\n");

	Lc=mxGetM(ch_if);
	Ly=mxGetM(y_if);
	nvar=mxGetScalar(nvar_in);
	vv=-1.0/nvar;
	M=Lc-1;
	Ls=1<<(3*M);
	stmax=Ls-1;

	yR=mxGetPr(y_if);
	yI=mxGetPi(y_if);
	chR=mxGetPr(ch_if);
	chI=mxGetPi(ch_if);
	if (nrhs==inputs)
		{
		llrin=(double*)mxGetPr(llrin_if);
		for(i=0;i<(Ly*3);i++) if (llrin[i]>LMAX) llrin[i]=LMAX; else if (llrin[i]<(-LMAX)) llrin[i]=-LMAX;
		}
	else
		{
		llrin=(double*)calloc(Ly*3,sizeof(double));
		for(i=0;i<(Ly*3);i++) llrin[i]=0;
		}
	llrout_if=mxCreateDoubleMatrix(Ly*3,1,mxREAL);
	llrout=(double*)mxGetPr(llrout_if);

	sig0=(double*)calloc(8,sizeof(double));
  	sig1=(double*)calloc(8,sizeof(double));
  	sig2=(double*)calloc(8,sizeof(double));
	alpha=(double*)calloc(Ls*Ly,sizeof(double));
	beta=(double*)calloc(Ls*Ly,sizeof(double));
	for(i=0;i<(Ly*Ls);i++) {alpha[i]=-INFINITY; beta[i]=-INFINITY;}
	for(i=0;i<Ls;i++) {beta[i*Ly+Ly-1]=0.0; alpha[i*Ly]=0.0;} // termination

 	CValloc(&tv,Ls<<3);
	for(i=0;i<(Ls<<3);i++)
		{
		d=0; dd=0;
		if (chI!=0) for(l=0;l<=M;l++)
			{
			d +=chR[l]*symR[(i>>(l*3))&7]-chI[l]*symI[(i>>(l*3))&7];
			dd+=chR[l]*symI[(i>>(l*3))&7]+chI[l]*symR[(i>>(l*3))&7];
			}
		else for(l=0;l<=M;l++)
			{
			d +=chR[l]*symR[(i>>(l*3))&7];
			dd+=chR[l]*symI[(i>>(l*3))&7];
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
				st2=j<<3;
				d0.r=yR[i];
				d0.i=yI[i];
				st2lim=st2&stmax;
				for(k=0;k<8;k++)
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
				st2=j<<3;
				d0.r=yR[i+1];
				d0.i=yI[i+1];
				st2lim=st2&stmax;
				for(k=0;k<8;k++)
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
			for(k=0;k<8;k++) {sig0[k]=-INFINITY; sig1[k]=-INFINITY; sig2[k]=-INFINITY;}
			for(j=0;j<Ls;j++)
				{
				st2=j<<3;
				d0.r=yR[i];
				d0.i=yI[i];
				st2lim=st2&stmax;
				for(k=0;k<8;k++)
			        {
					CNsub(&tv.c[st2+k],&d0,&d1); d=vv*CNabs(&d1);
			        d+=alpha[i+j*Ly]+beta[i+(st2lim+k)*Ly];
			        addplogs(&sig0[k],d+getProb0(llrin,i,k));
			        addplogs(&sig1[k],d+getProb1(llrin,i,k));
			        addplogs(&sig2[k],d+getProb2(llrin,i,k));
			        }
				}
			k=3*i; d=-INFINITY; dd=-INFINITY;
		    addplogs(&d ,sig0[0]); addplogs(&d ,sig0[2]); addplogs(&d ,sig0[4]); addplogs(&d ,sig0[6]);
		    addplogs(&dd,sig0[1]); addplogs(&dd,sig0[3]); addplogs(&dd,sig0[5]); addplogs(&dd,sig0[7]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
		    k++; d=-INFINITY; dd=-INFINITY;
		    addplogs(&d ,sig1[0]); addplogs(&d ,sig1[1]); addplogs(&d ,sig1[4]); addplogs(&d ,sig1[5]);
		    addplogs(&dd,sig1[2]); addplogs(&dd,sig1[3]); addplogs(&dd,sig1[6]); addplogs(&dd,sig1[7]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
		    k++; d=-INFINITY; dd=-INFINITY;
		    addplogs(&d ,sig2[0]); addplogs(&d ,sig2[1]); addplogs(&d ,sig2[2]); addplogs(&d ,sig2[3]);
		    addplogs(&dd,sig2[4]); addplogs(&dd,sig2[5]); addplogs(&dd,sig2[6]); addplogs(&dd,sig2[7]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
			}
		}
	else
		for(i=0;i<Ly;i++)
			{
			for(k=0;k<8;k++) {sig0[k]=-INFINITY; sig1[k]=-INFINITY; sig2[k]=-INFINITY;}
			d0.r=yR[i];
			d0.i=yI[i];
			for(k=0;k<8;k++)
		        {
				CNsub(&tv.c[k],&d0,&d1); d=vv*CNabs(&d1);
		        addplogs(&sig0[k],d+getProb0(llrin,i,k));
		        addplogs(&sig1[k],d+getProb1(llrin,i,k));
		        addplogs(&sig2[k],d+getProb2(llrin,i,k));
		        }
			k=3*i; d=-INFINITY; dd=-INFINITY;
		    addplogs(&d ,sig0[0]); addplogs(&d ,sig0[2]); addplogs(&d ,sig0[4]); addplogs(&d ,sig0[6]);
		    addplogs(&dd,sig0[1]); addplogs(&dd,sig0[3]); addplogs(&dd,sig0[5]); addplogs(&dd,sig0[7]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
		    k++; d=-INFINITY; dd=-INFINITY;
		    addplogs(&d ,sig1[0]); addplogs(&d ,sig1[1]); addplogs(&d ,sig1[4]); addplogs(&d ,sig1[5]);
		    addplogs(&dd,sig1[2]); addplogs(&dd,sig1[3]); addplogs(&dd,sig1[6]); addplogs(&dd,sig1[7]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
		    k++; d=-INFINITY; dd=-INFINITY;
		    addplogs(&d ,sig2[0]); addplogs(&d ,sig2[1]); addplogs(&d ,sig2[2]); addplogs(&d ,sig2[3]);
		    addplogs(&dd,sig2[4]); addplogs(&dd,sig2[5]); addplogs(&dd,sig2[6]); addplogs(&dd,sig2[7]);
		    d-=dd; if (d<(-LMAX)) llrout[k]=-LMAX; else if (d>LMAX) llrout[k]=LMAX; else llrout[k]=d;
			}

	free(alpha);
	free(beta);
	free(sig0);
  	free(sig1);
  	free(sig2);
  	if (nrhs!=inputs) free(llrin);
	CVfree(&tv);
	}