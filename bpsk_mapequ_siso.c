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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{
	int Lc,Ly,M,Ls,i,j,k,l,st2lim,stmax,st2;
	double *yR,*yI,*chR,*chI,*llrout,*llrin,*alpha,*beta,nvar,vv,d,sig0,sig1;
	double symR[2]={+1,-1};
	CV tv;
	CN d0,d1;

	if (!((nrhs==inputs)||(nrhs==inputs-1))||(nlhs!=outputs))
		mexErrMsgTxt("usage: llr_out=bpsk_mapequ_siso(y,nvar,h) or llr_out=bpsk_mapequ_siso(y,nvar,h,llr_in)\n");

	Lc=mxGetM(ch_if);
	Ly=mxGetM(y_if);
	nvar=mxGetScalar(nvar_in);
	vv=-1.0/nvar;
	M=Lc-1;
	Ls=1<<M;
	stmax=Ls-1;

	yR=mxGetPr(y_if);
	yI=mxGetPi(y_if);
	chR=mxGetPr(ch_if);
	chI=mxGetPi(ch_if);
	if (nrhs==inputs)
		{
		llrin=(double*)mxGetPr(llrin_if);
		for(i=0;i<Ly;i++) if (llrin[i]>LMAX) llrin[i]=LMAX; else if (llrin[i]<(-LMAX)) llrin[i]=-LMAX;
		}
	else
		{
		llrin=(double*)calloc(Ly,sizeof(double));
		for(i=0;i<Ly;i++) llrin[i]=0;
		}
	llrout_if=mxCreateDoubleMatrix(Ly,1,mxREAL);
	llrout=(double*)mxGetPr(llrout_if);

	CValloc(&tv,Ls<<1);
	alpha=(double*)calloc(Ls*Ly,sizeof(double));
	beta=(double*)calloc(Ls*Ly,sizeof(double));
	for(i=0;i<(Ly*Ls);i++) {alpha[i]=-INFINITY; beta[i]=-INFINITY;}
	for(i=0;i<Ls;i++) {beta[i*Ly+Ly-1]=0.0; alpha[i*Ly]=0.0;} // termination

	for(i=0;i<(Ls<<1);i++)
		{
		tv.c[i].r=0.0; tv.c[i].i=0.0;
		if (chI!=0) for(l=0;l<=M;l++) {tv.c[i].r+=chR[l]*symR[(i>>l)&1]; tv.c[i].i+=chI[l]*symR[(i>>l)&1];}
		else for(l=0;l<=M;l++) tv.c[i].r+=chR[l]*symR[(i>>l)&1];
		}

	if (Lc>1)
		{
		// obtain alphas
		for(i=0;i<Ly-1;i++)
			{
			for(j=0;j<Ls;j++)
				{
				st2=j<<1;
				d0.r=yR[i];
				d0.i=yI[i];
				st2lim=st2&stmax;
				CNsub(&tv.c[st2],&d0,&d1);
				addplogs(&alpha[i+1+st2lim*Ly],alpha[i+j*Ly]+vv*CNabs(&d1)+0.5*llrin[i]);
				CNsub(&tv.c[st2+1],&d0,&d1);
				addplogs(&alpha[i+1+(st2lim+1)*Ly],alpha[i+j*Ly]+vv*CNabs(&d1)-0.5*llrin[i]);
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
				st2=j<<1;
				d0.r=yR[i+1];
				d0.i=yI[i+1];
				st2lim=st2&stmax;
				CNsub(&tv.c[st2],&d0,&d1);
				addplogs(&beta[i+j*Ly],beta[i+1+st2lim*Ly]+vv*CNabs(&d1)+0.5*llrin[i+1]);
				CNsub(&tv.c[st2+1],&d0,&d1);
				addplogs(&beta[i+j*Ly],beta[i+1+(st2lim+1)*Ly]+vv*CNabs(&d1)-0.5*llrin[i+1]);
				}
			d=-INFINITY;
			for(l=0;l<Ls;l++) if (beta[i+l*Ly]>d) d=beta[i+l*Ly];
			for(l=0;l<Ls;l++) beta[i+l*Ly]-=d;
			}

		// obtain output
		for(i=0;i<Ly;i++)
			{
			sig0=-INFINITY; sig1=-INFINITY;
			for(j=0;j<Ls;j++)
				{
				st2=j<<1;
				d0.r=yR[i];
				d0.i=yI[i];
				st2lim=st2&stmax;
				CNsub(&tv.c[st2],&d0,&d1);
				addplogs(&sig0,vv*CNabs(&d1)+alpha[i+j*Ly]+beta[i+st2lim*Ly]);
				CNsub(&tv.c[st2+1],&d0,&d1);
				addplogs(&sig1,vv*CNabs(&d1)+alpha[i+j*Ly]+beta[i+(st2lim+1)*Ly]);
				}
			sig0-=sig1;
			if (sig0<(-LMAX)) llrout[i]=-LMAX; else if (sig0>LMAX) llrout[i]=LMAX; else llrout[i]=sig0;
			}
		}
	else
		for(i=0;i<Ly;i++)
			{
			d0.r=yR[i];
			d0.i=yI[i];
			CNsub(&tv.c[0],&d0,&d1);
			sig0=vv*CNabs(&d1);
			CNsub(&tv.c[1],&d0,&d1);
			sig0-=vv*CNabs(&d1);
			if (sig0<(-LMAX)) llrout[i]=-LMAX; else if (sig0>LMAX) llrout[i]=LMAX; else llrout[i]=sig0;
			}

	free(alpha);
	free(beta);
	if (nrhs!=inputs) free(llrin);
	CVfree(&tv);
	}
