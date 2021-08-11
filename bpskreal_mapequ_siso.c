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
	double *y,*ch,*tv,*llrout,*llrin,*alpha,*beta,nvar,vv,d,sig0,sig1;
	double sym[2]={+1,-1};

	if (!((nrhs==inputs)||(nrhs==inputs-1))||(nlhs!=outputs))
		mexErrMsgTxt("usage: llr_out=bpsk_mapequ_siso(y,nvar,h) or llr_out=bpsk_mapequ_siso(y,nvar,h,llr_in)\n");

	Lc=mxGetM(ch_if);
	Ly=mxGetM(y_if);
	nvar=mxGetScalar(nvar_in);
	vv=-0.5/nvar;
	M=Lc-1;
	Ls=1<<M;
	stmax=Ls-1;

	y=mxGetPr(y_if);
	ch=mxGetPr(ch_if);
	if (nrhs==inputs)
		{
		llrin=(double*) mxGetPr(llrin_if);
		for(i=0;i<Ly;i++) if (llrin[i]>LMAX) llrin[i]=LMAX; else if (llrin[i]<(-LMAX)) llrin[i]=-LMAX;
		}
	else
		{
		llrin=(double*)calloc(Ly,sizeof(double));
		for(i=0;i<Ly;i++) llrin[i]=0;
		}
	llrout_if=mxCreateDoubleMatrix(Ly,1,mxREAL);
	llrout=(double*)mxGetPr(llrout_if);

	tv=(double*)calloc(Ls<<1,sizeof(double));
	alpha=(double*)calloc(Ls*Ly,sizeof(double));
	beta=(double*)calloc(Ls*Ly,sizeof(double));
	for(i=0;i<(Ly*Ls);i++) {alpha[i]=-INFINITY; beta[i]=-INFINITY;}
	for(i=0;i<Ls;i++) {beta[i*Ly+Ly-1]=0.0; alpha[i*Ly]=0.0;} // termination
	for(i=0;i<(Ls<<1);i++) {tv[i]=0.0; for(l=0;l<=M;l++) tv[i]+=ch[l]*sym[(i>>l)&1];}

	if (Lc>1)
		{
		// obtain alphas
		for(i=0;i<Ly-1;i++)
			{
			for(j=0;j<Ls;j++)
				{
				st2=j<<1;
				st2lim=st2&stmax;
				addplogs(&alpha[i+1+st2lim*Ly    ],alpha[i+j*Ly]+vv*SQR(tv[st2  ]-y[i])+0.5*llrin[i]);
				addplogs(&alpha[i+1+(st2lim+1)*Ly],alpha[i+j*Ly]+vv*SQR(tv[st2+1]-y[i])-0.5*llrin[i]);
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
				st2lim=st2&stmax;
				addplogs(&beta[i+j*Ly],beta[i+1+st2lim*Ly    ]+vv*SQR(tv[st2  ]-y[i+1])+0.5*llrin[i+1]);
				addplogs(&beta[i+j*Ly],beta[i+1+(st2lim+1)*Ly]+vv*SQR(tv[st2+1]-y[i+1])-0.5*llrin[i+1]);
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
				st2lim=st2&stmax;
				addplogs(&sig0,vv*SQR(tv[st2  ]-y[i])+alpha[i+j*Ly]+beta[i+st2lim*Ly    ]);
				addplogs(&sig1,vv*SQR(tv[st2+1]-y[i])+alpha[i+j*Ly]+beta[i+(st2lim+1)*Ly]);
				}
			sig0-=sig1;	if (sig0<(-LMAX)) llrout[i]=-LMAX; else if (sig0>LMAX) llrout[i]=LMAX; else llrout[i]=sig0;
			}
		}
	else
		for(i=0;i<Ly;i++)
			{
			sig0=vv*SQR(tv[0]-y[i])-vv*SQR(tv[1]-y[i]);
			if (sig0<(-LMAX)) llrout[i]=-LMAX; else if (sig0>LMAX) llrout[i]=LMAX; else llrout[i]=sig0;
			}

	free(alpha);
	free(beta);
	if (nrhs!=inputs) free(llrin);
	CVfree(&tv);
	}
