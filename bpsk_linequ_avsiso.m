function Lout=bpsk_linequ_avsiso(y,nvar,h,N,Lin)

[NM,H,s,latency]=ut_equinit(h,N);
L=length(y)-N;
xm=tanh(Lin*0.5);
av=mean(1-xm.*xm);
xm=[xm;zeros(N,1)];
y=y-filter(h,1,xm);
c=(H*H'*av+nvar/2*eye(N))\s;
k=c'*s;
z=filter(c,1,y);
z=z(1+latency:L+latency);
z=z+k*xm(1:L);
Lout=2/(1-av*k)*real(z);
