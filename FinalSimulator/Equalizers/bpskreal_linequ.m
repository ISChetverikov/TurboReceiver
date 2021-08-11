function Lout=bpskreal_linequ(y,nvar,h,N)

[NM,H,s,latency]=ut_equinit(h,N);
L=length(y)-N;
c=(H*H'+nvar*eye(N))\s;
z=filter(c,1,y);
z=z(1+latency:L+latency);
k=c'*s;
Lout=2*z/(1-k);