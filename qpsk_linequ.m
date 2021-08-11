function Lout=qpsk_linequ(y,nvar,h,N)

[NM,H,s,latency]=ut_equinit(h,N);
L=length(y)-N;
c=(H*H'+nvar*eye(N))\s;
z=filter(c,1,y);
z=z(1+latency:L+latency);
k=real(c'*s);
Lout=zeros(2*L,1);
Lout(1:2:2*L)=sqrt(8)*real(z)/(1-k);
Lout(2:2:2*L)=sqrt(8)*imag(z)/(1-k);