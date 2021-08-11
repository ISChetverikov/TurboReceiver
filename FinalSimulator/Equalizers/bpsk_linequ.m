function Lout=bpsk_linequ(y,nvar,h,N)

[NM,H,s,latency]=ut_equinit(h,N);
y=[y;y(1:latency)];
c=(H*H'+nvar/2*eye(N))\s;
z=filter(c,1,y);
z=z(1+latency:end);
k=real(c'*s);
Lout=2*real(z)/(1-k);