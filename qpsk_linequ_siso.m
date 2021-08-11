function Lout=qpsk_linequ_siso(y,nvar,h,N,Lin)

[NM,H,s,latency]=ut_equinit(h,N);
L=length(y)-N;
th=tanh(Lin*0.5);
xm=(th(1:2:2*L)+th(2:2:2*L)*i)/sqrt(2);
xm=[xm;zeros(N,1)];
y=y-filter(h,1,xm);
Lout=zeros(2*L,1); l=1;
for n=1+latency:L+latency
    if (n<NM) x=flipud([zeros(NM-n,1);xm(1:n)]); else x=flipud(xm(n-NM+1:n)); end
    x(1+latency)=0;
    c=(H*diag(1-x.*conj(x))*H'+nvar*eye(N))\s;
    cf=flipud(c);
    if (n<N) z=cf(N-n+1:N)'*y(1:n); else z=cf'*y(n-N+1:n); end
    k=real(c'*s);
    z=z+k*xm(n-latency);
    Lout(l  )=sqrt(8)*real(z)/(1-k);
    Lout(l+1)=sqrt(8)*imag(z)/(1-k);
    l=l+2;
end