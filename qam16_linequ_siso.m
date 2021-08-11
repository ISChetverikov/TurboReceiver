function Lout=qam16_linequ_siso(y,nvar,h,N,Lin)

load('signal_constellations');
[NM,H,s,latency]=ut_equinit(h,N);
L=length(y)-N;
Lin=[Lin;zeros(3*L-length(Lin),1)];
P0=1-1./(1+exp(Lin)); l=1;
xm=zeros(L,1); xv=zeros(L,1);
for n=1:L
    P00=P0(l); P01=P0(l+1); P02=P0(l+2); P03=P0(l+3); P10=1-P00; P11=1-P01; P12=1-P02; P13=1-P03;
    w=[P00*P01*P02*P03 P10*P01*P02*P03 P00*P11*P02*P03 P10*P11*P02*P03 ...
        P00*P01*P12*P03 P10*P01*P12*P03 P00*P11*P12*P03 P10*P11*P12*P03 ...
        P00*P01*P02*P13 P10*P01*P02*P13 P00*P11*P02*P13 P10*P11*P02*P13 ...
        P00*P01*P12*P13 P10*P01*P12*P13 P00*P11*P12*P13 P10*P11*P12*P13];
    xm(n)=sym_16qam*w'; xv(n)=abs(sym_16qam-xm(n)).^2*w'; l=l+4;
end
xm=[xm;zeros(N,1)];
xv=[xv;ones(N,1)];
y=y-filter(h,1,xm);
Lout=zeros(4*L,1); p=zeros(16,1); l=1;
for n=1+latency:L+latency
    if (n<NM) x=flipud([zeros(NM-n,1);xv(1:n)]); else x=flipud(xv(n-NM+1:n)); end
    x(1+latency)=1;
    c=(H*diag(x)*H'+nvar*eye(N))\s;
    cf=flipud(c);
    if (n<N) z=cf(N-n+1:N)'*y(1:n); else z=cf'*y(n-N+1:n); end
    k=real(c'*s);
    z=z+k*xm(n-latency);
    for m=1:16 p(m)=exp(-abs(z-k*sym_16qam(m))^2/(k-k^2))+1e-10; end
    P00=P0(l); P01=P0(l+1); P02=P0(l+2); P03=P0(l+3); P10=1-P00; P11=1-P01; P12=1-P02; P13=1-P03;
    Lout(l  )=log(p(1)*P01*P02*P03+p(3)*P11*P02*P03+p(5)*P01*P12*P03+p(7)*P11*P12*P03+p(9)*P01*P02*P13+p(11)*P11*P02*P13+p(13)*P01*P12*P13+p(15)*P11*P12*P13)...
        -log(p(2)*P01*P02*P03+p(4)*P11*P02*P03+p(6)*P01*P12*P03+p(8)*P11*P12*P03+p(10)*P01*P02*P13+p(12)*P11*P02*P13+p(14)*P01*P12*P13+p(16)*P11*P12*P13);
    Lout(l+1)=log(p(1)*P00*P02*P03+p(2)*P10*P02*P03+p(5)*P00*P12*P03+p(6)*P10*P12*P03+p(9)*P00*P02*P13+p(10)*P10*P02*P13+p(13)*P00*P12*P13+p(14)*P10*P12*P13)...
        -log(p(3)*P00*P02*P03+p(4)*P10*P02*P03+p(7)*P00*P12*P03+p(8)*P10*P12*P03+p(11)*P00*P02*P13+p(12)*P10*P02*P13+p(15)*P00*P12*P13+p(16)*P10*P12*P13);
    Lout(l+2)=log(p(1)*P00*P01*P03+p(2)*P10*P01*P03+p(3)*P00*P11*P03+p(4)*P10*P11*P03+p(9)*P00*P01*P13+p(10)*P10*P01*P13+p(11)*P00*P11*P13+p(12)*P10*P11*P13)...
        -log(p(5)*P00*P01*P03+p(6)*P10*P01*P03+p(7)*P00*P11*P03+p(8)*P10*P11*P03+p(13)*P00*P01*P13+p(14)*P10*P01*P13+p(15)*P00*P11*P13+p(16)*P10*P11*P13);
    Lout(l+3)=log(p(1)*P00*P01*P02+p(2)*P10*P01*P02+p(3)*P00*P11*P02+p(4)*P10*P11*P02+p(5)*P00*P01*P12+p(6)*P10*P01*P12+p(7)*P00*P11*P12+p(8)*P10*P11*P12)...
        -log(p(9)*P00*P01*P02+p(10)*P10*P01*P02+p(11)*P00*P11*P02+p(12)*P10*P11*P02+p(13)*P00*P01*P12+p(14)*P10*P01*P12+p(15)*P00*P11*P12+p(16)*P10*P11*P12);
    l=l+4;
end