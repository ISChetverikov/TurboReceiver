function Lout=psk8_linequ_avsiso(y,nvar,h,N,Lin)

load('signal_constellations');
[NM,H,s,latency]=ut_equinit(h,N);
L=length(y)-N;
Lin=[Lin;zeros(3*L-length(Lin),1)];
th=tanh(Lin*0.5);
xm=th(1:3:3*L)/2*(0.5*i-0.5+i/sqrt(2))+th(2:3:3*L)/2*(-i*0.5-0.5-1/sqrt(2))+...
    th(3:3:3*L)/2.*(th(1:3:3*L)/2*(1-sqrt(2)+i)+th(2:3:3*L)/2*(1-i+sqrt(2)*i));
av=mean(1-xm.*conj(xm));
xm=[xm;zeros(N,1)];
y=y-filter(h,1,xm);
c=(H*H'*av+(1-av)*s*s'+nvar*eye(N))\s;
k=real(c'*s);
p=real(H'*c);
q=real(nvar*c'*c);
z=filter(c,1,y);
z=z(1+latency:L+latency);
z=z+k*xm(1:L);
zv=zeros(L,1);
for n=1+latency:L+latency
    if (n<NM) x=flipud([zeros(NM-n,1);xm(1:n)]); else x=flipud(xm(n-NM+1:n)); end
    x(1+latency)=1; zv(n-latency)=p'*diag(1-x.*conj(x))*p+q;
end
P0=1-1./(1+exp(Lin))+1e-10; Lout=zeros(3*L,1); p=zeros(8,1); l=1;
for n=1:L
    for m=1:8 p(m)=exp(-abs(z(n)-k*sym_8psk(m))^2/zv(n))+1e-10; end
    P00=P0(l); P01=P0(l+1); P02=P0(l+2); P10=1-P00; P11=1-P01; P12=1-P02;
    Lout(l  )=log(p(1)*P01*P02+p(3)*P11*P02+p(5)*P01*P12+p(7)*P11*P12)-log(p(2)*P01*P02+p(4)*P11*P02+p(6)*P01*P12+p(8)*P11*P12);
    Lout(l+1)=log(p(1)*P00*P02+p(2)*P10*P02+p(5)*P00*P12+p(6)*P10*P12)-log(p(3)*P00*P02+p(4)*P10*P02+p(7)*P00*P12+p(8)*P10*P12);
    Lout(l+2)=log(p(1)*P00*P01+p(2)*P10*P01+p(3)*P00*P11+p(4)*P10*P11)-log(p(5)*P00*P01+p(6)*P10*P01+p(7)*P00*P11+p(8)*P10*P11);
    l=l+3;
end