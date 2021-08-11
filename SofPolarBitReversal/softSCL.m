function [L_inf,R_prime] = softSCL(y,I,L,crc_poly,extr)
N=length(y);
n=log2(N);
M = index_Matrix(N);
fr(I) = 1/2;
[u,~,~,~] = polar_decode_list_extr(y,fr,L,0,[],extr);
crcdetector = comm.CRCDetector(crc_poly);

cands=[];
for l=1:L
    uu=u(l,1:end);
    [~,err]=crcdetector(uu(logical(fr))');
    if err == 0
        cands=[cands; uu];
    end
end
if isempty(cands)
    decoded_bits=u(1,:);
else
    decoded_bits=cands(1,:);
end


[L,v] = polar_decode_L2R(y, decoded_bits);

L=L';
LL=L(:,1);
L_inf=LL(fr>0)';

R=zeros(N,n+1);
R(:,1)=1e6*(~fr);

if ~isempty(cands)

    for j = 1:n
        for i = 1:N/2
            R(M(i, j),j+1) = f(R(2*i,j) + L(M(i,j)+N/2^j,j+1), R(2*i-1,j));
            R(M(i,j)+N/2^j,j+1) = f(R(2*i-1, j), L(M(i,j), j+1)) + R(2*i,j);
        end
    end

    R_prime=R(:,end);

    Inv = sign(R(:,end))==2*v(1,:)'-1;
    R_prime(Inv)=-R_prime(Inv);

    R_prime=R_prime';
else
    R_prime=y;
end
end

function result = f(x, y)
result = sign(x) * sign(y) * min(abs(x), abs(y));
end

function result = g(x, y, u)
if u == 0
    result = x + y;
elseif u == 1
    result = y - x;
else
    error('wrong u value.');
end
end

