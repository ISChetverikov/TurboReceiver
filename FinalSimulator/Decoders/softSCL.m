function [L_inf, R_prime, msg, v] = softSCL(y,I,L,crc_poly, bp_type)
N=length(y);
% n=log2(N);
% M = index_Matrix(N);
% M = flip(M, 2);
fr = zeros(1,N);
fr(I) = 1;
[u,~,~,~] = polar_decode_list(y,fr,L,0,[]);
crcdetector = comm.CRCDetector(crc_poly);

cands=[];
for l=1:L
    uu=u(l,1:end);
    words = uu(I)';
    [~,err]=crcdetector(words);
    if err == 0
        cands=[cands; uu];
    end
end
if isempty(cands)
    decoded_bits=u(1,:);
else
    decoded_bits=cands(1,:);
end
[msg, ~] = crcdetector(decoded_bits(I)');

[L, v] = polar_decode_L2R(y, decoded_bits);

L=L';
LL=L(:,1);
L_inf=LL(fr>0)';

% R=zeros(N,n+1);
% R(:,1)=1e6 * (~fr);
% R(fr>0,1) = LL(fr>0);

switch (bp_type)
    case 0
        bp_func = @bp_0;
    case 1
        bp_func = @bp_1;
    case 2
        bp_func = @bp_2;
    case 9
        bp_func = @bp_9;
    otherwise
        error('Unknown belief propagation for the soft SCL decoder');
end

if ~isempty(cands)
    R_prime = bp_func(L, v, fr);
else
    R_prime=y;
end
end

% belief propagation from the article about soft SCL
function R_prime = bp_0(L, v, fr)
    N=size(L, 1);
    n=log2(N);
    M = index_Matrix(N);
    M = flip(M, 2);
    
    R=zeros(N, n + 1);
    R(:, 1) = 1e6 * (~fr);
    
    for j = 1:n
        for i = 1:N/2
            m = M(i, j);
            shift = 2 ^ (j - 1);
            R_L1 = R(m + shift, j) + L(m + shift, j + 1);
            R(m, j + 1) = f(R(m, j), R_L1);
            
            f2 = f(R(m, j), L(m, j + 1));
            R(m + shift, j + 1) = f2 + R(m + shift, j);
            
            % R(M(i, j), j + 1) = f(R(M(i, j), j), R(M(i, j) + 2 ^ (j - 1), j)) + L( M(i, j) + 2 ^ (j - 1), j + 1);
            % R(M(i, j) + 2^(j - 1), j+1) = f( R(M(i, j), j), L(M(i, j), j + 1)) + R(M(i,j) + 2^(j - 1), j);
            % R(M(i, j), j+1) = f(R(2*i,j) + L(M(i,j)+N/2^j,j+1), R(2*i-1,j));
            % R(M(i,j)+N/2^j,j+1) = f(R(2*i-1, j), L(M(i,j), j+1)) + R(2*i,j);
        end
    end

    R_prime=R(:,end);

    Inv = sign(R(:,end))==2*v(1,:)'-1;
    R_prime(Inv)=-R_prime(Inv);

    R_prime=R_prime';
end

% belief propagation: use codeword with flip
function R_prime = bp_1(L, v, fr)
    N=size(L, 1);
    n=log2(N);
    M = index_Matrix(N);
    M = flip(M, 2);
    
    R = L(:, end);

    Inv = sign(R)==2*v(1,:)'-1;
    R(Inv)=-R(Inv);

    R_prime=R';
end

% belief propagation: use codeword with RELU
function R_prime = bp_2(L, v, fr)
    N=size(L, 1);
    n=log2(N);
    
    R = L(:, end);

    Inv = sign(R)==2*v(1,:)'-1;
    R(Inv)=0.0;

    R_prime=R';
end

% belief propagation: random llrs
function R_prime = bp_9(L, v, fr)
    N=size(L, 1);

    R_prime=randn(1, N);
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

