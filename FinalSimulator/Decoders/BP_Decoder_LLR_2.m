function [u_llr, c_llr, codeword] = BP_Decoder_LLR_2(frozen_bits, info_bits, llr, max_iter)
N = length(llr);
n = log2(N);
R = zeros(N, n+1);
L = zeros(N, n+1);

M=index_Matrix(N);

R(frozen_bits, n+1) = realmax;
L(:, 1) = llr;

for iter = 1:max_iter

    for j = 1:n
        for i = 1:N/2
            L(M(i, j),j+1) = f(L(M(i,j),j), L(M(i,j)+N/2^j,j) + R(M(i,j)+N/2^j,j+1));
            L(M(i,j)+N/2^j,j+1) = f(R(M(i,j),j+1), L(M(i,j),j)) + L(M(i,j)+N/2^j,j);
        end
    end
    
    for j = n:-1:1
        for i = 1:N/2
            R(M(i,j),j) = f(R(M(i,j),j+1), L(M(i,j)+N/2^j,j) + R(M(i,j)+N/2^j,j+1));
            R(M(i,j)+N/2^j,j) = f(R(M(i,j),j+1), L(M(i,j), j)) + R(M(i,j)+N/2^j,j+1);
        end
    end
   
end

u_hat = L(:,1) < 0;

u_llr = L(:,n+1)+R(:,n+1); u_llr = double(u_llr(info_bits) < 0);
c_llr = R(:,1);
codeword = c_llr < 0;

info_esti = u_hat(info_bits);
end
