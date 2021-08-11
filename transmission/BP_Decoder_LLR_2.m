
function [u_llr, c_llr] = BP_Decoder_LLR_2(frozen_bits, info_bits, llr, max_iter, M, msg,cwd,u)
N = length(llr);
n = log2(N);
R = zeros(N, n+1);
L = zeros(N, n+1);
Q = zeros(N, n+1);
%іхКј»ЇR

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

% u_llr = L(:,1)+R(:,1);
% c_llr = R(:,n+1);

u_llr = L(:,n+1)+R(:,n+1); u_llr = double(u_llr(info_bits) < 0);
c_llr = R(:,1);

info_esti = u_hat(info_bits);
end


% function [u_llr, c_llr] = BP_Decoder_LLR_2(info_bits, frozen_bits, llr, max_iter, M, msg,cwd,u)
% N = length(llr);
% n = log2(N);
% R = zeros(N, n+1);
% L = zeros(N, n+1);
% Q = zeros(N, n+1);
% %іхКј»ЇR
% for i = 1:N
%     if frozen_bits(i) == 1
%         R(i, n+1) = realmax;
%     end
% end
% 
% L(:, 1) = llr;
% 
% for iter = 1:max_iter
% 
%     for j = 1:n
%         for i = 1:N/2
%             L(M(i, j),j+1) = f(L(M(i,j),j), L(M(i,j)+N/2^j,j) + R(M(i,j)+N/2^j,j+1));
%             L(M(i,j)+N/2^j,j+1) = f(R(M(i,j),j+1), L(M(i,j),j)) + L(M(i,j)+N/2^j,j);
%         end
%     end
%     
%     for j = n:-1:1
%         for i = 1:N/2
%             R(M(i,j),j) = f(R(M(i,j),j+1), L(M(i,j)+N/2^j,j) + R(M(i,j)+N/2^j,j+1));
%             R(M(i,j)+N/2^j,j) = f(R(M(i,j),j+1), L(M(i,j), j)) + R(M(i,j)+N/2^j,j+1);
%         end
%     end
%     
% %     for j = 1:n
% %         for i = 1:N/2
% %             R(M(i, j),j+1) = f(R(2*i,j) + L(M(i,j)+N/2^j,j+1), R(2*i-1,j));
% %             R(M(i,j)+N/2^j,j+1) = f(R(2*i-1, j), L(M(i,j), j+1)) + R(2*i,j);
% %         end
% %     end
% %     
% %     for j = n:-1:1
% %         for i = 1:N/2
% %             L(2*i-1,j) = f(L(M(i,j),j+1), L(M(i,j)+N/2^j,j+1) + R(2*i,j));
% %             L(2*i,j) = L(M(i,j)+N/2^j,j+1) + f(L(M(i,j),j+1), R(2*i-1, j));
% %         end
% %     end
% 
%     
% end
% 
% 
% u_llr = L(:,n+1)+R(:,n+1);
% c_llr = R(:,1);
% 
% end