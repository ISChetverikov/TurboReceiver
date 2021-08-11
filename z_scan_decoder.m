function [msg_cap, cwd_soft] = z_scan_decoder(F, M, llr, ScanIter)

N = length(llr);
K = length(M);
n = log2(N);   %power of Kroneker product
Rate = K / N;

f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a),abs(b)); %minsum
g = @(a,b,c) b+(1-2*c).*a; %g function
        
%SC decoder
L = zeros(n+1,N);  %beliefs
ucap = zeros(n+1,N); %decisions
ns = zeros(1,2*N-1); %node state vector

beta = zeros(n,N); %decisions
ucap(n+1,F) = Inf;
L(1,:) = llr; %belief of root

for it = 1:ScanIter
    node = 0; depth = 0; %start at root
    done = 0; %decoder is finished or not
    ns = zeros(1,2*N-1);
    while(done == 0) %traverse till all bits are decoded
        %leaf or not
        if(depth == n)
            ucap(n+1,F) = Inf;
            ucap(n+1,M) = 0;
            if(node == (N-1))
                done = 1;
            else
                node = floor(node/2);
                depth = depth - 1;
                %changing the state of the leaf is not so important
            end
        else
            %nonleaf
            npos = (2^depth-1) + node + 1; %position of node in node state vector
            if(ns(npos) == 0) %step L and go to left child
                %disp('L');
                %disp([node depth]);
                temp = 2^(n-depth);
                Ln = L(depth+1,temp*node+1:temp*(node+1)); %incoming beliefs
                a = Ln(1:temp/2); b = Ln(temp/2+1:end); %split beliefs into 2

                rnode = 2 * node + 1; ldepth = depth + 1; %right child
                rtemp = temp/2; %2^(n-ldepth);
                ucapn = ucap(ldepth+1,rtemp*rnode+1:rtemp*(rnode+1)); %incoming decisions from right child

                node = node * 2; depth = depth + 1; %next node: left child
                temp = temp / 2; %incoming belief length for left child
                L(depth+1,temp*node+1:temp*(node+1)) = f(a,b+ucapn); %minsum and storage
                ucap(depth+1,temp*node+1:temp*(node+1)) = L(depth+1,temp*node+1:temp*(node+1));
                ns(npos) = 1; %change state
            else
                if(ns(npos) == 1) %step R and go to right child
                    %disp('R');
                    %disp([node depth]);
                    temp = 2^(n-depth);
                    Ln = L(depth+1,temp*node+1:temp*(node+1)); %incoming beliefs
                    a = Ln(1:temp/2); b = Ln(temp/2+1:end); %split beliefs into 2

                    lnode = 2 * node; ldepth = depth + 1; %left child
                    ltemp = temp/2; %2^(n-ldepth);
                    ucapn = ucap(ldepth+1,ltemp*lnode+1:ltemp*(lnode+1)); %incoming decisions from left child

                    node = node * 2 + 1; depth = depth + 1; %next node: right child
                    temp = temp/2; %incoming belief length for left child
                    L(depth+1,temp*node+1:temp*(node+1)) = f(a,ucapn) + b; %g and storage
                    ucap(depth+1,temp*node+1:temp*(node+1)) = L(depth+1,temp*node+1:temp*(node+1));
                    ns(npos) = 2; %change state
                else %ns(npos) == 2: step U and go to parent
                    temp = 2^(n-depth);

                    Ln = L(depth+1,temp*node+1:temp*(node+1)); %incoming beliefs
                    a = Ln(1:temp/2); b = Ln(temp/2+1:end); %split beliefs into 2

                    lnode = 2 * node; rnode = 2 * node + 1; cdepth = depth + 1; %left and rigth child
                    ctemp = temp/2; %2^(n-ldepth);
                    ucapl = ucap(cdepth+1,ctemp*lnode+1:ctemp*(lnode+1)); %incoming decisions from left child
                    ucapr = ucap(cdepth+1,ctemp*rnode+1:ctemp*(rnode+1)); %incoming decisions from right child

                    ucap(depth+1,temp*node+1:temp*(node+1)) = [f(ucapl, b+ucapr) ucapr+f(a,ucapl)];
                    %ucap(depth+1,temp*node+1:temp*(node+1)) = [mod(ucapl+ucapr,2) ucapr]; %combine
                    node = floor(node/2); depth = depth - 1;
                end
            end
        end
    end
    depth = 1; node = 0;
    temp = 2^(n-depth+1);
    Ln = L(depth,temp*node+1:temp*(node+1)); %incoming beliefs
    a = Ln(1:temp/2); b = Ln(temp/2+1:end); %split beliefs into 2
    
    lnode = 2 * node; rnode = 2 * node + 1;
    ctemp = temp/2; %2^(n-ldepth);
    ucapl = ucap(depth+1,ctemp*lnode+1:ctemp*(lnode+1)); %incoming decisions from left child
    ucapr = ucap(depth+1,ctemp*rnode+1:ctemp*(rnode+1)); %incoming decisions from right child
    
    ucap(depth,temp*node+1:temp*(node+1)) = [f(ucapl, b+ucapr) ucapr+f(a,ucapl)];
    %ucap(depth+1,temp*node+1:temp*(node+1)) = [mod(ucapl+ucapr,2) ucapr]; %combine
    node = floor(node/2); depth = depth - 1;
end
cwd_soft = L(1,:) + ucap(1,:);
usoft = L(n+1,:) + ucap(n+1,:);
uhard = usoft < 0;
msg_cap = double(uhard(M));
end

%hold on
%semilogy(EbNodB, FER_sim);
%figure();
%semilogy(EbNodB, FER_sim, EbNodB, BER_sim);
%grid on;
%legend('FER', 'BER');
%xlabel('E_b/N_0, dB'); ylabel('FER/BER');
