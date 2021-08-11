function [u, x, idx, PM_out] = polar_decode_list(y, f, L, PM, u)
% y = bit APP from channel in output order (#row = index in list)
% f = input a priori probs in input order
% L = maximum list size
% listprob = probabilities of each path in the list
% x = output hard decision in output order
% u = input hard decisions in input order
% idx = indices of the list items used
% Recurse down to length 1
if nargin < 5, u = []; end
if nargin < 4, PM_out = 0; end
if nargin < 3, L = 1; end
if ~iscell(f), f = num2cell(f); end
N = size(y,2);
if (N==1)
    step_PM=[(y<0).*y, -(y>=0).*y];  % metrics of all possible ways to continue paths,
                                % zeros - first column, ones - second column 
    L0 = size(y,1);             % number of input paths
    if (f{1} == 1)
        % Make list of all possible decisions (double the list size)
        % ones are first, zeros - second
        [PM_out, j] = sort([PM + step_PM(:,2); PM + step_PM(:,1)], 1, 'descend');
        % Trim list if it is too large
        if length(j) > L
            PM_out = PM_out(1:L);
            j = j(1:L);
        end
        % the first half of the list was due to x=1
        x = j <= L0;
        idx = mod(j-1,L0)+1;
        u = [u(idx,:) x];
        % determine which path each decision belongs to
    else
        x = repmat(f{1},size(y));
        u = [u x];
        idx = 1:L0;
        PM_out = PM + step_PM(:,f{1}+1); % Or better avoid modifing it?
    end
else
    % Compute soft mapping back one stage
    y1 = cnop(y(:,1:N/2), y(:,N/2+1:end));
    [u1, x1, idx1, PM1] = polar_decode_list(y1, f(1:(N/2)), L, PM, u);
%     if (size(y, 1) ~= size(x1, 1))
%         disp('True');
%     end
    y2 = vnop(y(idx1, 1:N/2), y(idx1, N/2+1:end), x1);
    [u2, x2, idx2, PM2] = polar_decode_list(y2, f((N/2+1):end), L, PM1, u1);
    
    % Tunnel u decisions back up. Compute and interleave x1,x2 hard decisions
    u = u2;
    x = zeros(length(idx2), size(x2, 2) * 2);
    x(:, 1:N/2) = bitxor(x1(idx2,:), x2);
    x(:, N/2+1:end) = x2;
    PM_out = PM2;
    idx = idx1(idx2);
end
end

function l = cnop(l1,l2)
l = (sign(l1).*sign(l2)).* min(abs(l1),abs(l2));
end

function l = vnop(l1, l2, bit)
l = (1 - 2*bit).* l1 + l2;
end

