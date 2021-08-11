clear all; clc;
%% Simulation
%
% We start by defining some common simulation parameters
nTx = 1;                  % Number of transmit antennas
nRx = 4;                  % Number of receive antennas
EbNoVec = -5:1:2; %-5:1:3; %-5:1:2;        % Eb/No in dB
modOrd = 2;             % constellation size = 2^modOrd
warning('off');

% Create a local random stream to be used by random number generators for repeatability.
stream = RandStream('mt19937ar');
if(modOrd == 1)
    % Create PSK modulator and demodulator System objects
    pskModulator   = comm.PSKModulator(...
                'ModulationOrder',  2^modOrd, ...
                'PhaseOffset',      0, ...
                'BitInput',         true);
    pskDemodulator = comm.PSKDemodulator( ...
                'ModulationOrder',  2^modOrd, ...
                'PhaseOffset',      0, ...
                'BitOutput',        true, ...
                'DecisionMethod', 'Log-likelihood ratio');
else%(modOrd == 2)
    % Create PSK modulator and demodulator System objects
    pskModulator   = comm.PSKModulator(...
                'ModulationOrder',  2^modOrd, ...
                'PhaseOffset',      pi/4, ...
               'BitInput',         true);
    pskDemodulator = comm.PSKDemodulator( ...
                'ModulationOrder',  2^modOrd, ...
                'PhaseOffset',      pi/4, ...
                'BitOutput',        true, ...
                'DecisionMethod', 'Log-likelihood ratio');
end
 
% Set a number of channel estimation iterations
iterChanEstim = 2;
% Create error rate calculation System objects for 3 different receivers
for i=1:iterChanEstim
    dec_CBERCalc{i} = comm.ErrorRate;
    eq_CBERCalc{i} = comm.ErrorRate;
end


% Pre-allocate variables to store BER results for speed
[BER_DEC, BER_EQ] = deal(zeros(3, length(EbNoVec), iterChanEstim));
[FER_DEC, FER_EQ] = deal(zeros(length(EbNoVec), iterChanEstim));


% Set up a figure for visualizing BER results
fig = figure;
grid on;
hold on;
ax = fig.CurrentAxes;
ax.YScale = 'log';
xlim([EbNoVec(1)-0.01 EbNoVec(end)]);
ylim([1e-3 1]);
xlabel('Eb/No (dB)');
ylabel('BER');
fig.NumberTitle = 'off';
fig.Renderer = 'zbuffer';
fig.Name = 'Spatial Multiplexing';
tit = nTx+"x"+nRx+" Coded QPSK System";
title(tit);
set(fig,'DefaultLegendAutoUpdate','off');

%Code settings
Ncode = 32; ncode = log2(Ncode);
Kcode = 16;
Rate = Kcode/Ncode;

% Code setup
indexMatrix = index_Matrix(Ncode);
Q = load('RelSeq_Q.mat');
Q = Q.Q;
Q1 = Q(Q<=Ncode);       %reliability sequence for N
F = Q1(1:Ncode-Kcode);     %frozen positions
M = Q1(Ncode-Kcode+1:end); %message positions

G2 = [1 0; 1 1]; Gn = G2;
for blk=1:ncode-1
    Gn = mod(kron(G2,Gn),2);
end

int = 1:Ncode; iint = 1:Ncode; iint(int) = iint; %1:Ncode (randperm(Ncode))

%% Pilots set up
pilotsNum = 4;
Nfull = Ncode/modOrd + pilotsNum;
if(mod(Nfull/pilotsNum,1) ~= 0)
    disp('Choose another number of pilots');
end
pilotsPositions = 1:Nfull/pilotsNum:Nfull;
dataPositions = setdiff(1:Nfull,pilotsPositions);

PP = [];
fprintf('Message length is %d\n', Kcode);


%% Loop over selected EbNo points
times = zeros(1,length(EbNoVec));
for idx = 1:length(EbNoVec)
    % Reset error rate calculation System objects
    for i=1:iterChanEstim
        reset(dec_CBERCalc{i});
        reset(eq_CBERCalc{i});
    end
    
    % Calculate SNR from EbNo for each independent transmission link
    snrIndB = EbNoVec(idx) + 10*log10(modOrd);
    snrLinear = 10^(0.1*snrIndB);
    sigma2 = 1/snrLinear;
    time1 = tic;
    maxWords = 5000; words = 0;
    while (words < maxWords)
        % Create random bit vector to modulate
        msg = randi(stream, [0 1], [Kcode, nTx]);
        % Modulate data
        %txSig = pskModulator(msg);
       
        for user = 1:nTx
            %create message
            u(:,user) = zeros(Ncode,1); 
            u(M,user) = msg(:,user);
            cwd(:,user) = mod(u(:,user)'*Gn,2)';
            cwdInter(:,user) = cwd(int,user);
            cwdMod(:,user) = pskModulator(cwdInter(:,user));
            
            %create pilots
            pilots(:,user) = randi(stream, [0 1], pilotsNum*modOrd, 1); %(randn(stream,pilotsNum*modOrd,1)>0); %randi(stream, [0 1], pilotsNum, 1);
            pilotsMod(:,user) = pskModulator(pilots(:,user));
            
            txSig(pilotsPositions,user) = pilotsMod(:,user);
            txSig(dataPositions,user) = cwdMod(:,user);
        end
        sMod = transpose(txSig);
        % Flat Rayleigh fading channel with independent links
        Hw = (randn(stream, nRx, nTx) +  1i*randn(stream, nRx, nTx))/sqrt(2);

        rxSig = awgn(Hw*sMod, snrIndB, 0, stream);
        
        
        % ZF-SIC receiver
        for iterCE = 1:iterChanEstim
            if(iterCE == 1)
                r = rxSig(:, dataPositions); 
                xx=sMod(:,pilotsPositions);
                yy=rxSig(:,pilotsPositions);

                H = yy/xx; 

            elseif(iterCE > 1)
                addPilots = PP;
                if(modOrd == 2)
                    addPilots = round(addPilots./2);
                end

                r = rxSig(:, dataPositions); 
                xx1 = sMod(:,pilotsPositions);
                yy1 = rxSig(:,pilotsPositions);
                
                
                xx2 = xmm(:,addPilots); 
                yy2 = r(:,addPilots);

                xx = [xx1 xx2];
                yy = [yy1 yy2];

                H = yy/xx; %Hw;
                PP = [];
            end
            
            % Initialization
            orderVec = 1:nTx;
            k = nTx+1;
            % Start ZF nulling loop
            for n = 1:nTx
                H = H(:, [1:k-1,k+1:end]);
                orderVec = orderVec(1, [1:k-1,k+1:end]);
                G = (H'*H + ((nTx-n+1)/snrLinear)*eye(nTx-n+1)) \ eye(nTx-n+1); % Same as inv(H'*H), but faster
                [~, k] = min(diag(G));
                symNum = orderVec(k);
                
                estimate = G(k,:) * H' * r;
                
                demodSoft = pskDemodulator(transpose(estimate));
                demodHard = double(demodSoft < 0);
   
                %deinterleaving
                msgEq = mod(demodHard(iint)'*Gn,2); msgEq = msgEq(M)';
                demodSoft = demodSoft(iint);
                
                %[msgDec, decSoft] = z_scan_decoder(F, M, demodSoft, 5); %msgData = msgData'; softData = softData';
                [msgDec, decSoft] = BP_Decoder_LLR_2(F, M, demodSoft, 10, indexMatrix,[],[],[]);
                %nL = 8; Rounds = 1;
                %[msgDec, decSoft] = scl_bp(F, M, nL, Rounds, indexMatrix, Gn, demodSoft,[]); msgDec = msgDec'; decSoft = decSoft';
                %softData = softData*0.1;
                
                decHard = double(decSoft < 0);
                xm = modulate(decSoft, modOrd);
                %Decoder feedback
                xmm(symNum,:) = pskModulator(decHard(int)); %Hard
                %xmm(symNum,:) = xm; %Soft
                [~,PP] = maxk(abs(decSoft(int)), 1); 
                %% Subtract the effect of the last decoded symbol from r
                if n < nTx
                    r = r - H(:, k) * pskModulator(decHard).';
                end
                % Update BER
                
                FER_DEC(idx, iterCE) = FER_DEC(idx, iterCE) + double(sum(msg(:,symNum) ~= msgDec)>0);
                FER_EQ(idx, iterCE) = FER_EQ(idx, iterCE) + double(sum(msg(:,symNum) ~= msgEq)>0);

                %BER_DEC(1, idx, iterCE) = BER_DEC(1, idx, iterCE)  + double(isequal(cwd(:,symNum), decHard)); %estZF
                BER_DEC(:, idx, iterCE) = dec_CBERCalc{iterCE}(cwd(:,symNum), decHard); %estZF
                BER_EQ(:, idx, iterCE) = eq_CBERCalc{iterCE}(cwd(:,symNum), demodHard); %estZF
            end
        end
        words = words + nTx;
    end
    times(idx) = toc(time1);
    FER_DEC(idx, :) = FER_DEC(idx, :) / words;
    FER_EQ(idx, :) = FER_EQ(idx, :) / words;
    % Plot results
    %semilogy(EbNoVec(1:idx), squeeze(FER_DEC(  1, 1:idx, :)), 'b*', EbNoVec(1:idx), squeeze(BER_EQ(  1, 1:idx, :)), 'r*');
    semilogy(EbNoVec(1:idx), squeeze(BER_DEC(1, 1:idx)), 'b*', EbNoVec(1:idx), squeeze(BER_EQ(1, 1:idx)), 'r*');
    %legend('ZF-SIC');
    
    %semilogy(EbNoVec(1:idx), squeeze(FER_DEC(1:idx, :)), 'b*');
    %semilogy(EbNoVec(1:idx), squeeze(FER_DEC(1:idx, :)), 'b*', EbNoVec(1:idx), squeeze(FER_EQ(1:idx, :)), 'r*');
    
    drawnow;
end

% Draw the lines
semilogy(EbNoVec, squeeze(FER_DEC(1:idx, :)));
hold off;

function xm = modulate(soft, modOrd)
    Ncode = length(soft);
    llr =  tanh(soft);
    %llr = soft;
    if(modOrd == 1)
        xm = llr;
    elseif(modOrd == 2)
        xm = (llr(2:2:Ncode)+llr(1:2:Ncode)*1i)/sqrt(2);
    end
end
