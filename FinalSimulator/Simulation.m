function FER = Simulation(ebn0, N, K, modulation, H, iterations, MINERR, MAXSIM, decoder)
% input:
% EbN0=0:0.5:20;  % SNR list in dB
% return:
% FER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% pilots

pilotsNum = 2;
modOrd = 2;
Nfull = N / modOrd + pilotsNum;
if(mod(Nfull/pilotsNum, 1) ~= 0)
    disp('Choose another number of pilots');
end
pilotsPositions = 1:Nfull/pilotsNum:Nfull;
dataPositions = setdiff(1:Nfull,pilotsPositions);


FER=zeros(1, iterations);

switch (modulation)
    % Create PSK modulator and demodulator System objects
    case 'bpsk'
    
        modulator   = comm.PSKModulator(...
                    'ModulationOrder',  2^modOrd, ...
                    'PhaseOffset',      0, ...
                    'BitInput',         true);
        demodulator = comm.PSKDemodulator( ...
                    'ModulationOrder',  2^modOrd, ...
                    'PhaseOffset',      0, ...
                    'BitOutput',        true, ...
                    'DecisionMethod', 'Log-likelihood ratio');
    case 'qpsk'
    % Create PSK modulator and demodulator System objects
        modulator   = comm.PSKModulator(...
                    'ModulationOrder',  2^modOrd, ...
                    'PhaseOffset',      pi/4, ...
                   'BitInput',         true);
        demodulator = comm.PSKDemodulator( ...
                    'ModulationOrder',  2^modOrd, ...
                    'PhaseOffset',      pi/4, ...
                    'BitOutput',        true, ...
                    'DecisionMethod', 'Log-likelihood ratio');
    otherwise
        error("Unknown modulation");
end

%Nq=log2(length(sym));   % bits per modulation symbol

snrIndB = ebn0 + 10*log10(modOrd);
snrLinear = 10 ^ (0.1 * snrIndB);
sigma2 = 1 / snrLinear;

Q = load('Data/RelSeq_Q.mat');
Q = Q.Q;
Q1 = Q(Q<=N);       %reliability sequence for N
F = Q1(1:N-K);     %frozen positions
M = Q1(N - K + 1:end); %message positions

% initialize random number generator
%rng(512);

%code settings

rate = K / N;
n = log2(N);
G2 = [1 0; 1 1]; G = G2;
for i=1:n-1
    G = mod(kron(G2,G),2);
end

isCrcUsed = isKey(decoder, 'CrcPoly');

% CRC
polyLength = 0;
if (isCrcUsed)
    crc_generator = comm.CRCGenerator(decoder('CrcPoly'));
    polyLength = decoder('PolyLength');
end

% set noise variance
ch = H;
nvar=ch'*ch/10^(ebn0/10)/ modOrd / rate;

blocksCount = zeros(1,iterations); % block count

% simulate a block transmission
while (blocksCount(end) < MAXSIM && FER(end) < MINERR)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % transmitter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % data bit generation
    msg = double(rand(K - polyLength, 1) > 0.5);
    msg_ext = zeros(1, N);
    if (isCrcUsed)
        msg_crc = crc_generator(msg);
        msg_ext(M) = msg_crc;
    else
        msg_ext(M) = msg;
    end

    % encode
    codeword = polar_encode(msg_ext, N)'; 
    
    %create pilots
    x = modulator(codeword);
    p = randi([0 1], pilotsNum * modOrd, 1); %(randn(stream,pilotsNum*modOrd,1)>0); %randi(stream, [0 1], pilotsNum, 1);
    pMod = modulator(p);

    xFull = zeros(Nfull, 1);
    xFull(pilotsPositions) = pMod;
    xFull(dataPositions) = x;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xFull = xFull';
    %y = awgn(H * xFull, snrIndB);
    temp = H * xFull;
    y = temp + randn(size(temp)).* (1 + 1i) * sqrt(nvar/2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % receiver
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    yData = y(:, dataPositions); 
    xPilots = xFull(:, pilotsPositions);
    yPilots = y(:, pilotsPositions);
    
    
    for it=1:iterations % loop by iteration
        if(it == 1)
            Hest = yPilots / xPilots;
            Lin = MseEqualize(yData, Hest, nvar, []);
            
        else
            addPilots = PP;
            if(modOrd == 2)
              addPilots = round(addPilots./2);
            end

            xData = xFull(:, dataPositions);
            xAddPilots = xData(:, addPilots); 
            yAddPilots = yData(:, addPilots);

            xPilots = [xPilots xAddPilots];
            yPilots = [yPilots yAddPilots];

            PP = [];
            Hest = yPilots / xPilots;
            Lin = MseEqualize(yData, Hest, nvar, []);
        end  

        switch (decoder('Type'))
            case 'BP'
                [msgDec, Lout, codewordDec] = BP_Decoder_LLR_2(F, M, Lin, decoder('Iterations'));
                Lout = Lout';
            case 'softSCL'
                [~, Lout, msgDec, codewordDec] = softSCL(Lin', M, decoder('L'), decoder('CrcPoly'), decoder('BpType'));
            otherwise
                error('Unknown decoder type');
        end
        Lout_hard = 1 - 2 * codewordDec;
        % Lout = Lout_hard
        [~, PP] = maxk(abs(Lout), 1);
        
        xOut = modulateLlr(Lout, modOrd);
        
        if (FER(it) >= MINERR)
            continue
        end
        blocksCount(it) = blocksCount(it) + 1;
        if(mean(abs(msgDec - msg)) > 0)
            FER(it) = FER(it) + 1;
        end
    end
end
FER = FER ./ blocksCount;

% display results
disp(['EbN0=' num2str(ebn0) ' ITeff=' num2str(iterations) ' blocks=' num2str(blocksCount) ' FER=' num2str(FER,'%1.2g ')]);

end

function xOut = modulateLlr(Lout, modOrd)
    N = length(Lout);
    llr =  tanh(Lout);
    %llr = soft;
    if(modOrd == 1)
        xOut = llr;
    elseif(modOrd == 2)
        xOut = (llr(2:2:N)+llr(1:2:N)*1i)/sqrt(2);
    else
        error('Unknown modulation');
    end
end

function L = demodulateX(xOut, modOrd)
    N = length(xOut);
    L = zeros(2 * N,1);
    L(1:2:2*N) = sqrt(8)*real(xOut) / (1-k);
    L(2:2:2*N) = sqrt(8)*imag(xOut) / (1-k);
end
