function [FER,BER,MIE,MID,B]=simulation(EbN0,CHA,EQU,MOD,BLO,PUN,NEQ,REC,IT,MINERR,MINSIM,MAXSIM,N,K,L,POLY,poly_length,bp_type)
% input:
% EbN0=0:0.5:20;  % SNR list in dB
% CHA:   'channel model' {'ch_none' 'ch_A' 'ch_B' 'ch_C'}
% EQU:   'define equalizer' {'eq_map_det' 'eq_exact_lin' 'eq_approx_lin'};   
% MOD:   'modulation type' {'mo_bpsk_real' 'mo_bpsk' 'mo_qpsk' 'mo_8psk' 'mo_16qam'}
% BLO:   'block interleaving type' {'bl_noint' 'bl_128' 'bl_1024' 'bl_8192' 'bl_65536'}
% PUN:   'punturing type' {'pu_1' 'pu_100' 'pu_1000000'}
% NEQ=15;     % equalizer length
% REC=1:      % recursive systematic encoding
% IT=9;       % number of iterations      
% MINERR=50;  % minimum number of errors to quit simulation   
% MINSIM=1e5; % minimum number of code bits per simulation
% MAXSIM=1e6; % maximum number of code bits per simulation
% return:
% BER:
% MIE:
% MID: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ch,sym,Nc,int,prate,punc]=initialize(CHA,MOD,BLO,PUN);
%return
%ch:    channel model; 
%sym:   signal constellation
%Nc:    interleaving size
%int:   interleaving pattern
%prate: 1/(punturing rate)
%punc:  punturing pattern
if(BLO)
iint(int)=1:Nc;         % deinterleaver
%rate=prate/length(gen); % code rate
%Nb=Nc*rate;         % block length (data bits)
Nq=log2(length(sym));   % bits per modulation symbol
BER=zeros(1,IT);      % data bit error rate
FER=zeros(1,IT);
MIE=zeros(1,IT);      % mutual information of equalizer output LLRs
MID=zeros(1,IT);      % mutual information of decoder output LLRs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcnt=0; % bit error count
blocnt=0; % block count
Leall=[]; % list of equalizer output LLRs
Ldall=[]; % list of decoder output LLRs

% initialize random number generator
%rand('state',sum(70245*clock)); 
%randn('state',sum(93286*clock));


rng(512);
B = [];

%code settings
%Nc = N;
Nb = K;
rate = K / N;
punc = 1:N; 
%int = 1:Nc;  iint = 1:Nc; iint(int)=1:Nc;
n = log2(N);
G2 = [1 0; 1 1]; G = G2;
for i=1:n-1
    G = mod(kron(G2,G),2);
end

indexMatrix = index_Matrix(Nc);
Q = load('RelSeq_Q.mat');
Q = Q.Q;
Q1 = Q(Q<=N);       %reliability sequence for N
F = Q1(1:N-K);     %frozen positions
M = Q1(N-K+1:end); %message positions

% set noise variance
nvar=ch'*ch/10^(EbN0/10)/Nq/rate;

% CRC generator
crc_generator = comm.CRCGenerator(POLY);

errcnt = zeros(1,IT); % bit error count
blocnt = zeros(1,IT); % block count

% simulate a block transmission
while (blocnt(end) < MAXSIM && FER(end) < MINERR)
        
    % why zeros here
    Leall=[zeros(N,IT);Leall];
    Ldall=[zeros(Nc,IT);Ldall];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % transmitter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % data bit generation
    b = double(rand(Nb-poly_length,1)>0.5);
    b_crc = crc_generator(b);
    uu = zeros(1,Nc); uu(M) = b_crc;

    % DEBUG 

    %c=conv_encoder(b,gen,REC); c=c(punc); c=c(int);
    cc=polar_encode(uu,Nc)'; 
    c=cc(punc); 
    c=c(int);

    % modulate
    Nx=ceil(Nc/Nq);
    cr=reshape([c;zeros(Nx*Nq-Nc,1)],Nq,Nx);
    x=transpose(sym(2.^(0:Nq-1)*cr+1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (strcmp(EQU,'eq_map_det'))
        y=filter(ch,1,x)+randn(Nx,2)*[1;1i]*sqrt(nvar/2);
    else
        y=filter(ch,1,[x;x(1:NEQ)])+randn(Nx+NEQ,2)*[1;1i]*sqrt(nvar/2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % receiver
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % equalization and deinterleaving
    for it=1:IT % loop by iteration


        if (it==1)
            switch EQU
                case 'eq_map_det',
                    switch (MOD)
                        case 'mo_bpsk_real', Le=bpskreal_mapequ_siso(real(y),nvar/2,ch);
                        case 'mo_bpsk',      Le=bpsk_mapequ_siso(y,nvar,ch);
                        case 'mo_qpsk',      Le=qpsk_mapequ_siso(y,nvar,ch);
                        case 'mo_8psk',      Le=psk8_mapequ_siso(y,nvar,ch);
                        case 'mo_16qam',     Le=qam16_mapequ_siso(y,nvar,ch);
                    end
                case {'eq_exact_lin','eq_approx_lin'},
                    switch (MOD)
                        case 'mo_bpsk_real', Le=bpskreal_linequ(real(y),nvar/2,ch,NEQ);
                        case 'mo_bpsk',      Le=bpsk_linequ(y,nvar,ch,NEQ);
                        case 'mo_qpsk',      Le=qpsk_linequ(y,nvar,ch,NEQ);
                        case 'mo_8psk',      Le=psk8_linequ(y,nvar,ch,NEQ);
                        case 'mo_16qam',     Le=qam16_linequ(y,nvar,ch,NEQ);
                    end
            end
        else
            switch EQU
                case 'eq_map_det',
                    switch (MOD)
                        case 'mo_bpsk_real', Le=bpskreal_mapequ_siso(real(y),nvar/2,ch,Ld);
                        case 'mo_bpsk',      Le=bpsk_mapequ_siso(y,nvar,ch,Ld);
                        case 'mo_qpsk',      Le=qpsk_mapequ_siso(y,nvar,ch,Ld);
                        case 'mo_8psk',      Le=psk8_mapequ_siso(y,nvar,ch,Ld);
                        case 'mo_16qam',     Le=qam16_mapequ_siso(y,nvar,ch,Ld);
                    end
                case 'eq_exact_lin',
                    switch (MOD)
                        case 'mo_bpsk_real', Le=bpskreal_linequ_siso(real(y),nvar/2,ch,NEQ,Ld);
                        case 'mo_bpsk',      Le=bpsk_linequ_siso(y,nvar,ch,NEQ,Ld);
                        case 'mo_qpsk',      Le=qpsk_linequ_siso(y,nvar,ch,NEQ,Ld);
                        case 'mo_8psk',      Le=psk8_linequ_siso(y,nvar,ch,NEQ,Ld);
                        case 'mo_16qam',     Le=qam16_linequ_siso(y,nvar,ch,NEQ,Ld);
                    end
                case 'eq_approx_lin',
                    switch (MOD)
                        case 'mo_bpsk_real', Le=bpskreal_linequ_avsiso(real(y),nvar/2,ch,NEQ,Ld);
                        case 'mo_bpsk',      Le=bpsk_linequ_avsiso(y,nvar,ch,NEQ,Ld);
                        case 'mo_qpsk',      Le=qpsk_linequ_avsiso(y,nvar,ch,NEQ,Ld);
                        case 'mo_8psk',      Le=psk8_linequ_avsiso(y,nvar,ch,NEQ,Ld);
                        case 'mo_16qam',     Le=qam16_linequ_avsiso(y,nvar,ch,NEQ,Ld);
                    end
            end
        end
    Leall(1:Nc,it)=Le(1:Nc).*(1-2*c); %Le=Le(iint);
    Leu=zeros(Nc*prate,1); Leu(punc)=Le;

    % decoding and interleaving
    %[bh, Ld] = z_scan_decoder(F, M, Leu, 2);

    % DEBUG
    %Leu = [2995.41336494294;-3013.35534387100;-3044.97136616593;-2997.33560358008;2946.57282761782;2997.97709825946;2977.83928017960;2958.05823010378;2942.86691749557;-3107.50244143779;-2942.69771009553;-2965.05754714817;2925.56742136829;3107.54188636563;-3003.08819043366;-2923.81156604213;-2971.15697300063;2919.27050225658;-3028.75280661210;3010.10814628387;-3133.55318745752;-2893.84767014759;-2952.86760652498;2954.26455989276;2910.91927110516;-3045.21603406567;3034.01199506910;-3051.03916987710;2940.26909546812;3013.22405912085;-2957.76711880748;2924.73999739717];

    [~, Ldl, bhl, v] = softSCL(Leu', M, L, POLY, bp_type);
    %[bhl, Ldl, codewordDec] = BP_Decoder_LLR_2(F, M, Leu, 20);
    %Ldl = Ldl';
    %bhl = bhl < 0;
    bh = bhl;
    Ld = Ldl';
    % hard
    %Ld = 1 - 2 * v';
    
    
    %Ld = Ld';
    %bh = bh';
    Ld=Ld(punc);
    Ld=Ld(int);
    Ldall(1:Nc,it)=Ld.*(1-2*c);

    if (FER(it) >= MINERR)
        continue
    end
    blocnt(it) = blocnt(it) + 1;
    %B = [B bh];
    % bit error rate measurement
    %BER(it+1)=(BER(it+1)*(blocnt-1)+mean(abs(bh-b)))/blocnt;
    if(mean(abs(bh-b)) > 0)
        FER(it) = FER(it) + 1;
    end

    %errcnt = errcnt + sum(abs(bh-b));

    end
end
FER = FER ./ blocnt;
% check for unstable iteration results
%{
for n=2:IT+1
    if (BER(n)>BER(n-1)*1.05)
        for k=n:IT+1
            BER(k)=BER(n-1);
            Leall(:,k)=Leall(:,n-1);
            Ldall(:,k)=Ldall(:,n-1);
        end
        break
    end
end
%}
% mutual information calculation
for it=1:IT
    MIE(it:end)=ut_sample2mi(Leall(:,it));
    MID(it:end)=ut_sample2mi(Ldall(:,it));
end

% display results
disp(['EbN0=' num2str(EbN0) ' ITeff=' num2str(IT) ' blocks=' num2str(blocnt) ...
    ' FER=' num2str(FER,'%1.2g ') ...
    ' MIE=' num2str(MIE,'%1.2g ') ' MID=' num2str(MID,'%1.2g ')]);
%m=min(find(BER==0)); if (~isempty(m)) IT=min(IT,m); end
end
