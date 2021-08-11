function [Tin,Tout]=decoder_transfer(BLO,COD,PUN,REC,SIM)
% input:
% BLO:   'block interleaving type' {'bl_noint' 'bl_128' 'bl_1024' 'bl_8192' 'bl_65536'}
% COD:   'convolutional code type' {'co_none' 'co_r2_m2' 'co_r2_m4' 'co_r2_m6'}
% PUN:   'punturing type' {'pu_1' 'pu_100' 'pu_1000000'}
% REC=1:           recursive systematic encoding
% SIM=MINSIM=1e5:  minimum number of code bits per simulation

% return:
% Tin:   transfer function input
% Tout:  transfer function output

% creating system's parameters for simulation
[ch,sym,mem,gen,Nc,int,prate,punc]=initialize('ch_none','mo_bpsk',BLO,COD,PUN);
% ch:    channel model; 
% sym:   signal constellation
% mem:   code constrain length
% gen:   hex-formed generating function
% Nc:    interleaving size
% int:   interleaving pattern
% prate: 1/(punturing rate)
% punc:  punturing pattern

% load LLR variance data set
load('llr_variances');
% acc=1.e-6    
% ia=[0.005 .025:.05:.975 .995]     
% steps=10001   
% ans='/home/micha/matlab/dist'     
% llrvar=[lookup table]  
Tin=ia;                 % transfer function input
Nt=length(Tin);         % number of transfer function points
rate=prate/length(gen); % code rate
Nb=Nc*rate-mem;         % block length (data bits)
Tout=zeros(1,Nt);       % transfer function output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transfer function calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nt=1:Nt % number of transfer function points
    Lall=[]; % list of decoder output LLRs
    
    % initialize random number generator
    rand('state',sum(93251*clock)); 
    randn('state',sum(59206*clock));
    
    for ns=1:ceil(SIM/Nc) % # code bits / interleaving size = # frames
        % data bit generation
        b=double(rand(Nb,1)>0.5); % Nb block length (data bits)
        
        % encoding and interleaving
        c=conv_encoder(b,gen,REC); 
        c=c(punc);
        
        polar_encode();
        
        % generate input LLRs, built-in lookup table llrvar(nt)
        Lin=(randn(Nc,1)*sqrt(llrvar(nt))+llrvar(nt)/2).*(1-2*c);
        %Nc interleaving size
        
        % decoder
        Linu=zeros(Nc*prate,1); 
        Linu(punc)=Lin;
        [bh Ld]=conv_decoder(Linu,gen,REC); 
        Ld=Ld(punc); 
        Lall=[Lall;Ld.*(1-2*c)];
    end
    
    % mutual information
    Tout(nt)=ut_sample2mi(Lall);
end

disp(['Tout=' num2str(Tout,'%1.2g ')]);






