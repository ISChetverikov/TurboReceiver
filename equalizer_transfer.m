function [Tin,Tout]=equalizer_transfer(EbN0,CHA,EQU,MOD,BLO,COD,PUN,NEQ,SIM)

%[ch,sym,~,gen,Nc,~,prate,~]
[ch,sym,mem,gen,Nc,int,prate,punc]=initialize(CHA,MOD,BLO,COD,PUN);
% ch:    channel model; 
% sym:   signal constellation
% mem:   code constrain length
% gen:   hex-formed generating function
% Nc:    interleaving size
% int:   interleaving pattern
% prate: 1/(punturing rate)
% punc:  punturing pattern

load('llr_variances');  % load LLR variance data set
% acc=1.e-6    
% ia=[0.005 .025:.05:.975 .995]     
% steps=10001   
% ans='/home/micha/matlab/dist'     
% llrvar=[lookup table] 
Tin=ia;                 % transfer function input
Nt=length(Tin);         % number of transfer function points
rate=prate/length(gen); % code rate
Nq=log2(length(sym));   % bits per modulation symbol
Tout=zeros(1,Nt);       % transfer function output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transfer function calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set noise variance
nvar=ch'*ch/10^(EbN0/10)/Nq/rate;

for nt=1:Nt
    Lall=[]; % list of equalizer output LLRs
    
    % initialize random number generator
    rand('state',sum(30924*clock)); 
    randn('state',sum(61094*clock));
    
    for ns=1:ceil(SIM/Nc)
        % data bit generation
        c=double(rand(Nc,1)>0.5);
        
        % modulate
        Nx=ceil(Nc/Nq); 
        cr=reshape([c;zeros(Nx*Nq-Nc,1)],Nq,Nx);
        x=transpose(sym(2.^(0:Nq-1)*cr+1));
        
        % channel
        if (strcmp(EQU,'eq_map_det'))
            y=filter(ch,1,x)+randn(Nx,2)*[1;1i]*sqrt(nvar/2);
        else
            y=filter(ch,1,[x;x(1:NEQ)])+randn(Nx+NEQ,2)*[1;1i]*sqrt(nvar/2);
        end
        
        % generate input LLRs
        Lin=(randn(Nc,1)*sqrt(llrvar(nt))+llrvar(nt)/2).*(1-2*c);
        
        % equalizer
        switch EQU
            case 'eq_map_det',
                switch (MOD)
                    case 'mo_bpsk_real', Le=bpskreal_mapequ_siso(real(y),nvar/2,ch,Lin);
                    case 'mo_bpsk',      Le=bpsk_mapequ_siso(y,nvar,ch,Lin);
                    case 'mo_qpsk',      Le=qpsk_mapequ_siso(y,nvar,ch,Lin);
                    case 'mo_8psk',      Le=psk8_mapequ_siso(y,nvar,ch,Lin);
                    case 'mo_16qam',     Le=qam16_mapequ_siso(y,nvar,ch,Lin);
                end
            case 'eq_exact_lin',
                switch (MOD)
                    case 'mo_bpsk_real', Le=bpskreal_linequ_siso(real(y),nvar/2,ch,NEQ,Lin);
                    case 'mo_bpsk',      Le=bpsk_linequ_siso(y,nvar,ch,NEQ,Lin);
                    case 'mo_qpsk',      Le=qpsk_linequ_siso(y,nvar,ch,NEQ,Lin);
                    case 'mo_8psk',      Le=psk8_linequ_siso(y,nvar,ch,NEQ,Lin);
                    case 'mo_16qam',     Le=qam16_linequ_siso(y,nvar,ch,NEQ,Lin);
                end
            case 'eq_approx_lin',
                switch (MOD)
                    case 'mo_bpsk_real', Le=bpskreal_linequ_avsiso(real(y),nvar/2,ch,NEQ,Lin);
                    case 'mo_bpsk',      Le=bpsk_linequ_avsiso(y,nvar,ch,NEQ,Lin);
                    case 'mo_qpsk',      Le=qpsk_linequ_avsiso(y,nvar,ch,NEQ,Lin);
                    case 'mo_8psk',      Le=psk8_linequ_avsiso(y,nvar,ch,NEQ,Lin);
                    case 'mo_16qam',     Le=qam16_linequ_avsiso(y,nvar,ch,NEQ,Lin);
                end
        end
        Lall=[Lall;Le(1:Nc).*(1-2*c)];        
    end
    
    % mutual information
    Tout(nt)=ut_sample2mi(Lall);
end

disp(['EbN0=' num2str(EbN0) ' Tout=' num2str(Tout,'%1.2g ')]);

