clc; 
clear all
disp('Turbo equalization simulation ------------------------------------');
%rng('default');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compile mex files if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Checking whether or not MEX/DLL-file (3 returned) exists 
%then compiling the MEX-function from C source code
%if (exist('conv_encoder')~=3) mex('conv_encoder.c'); end
%if (exist('conv_decoder')~=3) mex('conv_decoder.c','utilities.c'); end
if (exist('bpsk_mapequ_siso')~=3) mex('bpsk_mapequ_siso.c','utilities.c'); end
if (exist('qpsk_mapequ_siso')~=3) mex('qpsk_mapequ_siso.c','utilities.c'); end
if (exist('psk8_mapequ_siso')~=3) mex('psk8_mapequ_siso.c','utilities.c'); end
if (exist('qam16_mapequ_siso')~=3) mex('qam16_mapequ_siso.c','utilities.c'); end
if (exist('bpskreal_mapequ_siso')~=3) mex('bpskreal_mapequ_siso.c','utilities.c'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define channel
channels={'ch_none' 'ch_A' 'ch_B' 'ch_C' 'ch_Anna'};
CHA=channels{3}; %3

% define equalizer
equalizers={'eq_map_det' 'eq_exact_lin' 'eq_approx_lin'};
EQU=equalizers{3};

% define bit-to-symbol mapping
mappings={'mo_bpsk_real' 'mo_bpsk' 'mo_qpsk' 'mo_8psk' 'mo_16qam'};
MOD=mappings{3};

% define block length (coded bits) and interleaver
blocklengths={'bl_noint' 'bl_128' 'bl_1024' 'bl_8192' 'bl_65536'};
BLO=blocklengths{2};
BLO='bl_128';

% define convolutional code
% codes={'co_none' 'co_r2_m2' 'co_r2_m4' 'co_r2_m6'};
% COD=codes{2};

% define puncturing pattern
patterns={'pu_1' 'pu_100' 'pu_1000000'};
PUN=patterns{1};

N=128;
K=72;
L = 8;
poly_length = 8;
POLY = 'z^4 + z^3 + z^2 + z + 1'; 
POLY = 'z^8 + z^7 + z^6 + z^4 + z^2 + 1';
% POLY = 'z^24 + z^23 + z^14 + z^12 + z^8 + 1';

%dec_func0 = @(y,I,L,crc_poly)softSCL(y,I,L,crc_poly,0);
%dec_func1 = @(y,I,L,crc_poly)softSCL(y,I,L,crc_poly,1);
%dec_func2 = @(y,I,L,crc_poly)softSCL(y,I,L,crc_poly,2);

%[4:0.5:6]
%[6:0.5:7.5] %scan
EbN0 = [4:0.5:7];     % SNR list in dB
Dtf=4;              % SNR down sampling for transfer function calculation
IT=5;               % number of iterations
NEQ=5;             % equalizer length
REC=1;              % recursive systematic encoding
MULT = 1;
MINERR=MULT * 100;          % minimum number of errors to quit simulation
MINSIM=MULT * 1e4;         % minimum number of code bits per simulation
MAXSIM=MULT * 5e3;         % maximum number of code bits per simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=')_(';  
settings_str_short =['(' BLO s PUN ')'];
settings_str = ['(' CHA s EQU s MOD s BLO s PUN ')'];

disp(['Initiate calculation for the setting ' settings_str ':']); 
t0=cputime;

% results file names
fd=['Results/result_tfd_' settings_str_short '.mat']; %transfer function decode
fe=['Results/result_tfe_' settings_str '.mat']; %transfer function equalization
fs=['Results/result_sim_' settings_str '.mat'];

% decoder transfer function
%{
if (exist(fd)>0)&&(1)
    disp(['Decoder transfer function for the setting ' se1 ' already available.']);
    load(fd);
else
    disp(['Start calculation of decoder transfer function for the setting ' se1 ':']);
    [TID,TOD]=decoder_transfer(BLO,COD,PUN,REC,MINSIM); 
    %return transfer function input TID, output TOD
    save(fd,'TID','TOD');
end
%}

% equalizer transfer function
%{
EbN0tf=EbN0(1:Dtf:end); 
TOEall=[];
for ne=1:length(EbN0tf)
    ss=num2str(EbN0tf(ne)); 
    ss(find(ss=='.'))='_'; 
    ss=['(en_' ss ')_'];
    fes=[fe(1:11) ss fe(12:end)];
    if (exist(fes)>0)&&(1)
        disp(['Equalizer transfer function for the setting ' ss settings_str ' already available.']);
        load(fes);
    else
        disp(['Start calculation of equalizer transfer function for the setting ' ss se2 ':']);
        [TIE,TOE]=equalizer_transfer(EbN0tf(ne),CHA,EQU,MOD,BLO,COD,PUN,NEQ,MINSIM);
        save(fes,'TIE','TOE');
    end
    TOEall=[TOEall;TOE];
end
%}

% simulation
FER0all=[];
FER1all=[];
%BERall=[]; 
%MIDall=[]; 
%MIEall=[];
for ne=1:length(EbN0)
    ss=num2str(EbN0(ne)); 
    ss(find(ss=='.'))='_'; 
    ss=['(en_' ss ')_'];
    fss=[fs(1:11) ss fs(12:end)];
    if (exist(fss)>0)&&(0)
        disp(['Simulation data for the setting ' ss settings_str ' already available.']);
        load(fss);
    else
        disp(['Start simulation for the setting ' ss settings_str ':']);
        [FER0,MIE0,MID0, B0]=simulation(EbN0(ne),CHA,EQU,MOD,BLO,PUN,NEQ,REC,IT,MINERR,MINSIM,MAXSIM,N,K,L,POLY,poly_length,0);
        % [FER1,MIE1,MID1, B1]=simulation(EbN0(ne),CHA,EQU,MOD,BLO,PUN,NEQ,REC,IT,MINERR,MINSIM,MAXSIM,N,K,L,POLY,poly_length,2);
        %save(fss,'BER','MID','MIE');
    end
    FER0all=[FER0all; FER0];
    %FER1all=[FER1all; FER1];
    %B_diff = xor(B0, B1);
    %diff = sum(B_diff, 'all');
    %disp(['Diff in generating: ' diff]);
    %BER0all=[BER0all;BER]; 
    %MIDall=[MIDall;MID]; 
    %MIEall=[MIEall;MIE];
end

t1=cputime;
disp(['Calculation finished in ' datestr((t1-t0)/3600/24,'dd HH:MM:SS') '.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FERall = [FER0all(:,1) FER0all(:,end)]; %FER1all(:,1) FER1all(:,end)];

figure();
h=semilogy(EbN0, FERall, '.-','linewidth', 2); 
set(gca, 'fontsize', 12); 
grid on;

% for n=0:IT
%     lstr{n+1}=['iteration ' num2str(n)];
% end

lstr = {'BP-0_0', ['BP-0_' num2str(IT)], 'BP-1_0', ['BP-1_' num2str(IT)]}; 

legend(h,lstr,'fontsize',10,'location','southwest');
xlabel('Eb/N0 in dB');
ylabel('Data frame error rate');
h1=text(EbN0(1),3.0,'Frame error rate performance of Turbo equalization for the settings');
h2=text(EbN0(1),1.7, settings_str); 
axis([EbN0(1) EbN0(end)+1e-3 1e-3 1.2]);
set([h1 h2],'fontsize',12,'interpreter','none','fontweight','bold');
set(gcf,'position',[0 0 800 600]); 
print(gcf,'-djpeg',['result_sim_' settings_str]);

%{
figure();
for ne=1:length(EbN0tf)
    me=[]; 
    md=0; 
    nes=(ne-1)*Dtf+1;
    for l=1:IT+1 
        me=[me MIEall(nes,l) MIEall(nes,l)]; 
    end
    for l=1:IT+1 
        md=[md MIDall(nes,l) MIDall(nes,l)]; 
    end
    me=[me MIEall(nes,IT+1)];
    plot(TIE,TOEall(ne,:),'b-','linewidth',2);
    hold on
    plot(TOD,TID,'r-','linewidth',2);
    plot(md,me,'ko-','linewidth',2);
    hold off
    axis([0 1 0 1]); 
    set(gca,'fontsize',12); 
    grid on; 
    daspect([1 1 1]);
    h=legend('Equalizer','Decoder','Simulation','location','southeast'); 
    set(h,'fontsize',10);
    xlabel('Mutual information at the decoder output or equalizer input');
    ylabel('Mutual information at the decoder input or equalizer output');
    ss=num2str(EbN0tf(ne)); 
    ss(find(ss=='.'))='_';
    h1=text(0,1.058,'Decoder/Equalizer transfer function for the settings');
    h2=text(0,1.023,['(en_' ss ')_' settings_str]);
    set([h1 h2],'fontsize',12,'interpreter','none','fontweight','bold');
    set(gcf,'position',[800 0 700 700]);
    print(gcf,'-djpeg',['result_tfun_(en_' ss ')_' settings_str]);
end
%}




