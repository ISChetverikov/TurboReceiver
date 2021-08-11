function [ch,sym,Nc,int,prate,punc]=initialize(CHA,MOD,BLO,PUN)
%input
%CHA:   'channel model' {'ch_none' 'ch_A' 'ch_B' 'ch_C'}
%MOD:   'modulation type' {'mo_bpsk_real' 'mo_bpsk' 'mo_qpsk' 'mo_8psk' 'mo_16qam'}
%BLO:   'block interleaving type' {'bl_noint' 'bl_128' 'bl_1024' 'bl_8192' 'bl_65536'}
%PUN:   'punturing type' {'pu_1' 'pu_100' 'pu_1000000'}

%return
%ch:    channel model; 
%sym:   signal constellation
%Nc:    interleaving size
%int:   interleaving pattern
%prate: 1/(punturing rate)
%punc:  punturing pattern

% select channel and number of iterations
switch (CHA)
    case 'ch_none', ch=1;
    case 'ch_A',    ch=[0.04 -0.05 0.07 -0.21 -0.5 0.72 0.36 0 0.21 0.03 0.07]';
    case 'ch_B',    ch=[0.407 0.815 0.407]';
    case 'ch_C',    ch=[0.227 0.460 0.688 0.460 0.227]';
    case 'ch_Anna', ch=[0.2 + 0.1i];
    otherwise, disp('Unknown channel!');
end

% select bit-to-symbol mapping
load('signal_constellations');
switch (MOD)
    case 'mo_bpsk_real', sym=sym_bpsk;
    case 'mo_bpsk',      sym=sym_bpsk;
    case 'mo_qpsk',      sym=sym_qpsk;
    case 'mo_8psk',      sym=sym_8psk;
    case 'mo_16qam',     sym=sym_16qam;
    otherwise, disp('Unknown signal constellation!');
end

% select convolutional code
%{
switch (COD)
    case 'co_none',  mem=0; gen=1;
    case 'co_r2_m2', mem=2; gen=[5 7];
    case 'co_r2_m4', mem=4; gen=[19 29];
    case 'co_r2_m6', mem=6; gen=[115 93];
    otherwise, disp('Unknown convolutional code!');
end
%}

% select block length (coded bits) and interleaver
load('interleavers');
switch (BLO)
    case 'bl_noint', Nc=1024;  int=1:1024;
    case 'bl_32',    Nc=32;    int=1:32;
    case 'bl_8',     Nc=8;     int=1:8;
    case 'bl_4',     Nc=4;     int=1:4; 
    case 'bl_16',    Nc=16;    int=1:16;
    case 'bl_128',   Nc=128;   int=1:128; % int=int128;
    case 'bl_1024',  Nc=1024;  int=1:1024; %int=int1024;
    case 'bl_8192',  Nc=8192;  int=int8192;
    case 'bl_65536', Nc=65536; int=int65536;
    otherwise, disp('Unknown block length!');
end

% select puncturing pattern
%{
if (strcmp(COD,'co_none'))
    prate=1; 
    punc=1:Nc;
else
%}
switch (PUN)
    case 'pu_1',       prate=1;    punc=1:Nc;
    case 'pu_100',     prate=6/4;  punc=sort([1:2:Nc*prate 2:6:Nc*prate]);
                                   %rate 3/4; puncturing pattern 111/100
    case 'pu_1000000', prate=14/8; punc=sort([1:2:Nc*prate 2:14:Nc*prate]);
                                   %rate 7/8; puncturing pattern 1111111/1000000
    otherwise, disp('Unknown puncturing pattern!');
end
%end












