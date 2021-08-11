clc; 
clear all;
disp('Turbo equalization simulation ------------------------------------');

%% Paths
addpath('Decoders');
addpath('Data');
addpath('Equalizers');

%% define parameters
nTx = 1;                  % Number of transmit antennas
nRx = 4;                  % Number of receive antennas
rng(512);
H = (randn(nRx, nTx) +  1i*randn(nRx, nTx))/sqrt(2);

% equalizers={'eq_map_det' 'eq_exact_lin' 'eq_approx_lin'};
% equalizer=equalizers{3};

mappings={'bpsk', 'qpsk'};
modulation=mappings{2};

decoder1 = containers.Map('Type', 'softSCL', 'UniformValues', false);
decoder1('L') = 8;
decoder1('BpType') = 0;
decoder1('PolyLength') = 8;
decoder1('CrcPoly') = 'z^8 + z^7 + z^6 + z^4 + z^2 + 1';

decoder2 = containers.Map('Type', 'softSCL', 'UniformValues', false);
decoder2('L') = 8;
decoder2('BpType') = 2;
decoder2('PolyLength') = 8;
decoder2('CrcPoly') = 'z^8 + z^7 + z^6 + z^4 + z^2 + 1';

decoder3 = containers.Map('Type', 'BP', 'UniformValues', false);
decoder3('Iterations') = 20;

decoders = { decoder1 decoder2 };

iterations = 5;

N=128;
K=64;
Ebn0Array = (3:0.5:5);
        
MINERR=100;          % minimum number of errors to quit simulation
MAXSIM=1e3;         % maximum number of code bits per simulation

%% initialiaztion
FER = zeros(length(decoders), length(Ebn0Array), iterations);

%% simulation iterations over ebn0
setting_str = 'settings_dummy';
disp(['Initiate calculation for the setting ' setting_str ':']); 
t0=cputime;

for i = 1:length(Ebn0Array)
    ebn0 = Ebn0Array(i);
    disp(['Start simulation for ebno: ' num2str(ebn0)]);
    
    for j = 1:length(decoders)
        FER(j, i, :) = Simulation(ebn0, N, K, modulation, H, iterations, MINERR, MAXSIM, decoders{j});
    end
end

t1=cputime;
disp(['Calculation finished in ' datestr((t1-t0)/3600/24,'dd HH:MM:SS') '.']);

%% Plotting
iterationsToPlot = [ 1 iterations];
FerToPlot = [];
legendStr = cell(length(iterationsToPlot) * length(decoders), 1);
FerToPlot = [];
for j = 1:length(decoders)
    
    for i = 1:length(iterationsToPlot)
        FerToPlot = [FerToPlot; FER(j,:, iterationsToPlot(i))];
        decoder = decoders{j};
        legendStr{length(iterationsToPlot) * (j - 1) + i} = [decoder('Type') '_' num2str(iterationsToPlot(i))];
    end
end

legendStr = {'Soft SCL8, it=0', 'Soft SCL8, it=5', 'Soft SCL8, proposed approximation, it=0', 'Soft SCL8, proposed approximation, it=5'};
figure();
h = semilogy(Ebn0Array, FerToPlot, '.-', 'linewidth', 2); 
set(gca, 'fontsize', 12); 
grid on;

legend(h, legendStr, 'fontsize', 10, 'location', 'southwest');
xlabel('Eb/N0, dB');
ylabel('Frame error rate');
title('FER for scheme with pilots P(128, 64 + CRC8), QPSK, turbo-reciever(it=5)');
axis([Ebn0Array(1) Ebn0Array(end)+1e-3 1e-3 1.2]);

