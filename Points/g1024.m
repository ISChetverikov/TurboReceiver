
e = (4:0.5:7);
f0 = load('bph1024.mat').bph1024(:,1);
f1 = load('bph1024.mat').bph1024(:,end);
f2 = load('bps1024.mat').bps1024(:,end);

f3 = load('sclh1024.mat').sclh1024(:,1);
f4 = load('sclh1024.mat').sclh1024(:,end);
f5 = load('scls01024.mat').scls01024(:,end);
f6 = load('scls21024.mat').scls21024(:,end);

figure();
h=semilogy(e, [f0 f2 f3 f5 f6], '.-','linewidth', 2);
set(gca, 'fontsize', 12); 
grid on;

lstr = {'BP HD, It=20', 'BP soft, It=20', 'SCL8 HD', 'SCL8 soft', 'SCL8 soft, proposed approximation'};
legend(h,lstr,'fontsize',10,'location','southwest');
xlabel('Eb/N0 in dB');
ylabel('Frame error rate');
title('FER performance of Turbo equalization (It=5). P(1024, 512 + CRC24), QPSK');

axis([e(1) e(end)+1e-3 1e-2 1.2]);
%print(gcf,'-djpeg',['result_sim_' settings_str]);