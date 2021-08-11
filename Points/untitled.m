
e = (4:0.5:7);
f0 = load('bph128.mat').fer_bph128(:,1);
f1 = load('bph128.mat').fer_bph128(:,end);
f2 = load('bps128.mat').fer_bps128(:,end);

f3 = load('sch128.mat').fer_sch128(:,1);
f4 = load('sch128.mat').fer_sch128(:,end);
f5 = load('scls1128.mat').fer_scs1128(:,end);
f6 = load('scls2128.mat').fer_scs2128(:,end);

figure();
h=semilogy(e, [f0 f2 f3 f5 f6], '.-','linewidth', 2);
set(gca, 'fontsize', 12); 
grid on;

lstr = {'BP HD, It=20', 'BP soft, It=20', 'SCL8 HD', 'SCL8 soft', 'SCL8 soft, proposed approximation'};
legend(h,lstr,'fontsize',10,'location','southwest');
xlabel('Eb/N0 in dB');
ylabel('Frame error rate');
title('FER performance of Turbo equalization (It=5). P(128, 64 + CRC8), QPSK');

axis([e(1) e(end)+1e-3 1e-2 1.2]);
%print(gcf,'-djpeg',['result_sim_' settings_str]);