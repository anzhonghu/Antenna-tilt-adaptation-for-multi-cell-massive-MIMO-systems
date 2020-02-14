
SNR_s = [-5; 0; 5; 10; 15;];%rho_up, rho_down, in dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = figure;
set(h1,'PaperType','A4');
xx = axes('FontSize',16);
plot(SNR_s, rate_store1(:,2),'b:^','LineWidth',2,'MarkerSize',10)
hold on
plot(SNR_s, rate_store1(:,1),'k-*','LineWidth',2,'MarkerSize',10)
plot(SNR_s, rate_store1(:,3),'r:o','LineWidth',2,'MarkerSize',15)
grid on
plot(SNR_s, rate_store2(:,2),'b:^','LineWidth',2,'MarkerSize',10)
hold on
plot(SNR_s, rate_store2(:,1),'k-*','LineWidth',2,'MarkerSize',10)
plot(SNR_s, rate_store2(:,3),'r:o','LineWidth',2,'MarkerSize',15)
le = legend('Center of MSs', 'Proposed approach', 'Exhaustive search', 'Location','Northwest');
set(le,'Fontsize',16,'Fontname','Times')
set(gca,'XTick',SNR_s)
xlabel('SNR (dB)','Fontsize',20,'Fontname','Times')
ylabel('Spectral efficiency (bps/Hz)','Fontsize',20,'Fontname','Times')
print(h1,'-dpdf','fig_snr')

