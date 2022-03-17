load('BER_AWGN.mat');
semilogy(ebn0db,ber,'r','LineWidth',2)
hold on
load('BER_rayleign_known.mat');
semilogy(ebn0db,ber,'b','LineWidth',2)
hold on
load('BER_rayleign_unknown.mat');
semilogy(ebn0db,ber,'g','LineWidth',2)
hold on
legend(["AWGN", "Rayleigh with side information", "Rayleigh without side information"]);
axis([-3 6 1e-5 1])
grid on
xlabel('E_b/N_0 (dB)','FontSize',16)
ylabel('BER','FontSize',16)