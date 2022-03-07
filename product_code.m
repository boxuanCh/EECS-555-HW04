%=======================================================%
% %
% This program simulates the performance of a product %
% of a (7,4) Hamming codes. The overall code is a %
% (40,16) product code. %
% %
%=======================================================%
clear all
G=[1 0 0 0 1 0 1;
   0 1 0 0 1 1 1;
   0 0 1 0 1 1 0;
   0 0 0 1 0 1 1]; % Generator Matrix of the (7,4) Hamming code
cb = zeros(16,7);
ebn0db = ((1:17) - 7) * 0.5;
ber = zeros(1,17);
for i=1:16
    b=dec2bin(i-1,4);
    cb(i,:)=mod(b*G,2); % This is a 16 x7 matrix with all possible Hamming code codewords
end
for k=1:17
disp(ebn0db(k));
ebn0=10^(ebn0db(k)/10.0); % Eb/N0 not in dB
EN0=ebn0*16/40; % Energy per coded bit/N0 (not in dB)
E=1; % Arbitrary normalization
N0=1/EN0;
sigma=sqrt(N0/2); % Variance of noise samples
nerrors=0;
ntrials=0;
if (k<13) 
    threshold=2000; 
end %Count 2000 errors for low EbN0
if (k>12) 
    threshold=1000; 
end %Count 1000 errors for high EbN0
while nerrors < threshold
ntrials=ntrials+1;
b=sign(rand(4,4)-0.5); % Generate the data bits
ph=mod(-0.5*(b-1)*G,2); % Generate the horizontal parity
pv=mod(-0.5*(b'-1)*G,2); % Generate the vertical parity
LVext=zeros(4,4); %Initialize the likelihoods with the apriori information (if any)
LHext = zeros(4, 4); % Initialize horizontal likelihood estimation
%=========================================================%
% Generate as 4 by 10 arrary. First 4 cols are info,
% next three are horizontal parity checks and last
% three are vertical parity checks.
%=========================================================%
c=[b 1-2*ph(:,5:7) 1-2*pv(:,5:7)]; % This is the transmitted codeword (4 x 10)
r=sqrt(E)*c+sigma*randn(4,10); % This is the received vector (4 x 10)
LHRCVD=2*r(:,1:7)*sqrt(E)/sigma^2; % This is the extrinisic H information
LVRCVD=2*[r(:,1:4)' r(:,8:10)]*sqrt(E)/sigma^2; % This is the extrinisic V information
LHcode = cb * LHRCVD';
LVcode = cb * LVRCVD';
niterations=2;
for j=1:niterations
%=========================================================%
% Compute The likelihoods based on horizontal %
% parity checks (including extrinisic and apriori) %
%=========================================================%
for i = 1:4
   LHext(i,:) = LLR(LVext(:,i), LHcode(:,i))' - LHRCVD(i,1:4);
end
%=========================================================%
% Compute the likelihoods based on vertical %
% parity checks (including extrinisic and apriori) %
%=========================================================%
for i = 1:4
   LVext(i,:) = LLR(LHext(:,i), LVcode(:,i))' - LVRCVD(i,1:4);
end
end % End of loop for number of iteration
%=========================================================%
% %
% Compute the decisions based on likelihoods %
% %
%=========================================================%
bh=sign(LHRCVD(1:4,1:4)+LVext(1:4,1:4)'+LHext(1:4,1:4)); % This is the final decision
nerrors=nerrors+16-sum(sum(bh==b)); % Count the number of errors
if mod(ntrials,500)==0 % save results periodically in case of crashes
display([ntrials nerrors]); % Display results to see how fast the program runs
save product4 % Save results
semilogy(ebn0db,ber,'r') % Plot partial results
axis([-3 6 1e-5 1])
grid on
end
end % End of loop for the number of trials
disp(ntrials);
ber(k)=nerrors/ntrials/16;
disp(ber(k));
end % End of loop for E_b/N_0
semilogy(ebn0db,ber,'r','LineWidth',2)
axis([-3 6 1e-5 1])
grid on
xlabel('E_b/N_0 (dB)','FontSize',16)
ylabel('BER','FontSize',16)