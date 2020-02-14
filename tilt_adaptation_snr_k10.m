%Semiblind channel estimation for rapid fading multicell MIMO systems
%This file is for the simulation calculation
clear;
close all;
%%constants
L = 3;%cell number
K = 10;%user number2,3
tao = K;%pilot length
eta_updown = 3 / 7;
r = 1600;%center to edge distance(m)
rc = r * 0.8;
rh = 100;%minimum terminal radius of the cell(m)
ra = rc / rh - 1;
gamma_decay = 3.8;%decay exponent
height = 32;
i_ant = 196;
sigma = 8;%in dB
Num = 1e4;%iteration 2e2
Np = 900;
SNR_s = [-5; 0; 5; 10; 15;];%rho_up, rho_down, in dB
%%position of every base
base(1:3,1) = [0;(1i * 2 * rc);(sqrt(3) * rc + 1i * rc);];
Am = 25;
SLAV = 20;
phi3dB = 70/180*pi;
theta3dB = 7/180*pi;
D = zeros(K,K*L*L);
Dq = zeros(K,K*L*L);
As = zeros(K,K*L*L);
H = zeros(i_ant,K*L*L);
G = zeros(i_ant,K*L*L);
Hjjhat_s = zeros(i_ant,K*L);
Fl = zeros(i_ant,K*L);
phi = zeros(K,L*L);
theta = zeros(K,L*L);
pos = zeros(K, L);
D0 = zeros(K, 1);
upetheta = zeros(L, 1);
downetheta = zeros(L, 1);
upetheta_s = zeros(L, 3);
downetheta_s = zeros(L, 3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rate_store = zeros(length(SNR_s), 3);
%Here the iteration begins
for nn = 1 : length(SNR_s)
    sum_rate_t = zeros(1, 3);
    sum_rate_td = zeros(1, 3);
    SNR = SNR_s(nn, 1);
    amp = 10 ^ (SNR*0.05);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1 : Num
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%position of every terminal, unifrom distribute
        dis(1:K,1:L) = (rem(rand(K,L) * ra, ra) + 1) * rh;
        ang(1:K,1) = -rand(K,1) * 2 * pi / 3;
        ang(1:K,2) = rand(K,1) * 2 * pi / 3;
        ang(1:K,3) = (ones(K,1) + rand(K,1)) * 2 * pi / 3;
        pos(1:K,1:L) = dis .* (exp(1i * ang));
        for ll = 1 : L - 1
            pos(:,ll+1) = pos(:,ll+1) + base(ll+1,1);
        end
        shadow_amp = sqrt(10.^(randn(1,K*L) * sigma * 0.1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %channel matrix, large scale fading and small scale fading
        for l1 = 1 : L%BS
            for l2 = 1 : L%user
                H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = 1 / sqrt(2) * (randn(i_ant,K)+1i*randn(i_ant,K));
                x = ((abs(pos(:,l2)-base(l1,1))).^2 + height^2).^(0.5);
                for k = 1 : K
                    pos_temp(1, 1) = real(pos(k,l2) - base(l1,1));
                    pos_temp(1, 2) = imag(pos(k,l2) - base(l1,1));
                    pos_temp(1, 3) = sqrt((abs(pos(k,l2) - base(l1,1)))^2 + height^2);
                    phi(k, (l1-1)*L+l2) = asin(pos_temp(1, 2) / sqrt(pos_temp(1, 2)^2 + (pos_temp(1, 1))^2));%az
                    theta(k, (l1-1)*L+l2) = asin(height / pos_temp(1, 3));%el
                end
                D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) =  diag(((x*0.01).^(-0.5*gamma_decay))) * diag(sqrt(shadow_amp(:,(l2-1)*K+1:l2*K)));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tilt_adap_up;
        tilt_adap_down;
        tilt_adap_search;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for nni = 1 : 3
            for l1 = 1 : L%BS
                switch  nni
                    case 1
                        upetheta(l1, 1) = upetheta_s(l1, 1);%proposed
                        downetheta(l1, 1) = downetheta_s(l1, 1);
                    case 2
                        upetheta(l1, 1) = upetheta_s(l1, 2);%%%perfect CSI
                        downetheta(l1, 1) = downetheta_s(l1, 2);
                    case 3
                        upetheta(l1, 1) = upetheta_s(l1, 3);%%%perfect CSI
                        downetheta(l1, 1) = downetheta_s(l1, 3);
                    otherwise
                end
            end
            medq = 0;
            for l1 = 1 : L%BS
                for l2 = 1 : L%user
                    for k = 1 : K
                        Aaz = -min(12 * phi(k, (l1-1)*L+l2)^2 / phi3dB^2, Am);
                        Ael = -min(12 * (theta(k, (l1-1)*L+l2) - upetheta(l1, 1))^2 / theta3dB^2, SLAV);
                        D0(k, 1) = -min(-Aaz-Ael, Am);
                        D0(k, 1) = 10^(D0(k, 1)*0.1);
                    end
                    Dq(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) =  D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) / D(:,(l2-1)*L*K+(l2-1)*K+1:(l2-1)*L*K+l2*K) * diag(sqrt(D0));
                end
            end
            for k = 1 : K
                medq = medq + Dq(k,(l1-1)*L*K+(l1-1)*K+k);
            end
            amp_use = amp / medq * K * L;%sqrt(rho),or sqrt(rho_d)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %calculate error rate
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%uplink%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j = 1 : L
                Hjjhat = H(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K) + (1 / amp_use) * randn(i_ant, K) / Dq(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K);
                for l = 1 : L
                    if l == j
                    else
                        Hjjhat = Hjjhat + H(:,(j-1)*L*K+(l-1)*K+1:(j-1)*L*K+l*K) * Dq(:,(j-1)*L*K+(l-1)*K+1:(j-1)*L*K+l*K) / Dq(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K);
                    end
                end
                Hjjhat_s(:, (j-1)*K+1:j*K) = Hjjhat;
                T = ((Hjjhat * Dq(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))' * (Hjjhat * Dq(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)) + (1 / amp_use^2) * eye(K))...
                    \ (Hjjhat * Dq(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))';
                for k = 1 : K
                    Sjk_t = amp_use^2 * abs(T(k, :) * H(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K) * Dq(:,(j-1)*L*K+(j-1)*K+k))^2;
                    Injk_t = 0;
                    for ll = 1 : L
                        if ll == j
                            for  kk = 1 : K
                                if kk == k
                                else
                                    Injk_t = Injk_t + amp_use^2 * (abs(T(k, :) * H(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K) * Dq(:,(j-1)*L*K+(j-1)*K+kk)))^2;
                                end
                            end
                        else
                            Injk_t = Injk_t + amp_use^2 * (norm(T(k, :) * H(:,(j-1)*L*K+(ll-1)*K+1:(j-1)*L*K+ll*K) * Dq(:,(j-1)*L*K+(ll-1)*K+1:(j-1)*L*K+ll*K)))^2;
                        end
                    end
                    Injk_t = Injk_t + (norm(T(k, :)))^2;
                    sinr_t = Sjk_t / Injk_t;
                    sum_rate_t(1, nni) = sum_rate_t(1, nni) + log2(1 + sinr_t);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%downlink%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            delta_pow = 0;
            for l1 = 1 : L%BS
                for l2 = 1 : L%user
                    Aaz = -min(12 * phi(k, (l1-1)*L+l2)^2 / phi3dB^2, Am);
                    Ael = -min(12 * (theta(k, (l1-1)*L+l2) - downetheta(l1, 1))^2 / theta3dB^2, SLAV);
                    D0(k, 1) = -min(-Aaz-Ael, Am);
                    D0(k, 1) = 10^(D0(k, 1)*0.1);
                    As(:, (l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = diag(D0);
                    if l1 == l2
                        for k = 1 : K
                            delta_pow = delta_pow +  D(k,(l1-1)*L*K+(l1-1)*K+k) * sqrt(D0(k, 1));
                        end
                    else
                    end
                end
            end
            amp_use_d = amp / delta_pow * K * L;
            epsilong = zeros(L, 1);
            for j = 1 : L
                X_temp = Hjjhat_s(:, (j-1)*K+1:j*K) * sqrt(As(:, (j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)) * D(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K);
                Fl(:, (j-1)*K+1:j*K) = X_temp / (X_temp' * X_temp + K / amp_use_d^2 * eye(K));
                epsilong(j, 1) = amp_use_d / sqrt(abs(trace((As(:, (j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))^(-2) * (Fl(:, (j-1)*K+1:j*K))' * Fl(:, (j-1)*K+1:j*K))));
            end
            for j = 1 : L
                for k = 1 : K
                    Sjk_t = (abs(D(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K) * sqrt(As(:, (j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)) * H(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)' ...
                        * Fl(:, (j-1)*K+k) / As(k, (j-1)*L*K+(j-1)*K+k)))^2;
                    Injk_t = 0;
                    for kk = 1 : K
                        if kk == k
                        else
                            Injk_t = Injk_t + (abs(D(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K) * sqrt(As(:, (j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)) * H(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)' ...
                                * Fl(:, (j-1)*K+kk) / As(kk, (j-1)*L*K+(j-1)*K+kk)))^2;
                        end
                    end
                    for ll = 1 : L
                        if ll == j
                        else
                            Injk_t = Injk_t + (epsilong(ll, 1))^2 / (epsilong(j, 1))^2 * (norm(D(k,(ll-1)*L*K+(j-1)*K+1:(ll-1)*L*K+j*K) ...
                                * sqrt(As(:, (l1-1)*L*K+(j-1)*K+1:(l1-1)*L*K+j*K)) * H(:,(ll-1)*L*K+(j-1)*K+1:(ll-1)*L*K+j*K)' * Fl(:, (ll-1)*K+1:ll*K) ...
                                / As(:, (l1-1)*L*K+(ll-1)*K+1:(l1-1)*L*K+ll*K)))^2;
                        end
                    end
                    Injk_t = Injk_t + 1 / (epsilong(j, 1))^2;
                    sinr_t = Sjk_t / Injk_t;
                    sum_rate_td(1, nni) = sum_rate_td(1, nni) + log2(1 + sinr_t);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sprintf('%d,%d',nn,jj)
    end
    sum_rate_t = sum_rate_t / Num;
    sum_rate_td = sum_rate_td / Num;
    rate_store(nn, :) = (sum_rate_t + sum_rate_td) * eta_updown;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = figure;
set(h1,'PaperType','A4');
xx = axes('FontSize',16);
plot(SNR_s, rate_store(:,2),'b:d','LineWidth',2,'MarkerSize',15)
hold on
plot(SNR_s, rate_store(:,1),'k-*','LineWidth',2,'MarkerSize',10)
plot(SNR_s, rate_store(:,3),'r:s','LineWidth',2,'MarkerSize',15)
grid on
le = legend('Center of MSs', 'Proposed approach', 'Exhaustive search', 'Location','Southeast');
set(le,'Fontsize',16,'Fontname','Times')
set(gca,'XTick',SNR_s)
xlabel('SNR (dB)','Fontsize',20,'Fontname','Times')
ylabel('Spectral efficiency (bps/Hz)','Fontsize',20,'Fontname','Times')
print(h1,'-dpdf','fig_snr')

