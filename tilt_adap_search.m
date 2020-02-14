%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%uplink%%%%%%%%%%%%%%%%%%%%%%%%%%%
rate_storexx = zeros(Np, 2);
for nnxx = 1 : Np
    thetaxx = nnxx * 0.1 * pi / 180;
    sum_rate_txx = 0;
    sum_rate_tdxx = 0;
    for j = 1 : L
        fx = 0;
        for l = 1 : L
            for k = 1 : K
                Aaz = -min(12 * phi(k, (j-1)*L+l)^2 / phi3dB^2, Am);
                Ael = -min(12 * (theta(k, (j-1)*L+l) - thetaxx)^2 / theta3dB^2, SLAV);
                aforti((j-1)*L*K+(l-1)*K+k, 1) = -min(-Aaz-Ael, Am);
            end
        end
        sumab = zeros(K, 1);
        for k = 1 : K
            for l = 1 : L
                if l == j
                else
                    sumab(k, 1) = sumab(k, 1) + 10^(0.2 * (aforti((j-1)*L*K+(l-1)*K+k, 1)-aforti((j-1)*L*K+(j-1)*K+k, 1))) * betabeta((j-1)*L*K+(l-1)*K+k, 1);
                end
            end
        end
        for k = 1 : K
            fx = fx + log2(1 + 1 / sumab(k, 1));
        end
        sum_rate_txx = sum_rate_txx + fx;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%downlink%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for l = 1 : L
        for ll = 1 : L
            for k = 1 : K
                Aaz = -min(12 * phi(k, (l-1)*L+ll)^2 / phi3dB^2, Am);
                Ael = -min(12 * (theta(k, (l-1)*L+ll) - thetaxx)^2 / theta3dB^2, SLAV);
                aforti((l-1)*L*K+(ll-1)*K+k, 1) = -min(-Aaz-Ael, Am);
            end
        end
    end
    fx = 0;
    sumabb = zeros(K, L);
    for k = 1 : K
        for l = 1 : L%j
            for ll = 1 : L%l
                if l == ll
                else
                    sumabb(k, l) = sumabb(k, l) + 10^(0.1 * (aforti((ll-1)*L*K+(l-1)*K+k, 1)-aforti((ll-1)*L*K+(ll-1)*K+k, 1))) * betaq((ll-1)*L*K+(l-1)*K+k, 1);
                end
            end
        end
    end
    for l = 1 : L%j
        for k = 1 : K
            fx = fx + log2(1 + 1 / sumabb(k, l));
        end
    end
    sum_rate_tdxx = sum_rate_tdxx + fx;
    rate_storexx(nnxx, 1) = sum_rate_txx;
    rate_storexx(nnxx, 2) = sum_rate_tdxx;
end
[~, nnxx] = max(rate_storexx(:,1));
upetheta_s(:, 3) = nnxx * 0.1 * pi / 180;
[~, nnxx] = max(rate_storexx(:,2));
downetheta_s(:, 3) = nnxx * 0.1 * pi / 180;



