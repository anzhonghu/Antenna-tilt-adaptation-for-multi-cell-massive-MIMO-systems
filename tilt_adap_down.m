betaq = zeros(K*L*L, 1);
for l1 = 1 : L%BS
    for l2 = 1 : L%user
        for k = 1 : K
            Aaz = -min(12 * phi(k, (l1-1)*L+l2)^2 / phi3dB^2, Am);
            Ael = -min(12 * (theta(k, (l1-1)*L+l2) - upetheta_s(l1, 1))^2 / theta3dB^2, SLAV);
            D0(k, 1) = -min(-Aaz-Ael, Am);
            D0(k, 1) = 10^(D0(k, 1)*0.1);
            w((l1-1)*L*K+(l2-1)*K+k, 1) = theta3dB * sqrt(min(Aaz+Am, SLAV) / 12);
        end
        Dq(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) =  D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) / D(:,(l2-1)*L*K+(l2-1)*K+1:(l2-1)*L*K+l2*K) * diag(sqrt(D0));
    end
    thetan(l1, 1) = max(theta(1:K, (l1-1)*L+l1) - w((l1-1)*L*K+(l1-1)*K+1:(l1-1)*L*K+l1*K, 1));
    thetan_s1(l1, 1) = thetan(l1, 1);
end
for j = 1 : L
    for l = 1 : L
        for k = 1 : K
            betaq((j-1)*L*K+(l-1)*K+k, 1) = (Dq(k,(j-1)*L*K+(l-1)*K+k) * D(k,(j-1)*L*K+(l-1)*K+k) / Dq(k,(j-1)*L*K+(j-1)*K+k) / D(k,(j-1)*L*K+(j-1)*K+k))^2;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%j for l, l for j%%%%%%%%%%%%%%
for j = 1 : L
    wmin = min(w((j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K, 1));
    detathetan = wmin;
    flag = 1;
    while detathetan > 1e-1
        detathetan = 2 * detathetan;
        if flag > 0
            fx = 0;
            for l = 1 : L
                for ll = 1 : L
                    for k = 1 : K
                        Aaz = -min(12 * phi(k, (l-1)*L+ll)^2 / phi3dB^2, Am);
                        Ael = -min(12 * (theta(k, (l-1)*L+ll) - thetan(l, 1))^2 / theta3dB^2, SLAV);
                        aforti((l-1)*L*K+(ll-1)*K+k, 1) = -min(-Aaz-Ael, Am);
                    end
                end
            end
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
            flag = 0;
        else
        end
        deltar = -1;
        while deltar <= -1e-5
            detathetan = detathetan * 0.5;
            for l = 1 : L
                for ll = 1 : L
                    for k = 1 : K
                        Aaz = -min(12 * phi(k, (l-1)*L+ll)^2 / phi3dB^2, Am);
                        Ael = -min(12 * (theta(k, (l-1)*L+ll) - thetan(l, 1))^2 / theta3dB^2, SLAV);
                        aforti((l-1)*L*K+(ll-1)*K+k, 1) = -min(-Aaz-Ael, Am);
                    end
                end
            end
            sumabb = zeros(K, L);
            sumabbd = zeros(K, L);
            for k = 1 : K
                if abs(theta(k, (j-1)*L+j) - thetan(j, 1)) <= w((j-1)*L*K+(j-1)*K+k, 1)+1e-6
                    flagparti = 1;
                else
                    flagparti = 0;
                end
                for l = 1 : L%j
                    for ll = 1 : L%l
                        if l == ll
                        else
                            sumabb(k, l) = sumabb(k, l) + aforti((ll-1)*L*K+(l-1)*K+k, 1) * betaq((ll-1)*L*K+(l-1)*K+k, 1);
                            if ll == j
                                if abs(theta(k, (j-1)*L+l) - thetan(j, 1)) <= w((j-1)*L*K+(l-1)*K+k, 1)+1e-6
                                    sumabbd(k, 1) = sumabbd(k, 1) + aforti((j-1)*L*K+(l-1)*K+k, 1) * betaq((j-1)*L*K+(l-1)*K+k, 1) * 0.2 * log(10) * 24 / theta3dB^2 ...
                                        * ((theta(k, (j-1)*L+l)-thetan(j, 1)) - (theta(k, (j-1)*L+j)-thetan(j, 1)) * flagparti);
                                else
                                    sumabbd(k, 1) = sumabbd(k, 1) + aforti((j-1)*L*K+(l-1)*K+k, 1) * betaq((j-1)*L*K+(l-1)*K+k, 1) * 0.2 * log(10) * 24 / theta3dB^2 ...
                                        * (0 - (theta(k, (j-1)*L+j)-thetan(j, 1)) * flagparti);
                                end
                            else
                            end
                        end
                    end
                end
            end
            partid = 0;
            for k = 1 : K
                for l = 1 : L%j
                    partid = partid - 1 / (1 + 1 / sumabb(k, l)) / sumabb(k, l)^2 * sumabbd(k, l);
                end
            end
            partid = partid / log(2);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thetaxyz = thetan(j, 1) + sign(partid) * detathetan;
            flagxyz = 0;
            for k = 1 : K
                if abs(theta(k, (j-1)*L+j) - thetaxyz) <= w((j-1)*L*K+(j-1)*K+k, 1)+1e-6
                    flagxyz = 1;
                else
                end
                for l = 1 : L
                    if l == j
                    else
                        if abs(theta(k, (j-1)*L+l) - thetaxyz) <= w((j-1)*L*K+(l-1)*K+k, 1)+1e-6
                            flagxyz = 1;
                        else
                        end
                    end
                end
            end
            if abs(flagxyz-1) < 1e-5
            else
                continue;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thetanew = sign(partid) * detathetan + thetan(j, 1);
            for l = 1 : L
                for k = 1 : K
                    Aaz = -min(12 * phi(k, (j-1)*L+l)^2 / phi3dB^2, Am);
                    Ael = -min(12 * (theta(k, (j-1)*L+l) - thetanew)^2 / theta3dB^2, SLAV);
                    aforti((j-1)*L*K+(l-1)*K+k, 1) = -min(-Aaz-Ael, Am);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fxt = 0;
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
                    fxt = fxt + log2(1 + 1 / sumabb(k, l));
                end
            end
            deltar = fxt - fx - alpha_search * detathetan * abs(partid);
        end
        fx = fxt;
        thetan(j, 1) = thetan(j, 1) + sign(partid) * detathetan;
    end
    downetheta_s(j, 1) = thetan(j, 1);
end