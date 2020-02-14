%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%uplink%%%%%%%%%%%%%%%
thetan = zeros(L, 1);
thetan_s1 = zeros(L, 1);
for j = 1 : L
    thetan(j, 1) = mean(theta(:, (j-1)*L+j));
    upetheta_s(j, 2) = thetan(j, 1);
    downetheta_s(j, 2) = thetan(j, 1);
end
w = zeros(K*L*L, 1);
betabeta = zeros(K*L*L, 1);
aforti = zeros(K*L*L, 1);
alpha_search = 0.5;%%0.5 is chosen between 0 and 0.5, as cover's book,p464
for j = 1 : L
    for l = 1 : L
        for k = 1 : K
            Aaz = -min(12 * phi(k, (j-1)*L+l)^2 / phi3dB^2, Am);
            w((j-1)*L*K+(l-1)*K+k, 1) = theta3dB * sqrt(min(Aaz+Am, SLAV) / 12);
            betabeta((j-1)*L*K+(l-1)*K+k, 1) = (D(k,(j-1)*L*K+(l-1)*K+k) / D(k,(l-1)*L*K+(l-1)*K+k))^4;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%uplink%%%%%%%%%%%%%%%
for j = 1 : L
    wmin = min(w((j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K, 1));
    flag = 1;
    detathetan = wmin;
    while detathetan > 1e-1 * pi / 180
        detathetan = 2 * detathetan;
        if flag > 0
            thetan(j, 1) = max(theta(1:K, (j-1)*L+j) - w((j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K, 1));
            fx = 0;
            for l = 1 : L
                for k = 1 : K
                    Aaz = -min(12 * phi(k, (j-1)*L+l)^2 / phi3dB^2, Am);
                    Ael = -min(12 * (theta(k, (j-1)*L+l) - thetan(j, 1))^2 / theta3dB^2, SLAV);
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
            flag = 0;
        else
        end
        deltar = -1;
        while deltar <= -1e-5
            detathetan = detathetan * 0.5;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for l = 1 : L
                for k = 1 : K
                    Aaz = -min(12 * phi(k, (j-1)*L+l)^2 / phi3dB^2, Am);
                    Ael = -min(12 * (theta(k, (j-1)*L+l) - thetan(j, 1))^2 / theta3dB^2, SLAV);
                    aforti((j-1)*L*K+(l-1)*K+k, 1) = -min(-Aaz-Ael, Am);
                end
            end
            sumab = zeros(K, 1);
            sumabd = zeros(K, 1);
            for k = 1 : K
                if abs(theta(k, (j-1)*L+j) - thetan(j, 1)) <= w((j-1)*L*K+(j-1)*K+k, 1) + 1e-6
                    flagparti = 1;
                else
                    flagparti = 0;
                end
                for l = 1 : L
                    if l == j
                    else
                        sumab(k, 1) = sumab(k, 1) + aforti((j-1)*L*K+(l-1)*K+k, 1) * betabeta((j-1)*L*K+(l-1)*K+k, 1);
                        if abs(theta(k, (j-1)*L+l) - thetan(j, 1)) <= w((j-1)*L*K+(l-1)*K+k, 1) + 1e-6
                            sumabd(k, 1) = sumabd(k, 1) + aforti((j-1)*L*K+(l-1)*K+k, 1) * betabeta((j-1)*L*K+(l-1)*K+k, 1) * 0.2 * log(10) * 24 / theta3dB^2 ...
                                * ((theta(k, (j-1)*L+l)-thetan(j, 1)) - (theta(k, (j-1)*L+j)-thetan(j, 1)) * flagparti);
                        else
                            sumabd(k, 1) = sumabd(k, 1) + aforti((j-1)*L*K+(l-1)*K+k, 1) * betabeta((j-1)*L*K+(l-1)*K+k, 1) * 0.2 * log(10) * 24 / theta3dB^2 ...
                                * (0 - (theta(k, (j-1)*L+j)-thetan(j, 1)) * flagparti);
                        end
                    end
                end
            end
            partid = 0;
            for k = 1 : K
                partid = partid - 1 / (1 + 1 / sumab(k, 1)) / sumab(k, 1)^2 * sumabd(k, 1);
            end
            partid = partid / log(2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thetaxyz = thetan(j, 1) + sign(partid) * detathetan;
            flagxyz = 0;
            for k = 1 : K
                if abs(theta(k, (j-1)*L+j) - thetaxyz) <= w((j-1)*L*K+(j-1)*K+k, 1) + 1e-6
                    flagxyz = 1;
                else
                end
                for l = 1 : L
                    if l == j
                    else
                        if abs(theta(k, (j-1)*L+l) - thetaxyz) <= w((j-1)*L*K+(l-1)*K+k, 1) + 1e-6
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
            thetanew = sign(partid) * detathetan + thetan(j, 1);
            for l = 1 : L
                for k = 1 : K
                        Aaz = -min(12 * phi(k, (j-1)*L+l)^2 / phi3dB^2, Am);
                        Ael = -min(12 * (theta(k, (j-1)*L+l) - thetanew)^2 / theta3dB^2, SLAV);
                        aforti((j-1)*L*K+(l-1)*K+k, 1) = -min(-Aaz-Ael, Am);
                end
            end
            fxt = 0;
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
                fxt = fxt + log2(1 + 1 / sumab(k, 1));
            end
            deltar = fxt - fx - alpha_search * detathetan * abs(partid);
        end
        if abs(partid) < 1e-6
            break;
        else
        end
        thetan(j, 1) = thetan(j, 1) + sign(partid) * detathetan;
        fx = fxt;
    end
    upetheta_s(j, 1) = thetan(j, 1);
end