function [Phi,Psi,K_actual] = HRRIS_dynamic(T,N,Pa_max,K,Pt,Phi_passive,q,sigma2,setA,s_max)
%% Algorithm 2
% Note: 
% 1. setA is the set A* obtained from (32a)
% 2. s_max is the set of zeta_n/xi_n associated with Omg


Phi = Phi_passive;
Psi = zeros(N,N);
K_actual = 0;

if K == 1
    tn = T(setA,:)';
    pn_uni = Pa_max/(K*sigma2 + Pt*norm(tn)^2);
    Phi(setA,setA) = sqrt(pn_uni)*Phi(setA,setA);
    Psi(setA,setA) = Phi(setA,setA);
else % power allocation for multiple active elements
    
    %% Find mu satisfying (32b)
    % At first, set mu to its lower bound
    mu0 = K/(Pa_max + sum(s_max));
    mu = mu0;
    Pvec = zeros(N,1);
    for k = 1:length(setA)
        n = setA(k);
        s_n = s_max(k);
        Pn0 = 1/mu - 1/s_n;
        Pvec(n) = max(Pn0,0);
    end
    Ptot = sum(Pvec);
    
    % expanding the range until the equality in (32b) occurs
    mu_max = mu0; mu_min = mu0;
    while Ptot/Pa_max > 1.05 || Ptot/Pa_max < 0.98
        mu_max = mu_max * 1.02;
        mu_min = mu_min / 1.02;
        mu = (mu_min + mu_max)/2;
        Pvec = zeros(N,1);
        for k = 1:length(setA)
            n = setA(k);
            s_n = s_max(k);
            Pvec(n) = max(1/mu - 1/s_n,0);
        end
        Ptot = sum(Pvec);
        if Ptot == 0
            break;
        end
    end
    
    %% obtain solution based on (33)
    for k = 1:length(setA)
        n = setA(k); tn = T(n,:)';
        s_n = s_max(k);
        xi_n = sigma2 + Pt*norm(tn)^2;
        pn = 1/mu - 1/s_n;
        an2 = pn/xi_n;
        if an2 > 1
            Phi(n,n) = sqrt(an2)*Phi(n,n);
            Psi(n,n) = Phi(n,n);
            K_actual = K_actual + 1;
        end
    end
end

% quantization
Phi = quantize(Phi,q);

% check power constraint
if K > 0 && real(trace(Psi*(Pt*T*T' + sigma2*eye(N))*Psi'))/Pa_max > 1.05
    disp('wrong power');
end

end % eof