function [Phi,Psi] = HRRIS_fix(R,T,Nr,Nt,N,rho,P_random,Pa_max,K,q,sigma2,Pt,setA)

xi_vec = zeros(N,1);
for n = 1:N
    tn = T(n,:)';
    xi_vec(n) = sigma2 + Pt*norm(tn)^2;
end

Phi = P_random;
Psi = zeros(N,N);

pn_fraction = randi(100,[K,1]);
pn_fraction = pn_fraction./sum(pn_fraction);
pn_vec = Pa_max*pn_fraction;

stop = 0; count = 0;
Rate0 = 0;
while stop == 0 || count < 20    
    
    for n = 1:N
        
        %% Compute An, Bn, Cn, ...
        Rn = R; Rn(:,n) = []; rn = R(:,n);
        Tn = T; Tn(n,:) = []; tn = T(n,:)';
        Phi_n = Phi([1:n-1,n+1:end],[1:n-1,n+1:end]);
        i = 0;
        if ismember(n,setA)
            i = i + 1;
            cols = setA([1:i-1,i+1:K]);
            Rk = R(:,cols); rn = R(:,n);
            Phi_k = Phi(cols,cols);
            
            An = eye(Nr) + (Rk*Phi_k)*(Rk*Phi_k)' + rho*(Rn*Phi_n*Tn)*(Rn*Phi_n*Tn)';
            Bn = rn*rn' + rho*rn*tn'*tn*rn';
            Cn = rn* sum(Rk*Phi_k,2)' + rho*rn*tn' * (Rn*Phi_n*Tn)';
        else
            An = eye(Nr) + rho*(Rn*Phi_n*Tn)*(Rn*Phi_n*Tn)';
            Bn = rho*rn*tn'*tn*rn';
            Cn = rho*rn*tn' * (Rn*Phi_n*Tn)';
        end
        Dn = eye(Nr) + (An^-1)*Bn;
        En = An*Dn;
        Fn = En^(-1)*Cn;
        
        %% Get solutions based on (27)
        [U,S] = eig(Fn);
        [~,i_max] = max(abs(eig(Fn)));
        theta = angle(S(i_max,i_max));
        if ismember(n,setA)
            xi_tilde = sum(pn_vec([1:n-1,n+1:end]));
            pn_vec(n) = Pa_max - xi_tilde;
            Phi(n,n) = sqrt(pn_vec(n)/xi_vec(n))*exp(-1i*theta);
        else
            Phi(n,n) = exp(-1i*theta);
        end
    end
    
    % check convergence
    Rate = real(log2(det( eye(Nr) + rho*(R*Phi*T)*(R*Phi*T)')));
    if abs(Rate - Rate0) < 1e-4
        stop = 1;
    end
    Rate0 = Rate;
    count = count + 1;
end

% quantization
Phi = quantize(Phi,q);

for n = 1:length(setA)
    Psi(setA(n),setA(n)) = Phi(setA(n),setA(n));
end

% check power constraint
if K > 0 && trace(Psi*(Pt*T*T' + sigma2*eye(N))*Psi') - Pa_max > 1e-4
    disp('wrong power');
end

end % eof