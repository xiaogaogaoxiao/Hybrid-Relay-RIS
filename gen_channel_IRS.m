clear, clc
array_type = 'ULA';

Nt = 32; Nr = 2; N = 200;
xM_vec = 40;

for ii = 1:length(xM_vec)
    xM = xM_vec(ii);
    % distance between source, RIS/Relay, and destination
    %xM = 50;                                 % MS position on x-axis
    yM = 2;                                 % MS position on y-axis
    xH = 51;                                % HR-RIS position on x-axis
    dt = xH;                                % BS-MS distance
    dr = sqrt((xH - xM)^2 + yM^2);          % IRS-MS distance
    alpha_t = 2.2;
    alpha_r = 2.8;
    beta_0 = db2pow(-30);
    
    % antenna gains (linear scale) at the source, relay/IRS, and destination
    Gt = db2pow(0);
    GRIS = db2pow(0);
    
    % pathloss_3GPP_LOS = @(x) db2pow(-28-20*log10(fc)-22*log10(x));
    % pathloss_3GPP_NLOS = @(x) db2pow(-22.7-26*log10(fc)-36.7*log10(x));
    pathloss_t = @(x) beta_0*x^(-alpha_t);
    pathloss_r = @(x) beta_0*x^(-alpha_r);
    
    %Compute the channel gains using the 3GPP models and antenna gains
    beta_t = pathloss_t(dt)*Gt*GRIS;
    beta_r = pathloss_r(dr)*GRIS;
    
    Nchannel = 200;
    
    for n = 1:Nchannel
        % generate 1 channel
        % AP-IRS channel
        kt = 1000000;
        T_NLOS = 1/sqrt(2)*(randn(N,Nt) + 1i*randn(N,Nt)); % Rayleigh fading for NLOS
        T_LOS = LoS_channel(N, Nt, Nr, 1);
        T = sqrt(beta_t) * (sqrt(kt/(kt+1))*T_LOS + sqrt(1/(kt+1))*T_NLOS);
        
        % IRS-MS channel
        kr = 0;
        R_NLOS = 1/sqrt(2)*(randn(Nr,N) + 1i*randn(Nr,N)); % Rayleigh fading for NLOS
        R_LOS = LoS_channel(N, Nt, Nr, 2);
        R = sqrt(beta_r) * (sqrt(kr/(kr+1))*R_LOS + sqrt(1/(kr+1))*R_NLOS);
        
        R_all(:,:,n) = R; % IRS -> MS
        T_all(:,:,n) = T; % BS -> IRS
        P_random_all(:,:,n) = diag(exp(1i*2*pi.*rand(N,1)));
    end
    
    % save data for simulation
    file_name = strcat('.\sim_data\',num2str(Nt),'x',num2str(Nr),'x',num2str(N),'x',num2str(xM),'x',num2str(xH),'.mat');
    save(file_name,'R_all','T_all','P_random_all')
end
