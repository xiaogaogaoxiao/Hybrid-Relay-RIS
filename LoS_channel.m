function [H] = LoS_channel(N, Nt, Nr, channel)
Nx = 5;

if channel == 1
    AoD = 2*pi * rand;
    AoA_az = 2*pi * rand;
    AoA_el = pi*rand - pi/2;
    
    a_tx = exp(1i * [0:Nt-1] .* pi * sin(AoD)).';
    
    sinsin = sin(AoA_az)*sin(AoA_el);
    sincos = sin(AoA_az)*cos(AoA_el);
    n = [0:N-1];
    a_rx = exp(1i * pi * ( floor(n./Nx).*sinsin + (n-floor(n./Nx)*Nx) * sincos )).';
    
    H = a_rx*a_tx';
else
    H = zeros(Nr,N);
    AoD = 2*pi * rand;
    AoA_az = 2*pi * rand;
    AoA_el = pi*rand - pi/2;
    
    a_rx = exp(1i * [0:Nr-1] .* pi * sin(AoD)).';
    
    sinsin = sin(AoA_az)*sin(AoA_el);
    sincos = sin(AoA_az)*cos(AoA_el);
    n = [0:N-1];
    a_tx = exp(1i * pi * ( floor(n./Nx).*sinsin + (n-floor(n./Nx)*Nx) * sincos )).';
    
    H = a_rx*a_tx';
end
end