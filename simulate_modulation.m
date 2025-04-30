% BER vs SNR for BPSK in different channels
clear; clc; close all;

% Simulation parameters
SNR_dB = 0:2:30;       % SNR range in dB
numBits = 1e6;         % Number of bits for simulation
modType = 'BPSK';      % Modulation scheme

% Channel settings
channels = {'AWGN', 'Rayleigh', 'Rician'};
lineStyles = {'-k', '--b', '-.r'};  % Line styles
legendText = {'AWGN (No fading)', 'Rayleigh (Severe fading)', 'Rician (K=3, Moderate fading)'};

% Initialize BER matrix
BER_results = zeros(length(channels), length(SNR_dB));

% Simulate for each channel 
for c = 1:length(channels)
    BER_results(c,:) = simulate_BER(modType, channels{c}, SNR_dB, numBits);
end

% Plot results
figure;
for c = 1:length(channels)
    semilogy(SNR_dB, BER_results(c,:), lineStyles{c}, 'LineWidth', 2);
    hold on;
end
hold off;

grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('BPSK Performance in Different Channels', 'FontSize', 14);
legend(legendText, 'Location', 'southwest');
axis([min(SNR_dB) max(SNR_dB) 1e-6 1]);

% BER simulation function with proper power control
function BER = simulate_BER(modType, channel, SNR_dB, numBits)
    bits = randi([0 1], numBits, 1);
    symbols = 2*bits - 1;  % BPSK modulation
    
    % Normalize symbol energy to 1
    symbols = symbols / sqrt(mean(abs(symbols).^2));
    
    BER = zeros(size(SNR_dB));
    
    for i = 1:length(SNR_dB)
        SNR_linear = 10^(SNR_dB(i)/10);
        
        switch channel
            case 'AWGN'
                % AWGN channel - simple additive noise
                noise_power = 1/SNR_linear;
                noise = sqrt(noise_power/2)*(randn(size(symbols)) + 1i*randn(size(symbols)));
                rxSymbols = symbols + noise;
                rxBits = real(rxSymbols) > 0;
                
            case 'Rayleigh'
                % Rayleigh fading - no dominant path
                h = (randn(size(symbols)) + 1i*randn(size(symbols)))/sqrt(2);
                h = h / sqrt(mean(abs(h).^2)); % Normalize
                
                faded = h .* symbols;
                noise_power = 1/SNR_linear;
                noise = sqrt(noise_power/2)*(randn(size(symbols)) + 1i*randn(size(symbols)));
                rxSymbols = faded + noise;
                
                % Coherent detection with perfect channel estimation
                rxBits = real(conj(h).*rxSymbols) > 0;
                
            case 'Rician'
                % Rician fading - with dominant path (K=3)
                K = 3;
                los = sqrt(K/(K+1));
                scattered = sqrt(1/(2*(K+1)))*(randn(size(symbols)) + 1i*randn(size(symbols)));
                h = los + scattered;
                h = h / sqrt(mean(abs(h).^2)); % Normalize
                
                faded = h .* symbols;
                noise_power = 1/SNR_linear;
                noise = sqrt(noise_power/2)*(randn(size(symbols)) + 1i*randn(size(symbols)));
                rxSymbols = faded + noise;
                
                % Coherent detection with perfect channel estimation
                rxBits = real(conj(h).*rxSymbols) > 0;
        end
        
        [~, BER(i)] = biterr(bits, rxBits);
    end
end