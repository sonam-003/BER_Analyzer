% BER vs SNR for BPSK in AWGN, Rayleigh, Rician (Linear Scale)
clear; clc; close all;

% Simulation parameters
SNR_dB = 0:2:20;       % SNR range in dB
numBits = 1e6;         % Number of bits
modType = 'BPSK';      % Modulation scheme

% Channels: AWGN, Rayleigh, Rician
channels = {'AWGN', 'Rayleigh', 'Rician'};
lineStyles = {'-k', '--b', '-.r'};  % Line styles
legendText = cell(1, length(channels));

% Initialize BER matrix
BER_results = zeros(length(channels), length(SNR_dB));

% Simulate for each channel
for c = 1:length(channels)
    ch = channels{c};
    BER_results(c, :) = simulate_modulation_BPSK(modType, ch, SNR_dB, numBits);
    legendText{c} = ch;
end

% Clear previous plots and create new figure
clf;
figure(1);

% Plot BER vs SNR (LINEAR Y-AXIS)
for c = 1:length(channels)
    plot(SNR_dB, BER_results(c,:), lineStyles{c}, 'LineWidth', 2);
    hold on;
end
hold off;

grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Performance of BPSK (Linear Scale)');
legend(legendText, 'Location', 'northeast'); % Changed position for better visibility

% Adjust Y-axis limits (since we're not using log scale)
ylim([0 0.5]); % Adjust based on your expected BER range

% Local function for BPSK simulation
function BER = simulate_modulation_BPSK(modType, channel, SNR_dB, numBits)
    bits = randi([0 1], numBits, 1);  
    symbols = 2*bits - 1;  % BPSK modulation: 0→-1, 1→+1
    
    BER = zeros(size(SNR_dB));
    
    for i = 1:length(SNR_dB)
        switch channel
            case 'AWGN'
                rxSymbols = awgn(symbols, SNR_dB(i), 'measured');
            case 'Rayleigh'
                h = (randn(size(symbols)) + 1i*randn(size(symbols)))/sqrt(2);
                rxSymbols = h.*symbols + awgn(zeros(size(symbols)), SNR_dB(i), 'measured');
            case 'Rician'
                K = 3;  % Rician K-factor
                mean = sqrt(K/(K+1));
                std_dev = sqrt(1/(2*(K+1)));
                h = (mean + std_dev*randn(size(symbols))) + 1i*(mean + std_dev*randn(size(symbols)));
                rxSymbols = h.*symbols + awgn(zeros(size(symbols)), SNR_dB(i), 'measured');
        end
        
        rxBits = real(rxSymbols) > 0;  % BPSK demodulation
        [~, BER(i)] = biterr(bits, rxBits);
    end
end