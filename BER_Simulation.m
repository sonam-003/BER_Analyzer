clc;
%% Simulation Parameters
numBits = 1e6;              % Number of bits to transmit
SNR_dB = 0:2:20;            % SNR range in dB
modulationSchemes = {'BPSK', 'QPSK', '16QAM', '64QAM'};
colors = {'b', 'r', 'g', 'm'};
markers = {'o', 's', 'd', '^'};
lineStyles = {'-', '--'};

%% Initialize BER storage
ber_awgn_sim = zeros(length(modulationSchemes), length(SNR_dB));
ber_rayleigh_sim = zeros(length(modulationSchemes), length(SNR_dB));
ber_rician_sim = zeros(length(modulationSchemes), length(SNR_dB));
ber_awgn_theory = zeros(length(modulationSchemes), length(SNR_dB));
ber_rayleigh_theory = zeros(length(modulationSchemes), length(SNR_dB));
ber_rician_theory = zeros(length(modulationSchemes), length(SNR_dB));

%% Main Simulation Loop
for modIdx = 1:length(modulationSchemes)
    modScheme = modulationSchemes{modIdx};
    fprintf('Processing %s...\n', modScheme);
    
    % Generate random binary data
    data = randi([0 1], numBits, 1);
    
    for snrIdx = 1:length(SNR_dB)
        SNR = SNR_dB(snrIdx);
        
        %% Modulation
    switch modScheme
        case 'BPSK'
            modData = 2*data - 1; % BPSK modulation: 0->-1, 1->1
            M = 2;
            
        case 'QPSK'
            % Ensure even number of bits
            extra_bits = mod(length(data), 2);
            if extra_bits ~= 0
                data = data(1:end-extra_bits);
            end
            groupedData = reshape(data, 2, [])';
            symbols = bi2de(groupedData);
            modData = exp(1j*(2*pi*symbols/4 + pi/4));
            M = 4;
            
        case '16QAM'
            % Ensure multiple of 4 bits
            extra_bits = mod(length(data), 4);
            if extra_bits ~= 0
                data = data(1:end-extra_bits);
            end
            groupedData = reshape(data, 4, [])';
            symbols = bi2de(groupedData);
            modData = qammod(symbols, 16, 'UnitAveragePower', true);
            M = 16;
            
        case '64QAM'
            % Ensure multiple of 6 bits
            extra_bits = mod(length(data), 6);
            if extra_bits ~= 0
                data = data(1:end-extra_bits);
            end
            groupedData = reshape(data, 6, [])';
            symbols = bi2de(groupedData);
            modData = qammod(symbols, 64, 'UnitAveragePower', true);
            M = 64;
    end
        
        %% AWGN Channel
        rx_awgn = awgn(modData, SNR, 'measured');
        
        % Demodulation
        switch modScheme
            case 'BPSK'
                demodData = real(rx_awgn) > 0;
            case 'QPSK'
                demodSymbols = round(angle(rx_awgn)*4/(2*pi));
                demodSymbols = mod(demodSymbols,4);
                demodData = reshape(de2bi(demodSymbols)', [], 1);
            case {'16QAM', '64QAM'}
                demodSymbols = qamdemod(rx_awgn, M, 'UnitAveragePower', true);
                demodData = reshape(de2bi(demodSymbols)', [], 1);
        end
        
        % Calculate BER
        [~, ber_awgn_sim(modIdx, snrIdx)] = biterr(data, demodData);
        
        % Theoretical BER for AWGN
        if strcmp(modScheme, 'BPSK')
            ber_awgn_theory(modIdx, snrIdx) = qfunc(sqrt(2*10^(SNR/10)));
        elseif strcmp(modScheme, 'QPSK')
            ber_awgn_theory(modIdx, snrIdx) = qfunc(sqrt(2*10^(SNR/10)));
        else
            % Theoretical BER for AWGN
            switch modScheme
                case 'BPSK'
                    ber_awgn_theory(modIdx, snrIdx) = qfunc(sqrt(2*10^(SNR/10)));
                case 'QPSK'
                    % For QPSK, we can use the BPSK formula since it has the same BER vs Eb/N0
                    ber_awgn_theory(modIdx, snrIdx) = qfunc(sqrt(2*10^(SNR/10)));
                case '16QAM'
                    ber_awgn_theory(modIdx, snrIdx) = berawgn(SNR, 'qam', 16);
                case '64QAM'
                    ber_awgn_theory(modIdx, snrIdx) = berawgn(SNR, 'qam', 64);
            end
         end
        
        %% Rayleigh Fading Channel
        h_rayleigh = (randn(size(modData)) + 1j*randn(size(modData)))/sqrt(2);
        rx_rayleigh = h_rayleigh .* modData + awgn(zeros(size(modData)), SNR, 'measured');
        
        % Equalization (Zero Forcing)
        rx_rayleigh_eq = rx_rayleigh ./ h_rayleigh;
        
        % Demodulation
        switch modScheme
            case 'BPSK'
                demodData = real(rx_rayleigh_eq) > 0;
            case 'QPSK'
                demodSymbols = round(angle(rx_rayleigh_eq)*4/(2*pi));
                demodSymbols = mod(demodSymbols,4);
                demodData = reshape(de2bi(demodSymbols)', [], 1);
            case {'16QAM', '64QAM'}
                demodSymbols = qamdemod(rx_rayleigh_eq, M, 'UnitAveragePower', true);
                demodData = reshape(de2bi(demodSymbols)', [], 1);
        end
        
        % Calculate BER
        [~, ber_rayleigh_sim(modIdx, snrIdx)] = biterr(data, demodData);
        
        % Theoretical BER for Rayleigh
        if strcmp(modScheme, 'BPSK')
            gamma_b = 10^(SNR/10);
            ber_rayleigh_theory(modIdx, snrIdx) = 0.5*(1 - sqrt(gamma_b/(1+gamma_b)));
        else
            % Map modulation scheme names to what berfading expects
            switch modScheme
                case 'QPSK'
                    modType = 'psk';
                    M = 4;
                case '16QAM'
                    modType = 'qam';
                    M = 16;
                case '64QAM'
                    modType = 'qam';
                    M = 64;
            end
            ber_rayleigh_theory(modIdx, snrIdx) = berfading(SNR, modType, M, 1);
        end
        
        %% Rician Fading Channel
        K = 3; % Rician K-factor
        h_rician = sqrt(K/(K+1)) + sqrt(1/(K+1))*(randn(size(modData)) + 1j*randn(size(modData)))/sqrt(2);
        rx_rician = h_rician .* modData + awgn(zeros(size(modData)), SNR, 'measured');
        
        % Equalization (Zero Forcing)
        rx_rician_eq = rx_rician ./ h_rician;
        
        % Demodulation
        switch modScheme
            case 'BPSK'
                demodData = real(rx_rician_eq) > 0;
            case 'QPSK'
                demodSymbols = round(angle(rx_rician_eq)*4/(2*pi));
                demodSymbols = mod(demodSymbols,4);
                demodData = reshape(de2bi(demodSymbols)', [], 1);
            case {'16QAM', '64QAM'}
                demodSymbols = qamdemod(rx_rician_eq, M, 'UnitAveragePower', true);
                demodData = reshape(de2bi(demodSymbols)', [], 1);
        end
        
        % Calculate BER
        [~, ber_rician_sim(modIdx, snrIdx)] = biterr(data, demodData);
        
        % Theoretical BER for Rician (approximation)
        if strcmp(modScheme, 'BPSK')
            gamma_b = 10^(SNR/10);
            X = sqrt(K*gamma_b/(K + gamma_b));
            if exist('marcumq', 'file')
                ber_rician_theory(modIdx, snrIdx) = marcumq(sqrt(2*K), X*sqrt(2), 1) / 2;
            else
                ber_rician_theory(modIdx, snrIdx) = 0.5*(1 - sqrt(K*gamma_b/(1+K+gamma_b)));
            end
        else
            % Map modulation scheme names to what berfading expects
            switch modScheme
                case 'QPSK'
                    modType = 'psk';
                    M = 4;
                case '16QAM'
                    modType = 'qam';
                    M = 16;
                case '64QAM'
                    modType = 'qam';
                    M = 64;
            end
            if exist('berfading', 'file')
                ber_rician_theory(modIdx, snrIdx) = berfading(SNR, modType, M, K);
            else
                ber_rician_theory(modIdx, snrIdx) = ber_awgn_theory(modIdx, snrIdx) * (1 + K)/(1 + K + 10^(SNR/10));
            end
        end
    end
end

%% Plotting Results
figure('Position', [100, 100, 1200, 800]);

% AWGN Channel
subplot(1,3,1);
hold on; grid on;
for modIdx = 1:length(modulationSchemes)
    semilogy(SNR_dB, ber_awgn_sim(modIdx,:), [colors{modIdx} markers{modIdx} '-'], 'LineWidth', 1.5, 'MarkerSize', 8);
    semilogy(SNR_dB, ber_awgn_theory(modIdx,:), [colors{modIdx} '--'], 'LineWidth', 1.5);
end
title('AWGN Channel');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
set(gca, 'YScale', 'log');
ylim([1e-6 1]);
xlim([0 20]);

% Rayleigh Channel
subplot(1,3,2);
hold on; grid on;
for modIdx = 1:length(modulationSchemes)
    semilogy(SNR_dB, ber_rayleigh_sim(modIdx,:), [colors{modIdx} markers{modIdx} '-'], 'LineWidth', 1.5, 'MarkerSize', 8);
    semilogy(SNR_dB, ber_rayleigh_theory(modIdx,:), [colors{modIdx} '--'], 'LineWidth', 1.5);
end
title('Rayleigh Fading Channel');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
set(gca, 'YScale', 'log');
ylim([1e-6 1]);
xlim([0 20]);

% Rician Channel
subplot(1,3,3);
hold on; grid on;
for modIdx = 1:length(modulationSchemes)
    semilogy(SNR_dB, ber_rician_sim(modIdx,:), [colors{modIdx} markers{modIdx} '-'], 'LineWidth', 1.5, 'MarkerSize', 8);
    semilogy(SNR_dB, ber_rician_theory(modIdx,:), [colors{modIdx} '--'], 'LineWidth', 1.5);
end
title('Rician Fading Channel (K=3)');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
set(gca, 'YScale', 'log');
ylim([1e-6 1]);
xlim([0 20]);

% Create legend only on the first subplot
subplot(1,3,1);
legendEntries = cell(1, 2*length(modulationSchemes)); % Pre-allocate cell array
for modIdx = 1:length(modulationSchemes)
    legendEntries{2*modIdx-1} = [modulationSchemes{modIdx} ' Sim'];
    legendEntries{2*modIdx} = [modulationSchemes{modIdx} ' Theory'];
end
legend(legendEntries, 'Location', 'southwest', 'FontSize', 9);

% Adjust subplot spacing
set(gcf, 'Color', 'w');
ha = findobj(gcf, 'type', 'axes');
set(ha, 'YScale', 'log');
