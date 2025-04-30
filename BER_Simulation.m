clc;
%% Parameters
numBits = 1e6;
SNR_dB = 0:1:20;
modulationSchemes = {'BPSK', 'QPSK', '16QAM', '64QAM'};
colors = {'b', 'r', 'g', 'm'};
markers = {'o', 's', 'd', '^'};

%% Initialize BER arrays
ber_awgn_sim = zeros(length(modulationSchemes), length(SNR_dB));
ber_awgn_theory = zeros(length(modulationSchemes), length(SNR_dB));

%% Loop over modulations
for modIdx = 1:length(modulationSchemes)
    modScheme = modulationSchemes{modIdx};
    fprintf('Processing %s...\n', modScheme);
    
    for snrIdx = 1:length(SNR_dB)
        SNR = SNR_dB(snrIdx);
        
        switch modScheme
            case 'BPSK'
                data = randi([0 1], numBits, 1);
                tx = 2*data - 1;
                rx = awgn(tx, SNR, 'measured');
                rxBits = real(rx) > 0;
                ber_awgn_theory(modIdx, snrIdx) = qfunc(sqrt(2*10^(SNR/10)));
            case 'QPSK'
                data = randi([0 1], numBits, 1);
                if mod(length(data),2) ~= 0
                    data = [data; 0];
                end
                symbols = bi2de(reshape(data, 2, []).');
                tx = exp(1j * (pi/4 + pi/2 * symbols));
                rx = awgn(tx, SNR, 'measured');
                angles = angle(rx);
                angles(angles < 0) = angles(angles < 0) + 2*pi;
                receivedSymbols = mod(round((angles - pi/4)/(pi/2)), 4);
                rxBits = reshape(de2bi(receivedSymbols, 2).', [], 1);
                ber_awgn_theory(modIdx, snrIdx) = qfunc(sqrt(2*10^(SNR/10)));
            case '16QAM'
                M = 16;
                bitsPerSymbol = log2(M);
                data = randi([0 1], numBits, 1);
                if mod(length(data), bitsPerSymbol) ~= 0
                    data = [data; zeros(bitsPerSymbol - mod(length(data), bitsPerSymbol), 1)];
                end
                symbols = bi2de(reshape(data, bitsPerSymbol, []).');
                tx = qammod(symbols, M, 'UnitAveragePower', true);
                rx = awgn(tx, SNR, 'measured');
                rxSymbols = qamdemod(rx, M, 'UnitAveragePower', true);
                rxBits = reshape(de2bi(rxSymbols, bitsPerSymbol).', [], 1);
                ber_awgn_theory(modIdx, snrIdx) = berawgn(SNR, 'qam', M);
            case '64QAM'
                M = 64;
                bitsPerSymbol = log2(M);
                data = randi([0 1], numBits, 1);
                if mod(length(data), bitsPerSymbol) ~= 0
                    data = [data; zeros(bitsPerSymbol - mod(length(data), bitsPerSymbol), 1)];
                end
                symbols = bi2de(reshape(data, bitsPerSymbol, []).');
                tx = qammod(symbols, M, 'UnitAveragePower', true);
                rx = awgn(tx, SNR, 'measured');
                rxSymbols = qamdemod(rx, M, 'UnitAveragePower', true);
                rxBits = reshape(de2bi(rxSymbols, bitsPerSymbol, 'left-msb').', [], 1);
                ber_awgn_theory(modIdx, snrIdx) = berawgn(SNR, 'qam', M);
        end

        rxBits = rxBits(1:length(data));
        [~, ber_awgn_sim(modIdx, snrIdx)] = biterr(data, rxBits);
    end
end

%% Plot (Linear Y-axis to match the sample image)
figure('Position', [100 100 800 600]); grid on; hold on;
title('BER Performance in AWGN Channel');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
ylim([0 1]);
xlim([min(SNR_dB), max(SNR_dB)]);

for modIdx = 1:length(modulationSchemes)
    M = 2^modIdx; % BPSK=2, QPSK=4, etc.
    plot(SNR_dB, ber_awgn_sim(modIdx,:), [colors{modIdx} markers{modIdx} '-'], ...
         'LineWidth', 1.5, 'MarkerSize', 6, ...
         'DisplayName', ['Sim ' modulationSchemes{modIdx} ' M=' num2str(M)]);
    plot(SNR_dB, ber_awgn_theory(modIdx,:), [colors{modIdx} '--'], ...
         'LineWidth', 1.5, ...
         'DisplayName', ['Theory ' modulationSchemes{modIdx} ' M=' num2str(M)]);
end

legend('Location', 'northeast', 'FontSize', 9);
set(gcf, 'Color', 'w');
