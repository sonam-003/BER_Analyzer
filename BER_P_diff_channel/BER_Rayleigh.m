 clc;

% Parameters
SNR_dB = 0:2:20;
numBits = 1e6;
modulations = {'BPSK', 'QPSK', '16QAM', '64QAM'};
M_vals = [2, 4, 16, 64];
colors = {'b', 'r', 'g', 'm'};
line_styles = {'-o', '-s', '-d', '-^'};

figure; hold on; grid on;
title('BER Performance over Rayleigh Channel');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
ylim([0 0.5]);

for mod = 1:length(modulations)
    M = M_vals(mod);
    k = log2(M);
    ber_sim = zeros(size(SNR_dB));
    ber_theory = zeros(size(SNR_dB));

    for idx = 1:length(SNR_dB)
        snr_dB = SNR_dB(idx);
        snr_linear = 10^(snr_dB/10);

        % Transmit Bits
        tx_bits = randi([0 1], numBits, 1);

        % Modulation
        switch modulations{mod}
            case 'BPSK'
                tx_symbols = 2*tx_bits - 1;
            case 'QPSK'
                bits_reshaped = reshape(tx_bits(1:floor(numBits/2)*2), [], 2);
                I = 2*bits_reshaped(:,1) - 1;
                Q = 2*bits_reshaped(:,2) - 1;
                tx_symbols = I + 1i*Q;
            case '16QAM'
                bits_reshaped = reshape(tx_bits(1:floor(numBits/4)*4), [], 4);
                symbols = bi2de(bits_reshaped);
                tx_symbols = qammod(symbols, 16, 'UnitAveragePower', true);
            case '64QAM'
                bits_reshaped = reshape(tx_bits(1:floor(numBits/6)*6), [], 6);
                symbols = bi2de(bits_reshaped);
                tx_symbols = qammod(symbols, 64, 'UnitAveragePower', true);
        end

        if strcmp(modulations{mod}, 'BPSK') || strcmp(modulations{mod}, 'QPSK')
            tx_symbols = tx_symbols / sqrt(mean(abs(tx_symbols).^2));
        end

        % Rayleigh Channel + Noise
        h = (randn(size(tx_symbols)) + 1i*randn(size(tx_symbols))) / sqrt(2);
        noise = sqrt(1/(snr_linear*k))*(randn(size(tx_symbols)) + 1i*randn(size(tx_symbols)));
        rx_symbols = h .* tx_symbols + noise;
        rx_symbols = rx_symbols ./ h;

        % Demodulation
        switch modulations{mod}
            case 'BPSK'
                rx_bits = real(rx_symbols) > 0;
            case 'QPSK'
                rx_bits_i = real(rx_symbols) > 0;
                rx_bits_q = imag(rx_symbols) > 0;
                rx_bits = reshape([rx_bits_i, rx_bits_q].', [], 1);
            case '16QAM'
                rx_symbols_demod = qamdemod(rx_symbols, 16, 'UnitAveragePower', true);
                rx_bits = reshape(de2bi(rx_symbols_demod, 4).', [], 1);
            case '64QAM'
                rx_symbols_demod = qamdemod(rx_symbols, 64, 'UnitAveragePower', true);
                rx_bits = reshape(de2bi(rx_symbols_demod, 6).', [], 1);
        end

        % BER
        N = min(length(tx_bits), length(rx_bits));
        ber_sim(idx) = sum(tx_bits(1:N) ~= rx_bits(1:N)) / N;

        % Theoretical BER (Rayleigh)
        switch M
            case 2
                ber_theory(idx) = 0.5 * (1 - sqrt(snr_linear / (1 + snr_linear)));
            case 4
                ber_theory(idx) = 0.5 * (1 - sqrt(snr_linear / (2 + snr_linear)));
            case 16
                ber_theory(idx) = 3/8 * (1 - sqrt(snr_linear / (10 + snr_linear)));
            case 64
                ber_theory(idx) = 7/24 * (1 - sqrt(snr_linear / (42 + snr_linear)));
        end
    end

    % Plot
    plot(SNR_dB, ber_sim, line_styles{mod}, 'Color', colors{mod}, ...
         'DisplayName', ['Sim ' modulations{mod}]);
    plot(SNR_dB, ber_theory, '--', 'Color', colors{mod}, ...
         'DisplayName', ['Theory ' modulations{mod}]);
end

legend('Location', 'northeast');
set(gcf, 'Color', 'w');
