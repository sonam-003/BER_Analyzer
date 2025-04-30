% BER Analysis with BPSK, QPSK, 16-QAM, 64-QAM over AWGN, Rayleigh, Rician Channels
clear; clc;

SNR_dB = 0:2:20;
numBits = 1e5;
channels = {'AWGN', 'Rayleigh', 'Rician'};
modSchemes = {'BPSK', 'QPSK', '16QAM', '64QAM'};
M_vals = [2, 4, 16, 64];

figure('Position', [100 100 1400 500]);
sgtitle('BER Performance Across Channels and Modulations');

for ch = 1:length(channels)
    subplot(1,3,ch); hold on; grid on;

    for mod = 1:length(modSchemes)
        M = M_vals(mod);
        k = log2(M);
        ber_sim = zeros(size(SNR_dB));
        ber_theory = zeros(size(SNR_dB));

        for idx = 1:length(SNR_dB)
            snr_db = SNR_dB(idx);
            snr_linear = 10^(snr_db/10);

            % --- Transmit Bits ---
            bits = randi([0 1], numBits, 1);

            % --- Modulation ---
            switch M
                case 2  % BPSK
                    tx = 2*bits - 1;
                case 4  % QPSK
                    bits_rs = reshape(bits(1:floor(numBits/2)*2), [], 2);
                    tx = (2*bits_rs(:,1) - 1) + 1i*(2*bits_rs(:,2) - 1);
                otherwise  % QAM
                    bits_rs = reshape(bits(1:floor(numBits/k)*k), [], k);
                    symbols = bi2de(bits_rs);
                    tx = qammod(symbols, M, 'InputType','integer', 'UnitAveragePower',true);
            end

            % --- Normalize (for BPSK, QPSK) ---
            if M <= 4
                tx = tx / sqrt(mean(abs(tx).^2));
            end

            % --- Channel + Noise ---
            EbNo = snr_linear / k;  % Adjust for bits/symbol
            N0 = 1 / EbNo;
            noise = sqrt(N0/2) * (randn(size(tx)) + 1i*randn(size(tx)));

            switch channels{ch}
                case 'AWGN'
                    rx = tx + noise;
                case 'Rayleigh'
                    h = (randn(size(tx)) + 1i*randn(size(tx))) / sqrt(2);
                    rx = h .* tx + noise;
                    rx = rx ./ h;
                case 'Rician'
                    K = 5;
                    s = sqrt(K / (K + 1));
                    sigma = sqrt(1 / (2*(K + 1)));
                    h = s + sigma * (randn(size(tx)) + 1i*randn(size(tx)));
                    rx = h .* tx + noise;
                    rx = rx ./ h;
            end

            % --- Demodulation ---
            switch M
                case 2  % BPSK
                    bits_rx = real(rx) > 0;
                case 4  % QPSK
                    bits_i = real(rx) > 0;
                    bits_q = imag(rx) > 0;
                    bits_rx = reshape([bits_i bits_q].', [], 1);
                otherwise  % QAM
                    symbols_rx = qamdemod(rx, M, 'OutputType','integer', 'UnitAveragePower',true);
                    bits_rx = reshape(de2bi(symbols_rx, k), [], 1);
            end

            % --- BER Calculation ---
            N = min(length(bits), length(bits_rx));
            ber_sim(idx) = sum(bits(1:N) ~= bits_rx(1:N)) / N;

            % --- Theoretical BER ---
            switch channels{ch}
                case 'AWGN'
                    if M == 2
                        ber_theory(idx) = qfunc(sqrt(2*EbNo));
                    elseif M == 4
                        ber_theory(idx) = qfunc(sqrt(2*EbNo));
                    elseif M == 16
                        ber_theory(idx) = (3/8)*erfc(sqrt(0.1*snr_linear));
                    elseif M == 64
                        ber_theory(idx) = (7/24)*erfc(sqrt(0.1*snr_linear));
                    end
                case 'Rayleigh'
                    if M == 2
                        ber_theory(idx) = 0.5*(1 - sqrt(EbNo/(1 + EbNo)));
                    elseif M == 4
                        ber_theory(idx) = 0.5*(1 - sqrt(EbNo/(1 + EbNo)));
                    elseif M == 16
                        ber_theory(idx) = 0.5*(1 - sqrt(EbNo/(10 + EbNo)));
                    elseif M == 64
                        ber_theory(idx) = 0.5*(1 - sqrt(EbNo/(42 + EbNo)));
                    end
                case 'Rician'
                    K = 5;
                    ber_theory(idx) = 0.5 * erfc(sqrt((K*EbNo)/(K+1)));
            end
        end

        % Plotting
        style = {'-o','-s','-d','-^'};
        plot(SNR_dB, ber_sim, style{mod}, 'DisplayName',['Sim ' modSchemes{mod} ' M=' num2str(M)]);
        plot(SNR_dB, ber_theory, '--', 'DisplayName',['Theory ' modSchemes{mod} ' M=' num2str(M)]);
    end

    title(['Channel: ' channels{ch}]);
    xlabel('SNR (dB)');
    ylabel('Bit Error Rate (BER)');
    legend('Location','southwest');
    ylim([0 0.6]);  % Adjusted for linear scale
    % set(gca, 'YScale','log'); % <-- Commented out for linear Y-axis
end
