% BER Analysis for BPSK, QPSK, 16-QAM, 64-QAM over AWGN Channel
clear; clc;

SNR_dB = 0:2:20;
numBits = 1e5;
modSchemes = {'BPSK', 'QPSK', '16QAM', '64QAM'};
M_vals = [2, 4, 16, 64];

figure('Name', 'AWGN Channel'); hold on; grid on;
title('BER Performance over AWGN Channel');
xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
ylim([0 0.6]);  % Linear scale
legendEntries = {};

for mod = 1:length(modSchemes)
    M = M_vals(mod);
    k = log2(M);
    ber_sim = zeros(size(SNR_dB));
    ber_theory = zeros(size(SNR_dB));

    for idx = 1:length(SNR_dB)
        snr_db = SNR_dB(idx);
        snr_linear = 10^(snr_db/10);

        bits = randi([0 1], numBits, 1);

        switch M
            case 2
                tx = 2*bits - 1;
            case 4
                bits_rs = reshape(bits(1:floor(numBits/2)*2), [], 2);
                tx = (2*bits_rs(:,1) - 1) + 1i*(2*bits_rs(:,2) - 1);
            otherwise
                bits_rs = reshape(bits(1:floor(numBits/k)*k), [], k);
                symbols = bi2de(bits_rs);
                tx = qammod(symbols, M, 'InputType','integer', 'UnitAveragePower',true);
        end

        if M <= 4
            tx = tx / sqrt(mean(abs(tx).^2));
        end

        EbNo = snr_linear / k;
        N0 = 1 / EbNo;
        noise = sqrt(N0/2) * (randn(size(tx)) + 1i*randn(size(tx)));

        rx = tx + noise;

        switch M
            case 2
                bits_rx = real(rx) > 0;
            case 4
                bits_i = real(rx) > 0;
                bits_q = imag(rx) > 0;
                bits_rx = reshape([bits_i bits_q].', [], 1);
            otherwise
                symbols_rx = qamdemod(rx, M, 'OutputType','integer', 'UnitAveragePower',true);
                bits_rx = reshape(de2bi(symbols_rx, k), [], 1);
        end

        N = min(length(bits), length(bits_rx));
        ber_sim(idx) = sum(bits(1:N) ~= bits_rx(1:N)) / N;

        if M == 2 || M == 4
            ber_theory(idx) = qfunc(sqrt(2*EbNo));
        elseif M == 16
            ber_theory(idx) = (3/8)*erfc(sqrt(0.1*snr_linear));
        elseif M == 64
            ber_theory(idx) = (7/24)*erfc(sqrt(0.1*snr_linear));
        end
    end

    plot(SNR_dB, ber_sim, '-o');
    plot(SNR_dB, ber_theory, '--');
    legendEntries{end+1} = ['Sim ' modSchemes{mod}];
    legendEntries{end+1} = ['Theory ' modSchemes{mod}];
end

legend(legendEntries, 'Location','southwest');
