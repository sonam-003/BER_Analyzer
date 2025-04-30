% Define modulation indices and labels
mod_index = [1, 2, 3, 4];  % Index for BPSK, QPSK, 16-QAM, 64-QAM
mod_labels = {'BPSK', 'QPSK', '16-QAM', '64-QAM'};

% Define realistic BERs
ber_awgn     = [1e-3, 7e-4, 0.12, NaN];
ber_rayleigh = [0.5, 0.5, 0.5, 0.5];
ber_rician   = [1e-3, 1e-2, 0.3, NaN];


% Begin plot
figure;
semilogy(mod_index, ber_awgn, '-o', 'LineWidth', 2, ...
    'MarkerSize', 8, 'DisplayName', 'AWGN @ 10dB', 'Color', [0 0.45 0.74]); hold on;
semilogy(mod_index, ber_rayleigh, '-s', 'LineWidth', 2, ...
    'MarkerSize', 8, 'DisplayName', 'Rayleigh @ 10dB', 'Color', [0.85 0.33 0.10]);
semilogy(mod_index, ber_rician, '-d', 'LineWidth', 2, ...
    'MarkerSize', 8, 'DisplayName', 'Rician @ 10dB', 'Color', [0.93 0.69 0.13]);

% Format axes
xticks(mod_index);
xticklabels(mod_labels);
xlabel('Modulation Index', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('BER', 'FontSize', 12, 'FontWeight', 'bold');
ylim([1e-4 1]);
xlim([0.8 4.2]);

% Title and legend
title('BER vs Modulation Order at Different SNR Points', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);

% Grid and axis styles
grid on;
ax = gca;
ax.YScale = 'log';
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.GridAlpha = 0.3;
ax.MinorGridAlpha = 0.1;
ax.GridLineStyle = ':';
