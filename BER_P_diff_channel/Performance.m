% List of modulation schemes
modSchemes = {'BPSK', 'QPSK', '16QAM', '64QAM'};
SNR_dB = 0:2:20; % SNR range
numBits = 1e5;

% Preallocate for BERs
BER_results = zeros(length(modSchemes), length(SNR_dB));

% Loop over channels
channels = {'AWGN', 'Rayleigh', 'Rician'};
bestMod = cell(1, length(channels));

for c = 1:length(channels)
    ch = channels{c};
    fprintf('\nChannel: %s\n', ch);

    for m = 1:length(modSchemes)
        modType = modSchemes{m};
        % Get BER array for this modulation in this channel
        BER = simulate_modulation(modType, ch, SNR_dB, numBits);
        BER_results(m, :) = BER;

        fprintf('%s Final BER: %.4e\n', modType, BER(end));
    end

    % Find best modulation (lowest BER at highest SNR)
    [~, idx] = min(BER_results(:, end));
    bestMod{c} = modSchemes{idx};
    fprintf('=> Best Modulation for %s Channel: %s\n', ch, modSchemes{idx});
end

% Display Summary
fprintf('\n===== Best Modulation per Channel =====\n');
for c = 1:length(channels)
    fprintf('%s: %s\n', channels{c}, bestMod{c});
end
