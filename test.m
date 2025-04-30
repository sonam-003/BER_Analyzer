% Test parameters
test_mod = 'QPSK';
test_channel = 'AWGN';
test_SNR = 0:5:20;
test_bits = 1e4;

% Run simulation
BER_results = simulate_modulation(test_mod, test_channel, test_SNR, test_bits);

% Display results
disp('SNR (dB) | BER');
disp('----------------');
for i = 1:length(test_SNR)
    fprintf('%5d    | %.4e\n', test_SNR(i), BER_results(i));
end