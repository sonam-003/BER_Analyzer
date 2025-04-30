% run_BER.m

% Parameters
modType = 'QAM';          % Choose: 'BPSK', 'QPSK', or 'QAM'
M = 16;                   % Modulation order (only for QAM; ignored for BPSK/QPSK)
SNRdB = 0:2:20;           % SNR range in dB
numBits = 1e5;            % Number of bits to transmit
channelType = 'AWGN';     % Choose: 'AWGN', 'Rayleigh', or 'Rician'

% Create object
berAnalyzer = BERAnalyzer(modType, M, SNRdB, numBits, channelType);

% Run analysis
[ber, theoretical] = berAnalyzer.analyzeBER();

% Plot results
berAnalyzer.plotResults(ber, theoretical);
