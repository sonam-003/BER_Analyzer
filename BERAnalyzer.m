classdef BERAnalyzer
    properties
        modType        % Modulation type (BPSK, QPSK, QAM)
        M             % Modulation order
        SNRdB         % SNR range in dB
        numBits       % Number of bits to simulate
        channelType   % Channel type (AWGN, Rayleigh, Rician)
    end
    
    methods
        function obj = BERAnalyzer(modType, M, SNRdB, numBits, channelType)
            obj.modType = modType;
            obj.M = M;
            obj.SNRdB = SNRdB;
            obj.numBits = numBits;
            obj.channelType = channelType;
        end
        
        % Custom Q-function implementation
        function q = qfunction(~, x)
            q = 0.5 * erfc(x/sqrt(2));
        end
        
        % Binary to decimal conversion
        function d = binary2decimal(~, b)
            d = sum(b .* (2.^(size(b, 2)-1:-1:0)), 2);
        end
        
        % Custom QAM constellation mapping
        function symbols = qamMapping(~, indices, M)
            % Calculate number of points per dimension
            L = sqrt(M);
            
            % Generate constellation points
            alphaI = -(L-1):2:L-1;
            alphaQ = -(L-1):2:L-1;
            
            % Create constellation mapping
            [I, Q] = meshgrid(alphaI, alphaQ);
            constellation = (I(:) + 1j*Q(:))/sqrt(10); % Normalize power
            
            % Map indices to constellation points
            symbols = constellation(indices + 1);
        end
        
        function [txSymbols, bits] = generateSignal(obj)
            bits = randi([0 1], obj.numBits, 1);
            
            switch obj.modType
                case 'BPSK'
                    txSymbols = 2*bits - 1;
                    
                case 'QPSK'
                    bitsReshaped = reshape(bits, 2, [])';
                    txSymbols = (1/sqrt(2))*(2*bitsReshaped(:,1)-1 + 1j*(2*bitsReshaped(:,2)-1));
                    
                case 'QAM'
                    % Reshape bits into groups for QAM
                    k = log2(obj.M);
                    bitsReshaped = reshape(bits(1:floor(length(bits)/k)*k), k, [])';
                    
                    % Convert binary to decimal
                    symbolIndices = obj.binary2decimal(bitsReshaped);
                    
                    % Map to QAM constellation
                    txSymbols = obj.qamMapping(symbolIndices, obj.M);
            end
        end
        
        % Custom QAM demodulation
        function indices = qamDemapping(~, symbols, M)
            % Calculate number of points per dimension
            L = sqrt(M);
            
            % Generate constellation points
            alphaI = -(L-1):2:L-1;
            alphaQ = -(L-1):2:L-1;
            
            % Create constellation mapping
            [I, Q] = meshgrid(alphaI, alphaQ);
            constellation = (I(:) + 1j*Q(:))/sqrt(10);
            
            % Find nearest constellation point
            [~, indices] = min(abs(repmat(symbols(:), 1, M) - ...
                             repmat(constellation.', length(symbols), 1)), [], 2);
            indices = indices - 1;
        end
        
        % Decimal to binary conversion
        function b = decimal2binary(~, d, width)
            b = zeros(length(d), width);
            for i = 1:width
                b(:,width-i+1) = mod(d, 2);
                d = floor(d/2);
            end
        end
        
        function rxSymbols = applyChannel(obj, txSymbols, SNR)
            SNRLinear = 10^(SNR/10);
            signalPower = mean(abs(txSymbols).^2);
            noisePower = signalPower/SNRLinear;
            
            switch obj.channelType
                case 'AWGN'
                    noise = sqrt(noisePower/2)*(randn(size(txSymbols)) + ...
                            1j*randn(size(txSymbols)));
                    rxSymbols = txSymbols + noise;
                    
                case 'Rayleigh'
                    h = sqrt(1/2)*(randn(size(txSymbols)) + ...
                        1j*randn(size(txSymbols)));
                    noise = sqrt(noisePower/2)*(randn(size(txSymbols)) + ...
                            1j*randn(size(txSymbols)));
                    rxSymbols = h.*txSymbols + noise;
                    
                case 'Rician'
                    K = 10;
                    direct = sqrt(K/(K+1));
                    diffuse = sqrt(1/(2*(K+1)))*(randn(size(txSymbols)) + ...
                              1j*randn(size(txSymbols)));
                    h = direct + diffuse;
                    noise = sqrt(noisePower/2)*(randn(size(txSymbols)) + ...
                            1j*randn(size(txSymbols)));
                    rxSymbols = h.*txSymbols + noise;
            end
        end
        
        function rxBits = demodulate(obj, rxSymbols)
            switch obj.modType
                case 'BPSK'
                    rxBits = real(rxSymbols) > 0;
                    
                case 'QPSK'
                    rxBits = zeros(2*length(rxSymbols), 1);
                    rxBits(1:2:end) = real(rxSymbols) > 0;
                    rxBits(2:2:end) = imag(rxSymbols) > 0;
                    
                case 'QAM'
                    % Demap symbols to indices
                    symbolIndices = obj.qamDemapping(rxSymbols, obj.M);
                    
                    % Convert indices to bits
                    k = log2(obj.M);
                    rxBitsMatrix = obj.decimal2binary(symbolIndices, k);
                    rxBits = rxBitsMatrix(:);
            end
        end
        
        function [ber, theoretical] = analyzeBER(obj)
            ber = zeros(size(obj.SNRdB));
            theoretical = zeros(size(obj.SNRdB));
            
            for i = 1:length(obj.SNRdB)
                [txSymbols, bits] = obj.generateSignal();
                rxSymbols = obj.applyChannel(txSymbols, obj.SNRdB(i));
                rxBits = obj.demodulate(rxSymbols);
                
                % Ensure bits and rxBits are same length for comparison
                minLen = min(length(bits), length(rxBits));
                ber(i) = sum(bits(1:minLen) ~= rxBits(1:minLen))/minLen;
                
                snr = 10^(obj.SNRdB(i)/10);
                switch obj.modType
                    case 'BPSK'
                        theoretical(i) = obj.qfunction(sqrt(2*snr));
                    case 'QPSK'
                        theoretical(i) = obj.qfunction(sqrt(snr));
                    case 'QAM'
                        theoretical(i) = 4*(1-1/sqrt(obj.M))* ...
                            obj.qfunction(sqrt(3*snr/(obj.M-1)))/ ...
                            log2(obj.M);
                end
            end
        end
        
        function plotResults(obj, ber, theoretical)
            figure;
            semilogy(obj.SNRdB, ber, 'o-', 'LineWidth', 2);
            hold on;
            semilogy(obj.SNRdB, theoretical, 'r--', 'LineWidth', 2);
            grid on;
            xlabel('SNR (dB)');
            ylabel('Bit Error Rate (BER)');
            title(['BER Performance Analysis - ' obj.modType ' in ' ...
                   obj.channelType ' Channel']);
            legend('Simulation', 'Theoretical');
            ylim([1e-5 1]);
        end
    end
end