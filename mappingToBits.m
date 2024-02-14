if exist('mappingToBits.m', 'file')
    disp('The file is in the current directory.');
else
    disp('The file is not in the current directory.');
end

% Convert each codeword symbol into its binary representation

function bitVectors = mappingToBits(codeword, b)

  bitVectors = zeros(length(codeword), b);

  for i = 1:length(codeword)
    binaryString = dec2bin(codeword(i), b);
    bitVectors(i, :) = arrayfun(@(x) str2double(x), binaryString);
  endfor

end

function gfSymbols = bitsToGF4(demodulatedBits)
  % Ensure demodulatedBits is a row vector for the processing
  if iscolumn(demodulatedBits)
    demodulatedBits = demodulatedBits';
  end

  % Calculate the number of GF(4) symbols
  numSymbols = length(demodulatedBits) / 2;

  % Reshape the bit vector to have 2 rows, each row represents a bit of the GF(4) symbol
  bitMatrix = reshape(demodulatedBits, 2, numSymbols);

  % Convert binary to decimal: bitMatrix(1, :) is the MSB, bitMatrix(2, :) is the LSB
  gfSymbols = bitMatrix(1, :) * 2 + bitMatrix(2, :);
end


function y_hard_limited = bpsk_ldpc_decode(codeword, Eb_No, N, K, b)

    bitVectors = mappingToBits(codeword, b);

    % Flatten bitVectors to a single row for BPSK modulation
    bpskInput = reshape(bitVectors.', 1, []);
    disp('mapped codeword to bits')
    disp(bpskInput);

    % BPSK Modulation: Map 0 to -1, and 1 to 1
    t = 2 * bpskInput - 1; % t = bpsksymbol -1 and 1

    disp('The received codeword is now mapped to a bipolar vector t')
    disp(t)

    % Add Gaussian noise
    sigma2 = 1 / (2 * K/N * Eb_No); % Noise variance
    n = sqrt(sigma2) * randn(1, length(t)); % AWGN
    r = t + n; % Received signal after adding noise

    disp('Received signal r = t + n')
    disp(r);

    % BPSK Demodulation and Hard Decision
    % Convert received values to hard bits (0 or 1)
    demodulatedBits = r >= 0; % 1 if r >= 0, else 0

    disp('demodulatedBits')
    disp(demodulatedBits)

    % Convert demodulated bits back to GF(4) symbols
    gfSymbols = bitsToGF4(demodulatedBits);
    y_hard_limited=gfSymbols;

end

% Define parameters
N = 6;    % Length of the codeword
K = 3;    % Dimension of the codeword
b = 2;    % Number of bits per symbol (q = 2^b)
Eb_No = 10; % Energy per bit to noise power spectral density ratio SNR

codeword = [0, 1, 2];

disp('received codeword')
disp(codeword)


y_hard_limited = bpsk_ldpc_decode(codeword, Eb_No, N, K, b);

disp('Hard-limited symbol decision y is:');
disp(y_hard_limited);

