
% Convert each codeword symbol into its binary representation

function bitVectors = mappingToBits(received_vector, b)

  bitVectors = zeros(length(received_vector), b);

  for i = 1:length(received_vector)
    binaryString = dec2bin(received_vector(i), b);
    bitVectors(i, :) = arrayfun(@(x) str2double(x), binaryString);
  endfor

end

function gfSymbols = bitsToGF4(bit_vector,b)
  % if bit_vector is a column vector transpose it to roe vector
  if iscolumn(bit_vector)
    bit_vector = bit_vector';
  end

  % Calculate the number of GF(4) symbols
  numSymbols = length(bit_vector) / b;

  % reshape the bit vector to a matrix
  % m=2 rows, each row represents a bit of the GF(4) symbols, column=numSymbols
  % first row contains the first bit of each symbol, and the second row contains the second bit of each symbol
  bitMatrix = reshape(bit_vector, b, numSymbols);

  % convert binary to decimal: bitMatrix(1, :) is the MSB, bitMatrix(2, :) is the LSB
  % gfSymbols = bitMatrix(1, :) * 2 + bitMatrix(2, :);
  % calculate weights for each binary position in descending order
  binaryWeights = 2.^(b-1:-1:0);

  % Convert binary to decimal for each symbol
  gfSymbols = (binaryWeights* bitMatrix);
end


function r = bpsk_output(received_vector, Eb_No, N, K, b)

    bitVectors = mappingToBits(received_vector, b);

    % convert bitVectors to a single row for BPSK modulation
    bpskInput = reshape(bitVectors.', 1, []);
    disp('mapped received_vector to bits')
    disp(bpskInput);

    % BPSK Modulation: Map 0 to -1, and 1 to 1
    t = 2 * bpskInput - 1; % t = bpsksymbol -1 and 1

    disp('The received codeword is now mapped to a bipolar vector t')
    disp(t)

    % Add Gaussian noise
    sigma2 = 1 / (2 * K/N * Eb_No); % Noise variance
    n = sqrt(sigma2) * randn(1, length(t)); % AWGN
    r = t + n; % Received signal after adding noise

end

N = 6;    % Length of the codeword
K = 3;    % Length of the recieved messege
b = 2;    % Number of bits per symbol (q = 2^b)
Eb_No = 10; % SNR


received_vector = [2, 1, 2];

disp('received vector')
disp(received_vector)

disp('channel output r = t + n')
r = bpsk_output(received_vector, Eb_No, N, K, b);

% convert channel output r to hard decision bit vector x

x = r >= 0; % 1 if r >= 0, else 0, x_i=sign(r_i),  sign(r)=1 if r>=1 else 0

disp('hard decision bit vector x derived from r')
disp(x)

%Convert x to GF(4) symbols
% the hard-limited symbol decision vector y
y_hard_decision = bitsToGF4(x,b);

disp('Hard-limited symbol decision y is:');
disp(y_hard_decision);
