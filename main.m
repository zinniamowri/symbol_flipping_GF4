%{
Title: Octave script for symbol flipping algorithm.

Reference paper: Bing Liu, Wei Tao, "Weighted Symbol-Flipping Decoding for Nonbinary
LDPC Codes", 2010, (https://ieeexplore.ieee.org/abstract/document/5480573)

Description: This script calculates the LLR value of each symbol over GF_0(q) and calls the
funciton calculate_E_n_k which calculates the flipping metric E_n_k for symbol over GF_0(q)

Notation and variables:

N = code length
K = message length
b = number of bits for each symbol in GF(q) field
q = number of GF field elements
t = bipolar vector after mapping to bpsksymbol {-1,+1}
r = channel output after adding noise n; r=t+n
LLR_n_a = vector of log-likelihood ratio of nonzero element a over GF(q)
H = prity check matrix
rec_vector = received code vector of length N
GF_0 = Vector of Galois feild elements except zero
alpha = weighting factor, not exactly defined in the paper, assumed to be 0.5 in this code
Imax = maximum number of iteration
r_idx = indices of r vector for corresponding nth symbol in the received vector

%}

N = 3;    % Length of the codeword or received vector
K = 2;    % Length of the recieved messege
b = 2;    % Number of bits per symbol (q = 2^b)
Eb_No = 10; % SNR
q=4; % number of GF(4) elements
GF_0 = [1,2,3]; %GF(4) elements except zero

bitVectors = []; % contains the bit representation of each symbol in the GF field
LLR_n_a = [];

  for i = 1:(q-1)
    binaryString = dec2bin(GF_0(i), b);
    numericVector = binaryString - '0';
    bitVectors = [bitVectors, numericVector];
  end

disp(bitVectors)

% BPSK Modulation: Map 0 to -1, and 1 to 1
t = 2 * bitVectors - 1; % t = bpsksymbol {-1, 1}

disp('bipolar vector t after BPSK modulation')
disp(t)

% Add Gaussian noise
sigma2 = 1 / (2 * K/N * Eb_No); % Noise variance
n = sqrt(sigma2) * randn(1, length(t)); % AWGN, noise
r = t + n; % channel output after adding noise

disp('channel output after adding noise, r = t + n')
disp(r)

  % Iterate over GF_0 and map each symbol to its corresponding channel output r
  for i = 1:length(GF_0)
      sum=0;
      r_idx = (i-1)*b + (1:b);
      %disp(r_indices)
      for j=1:b
        if bitVectors(r_idx(j))==1
          sum=sum+r(r_idx(j));
        end
      end
      LLR_n_a(i) = sum;
  end

disp('LLRs for each symbol L_n_a is ')
disp(LLR_n_a)

H = [1, 0, 2; 0, 1, 3];
M = rows(H);
rec_vector = [1, 2, 2];
alpha = 0.5;  % Weighting factor
Imax = 5;

% A higher absolute value indicates a higher confidence that flipping
% this symbol could lead to a valid codeword.
E_n_k = calculate_E_n_k(N, M, LLR_n_a, H, rec_vector, GF_0, alpha, Imax);

disp('E_n_k:');
disp(E_n_k);

[E_max, symbol_to_flip] = max(E_n_k);

% Display the index of the symbol to flip
disp('Symbol to flip');
disp(symbol_to_flip);
