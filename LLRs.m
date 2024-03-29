%calculation of LLR (LLR_n_a) for each symbol in GF field
%Excep the zero element

N = 6;    % Length of the codeword
K = 3;    % Length of the recieved messege
b = 2;    % Number of bits per symbol (q = 2^b)
Eb_No = 10; % SNR
q=4; % number of GF(4) elements

GF_0 = [1,2,3]; %GF(4) elements except zero

bitVectors = [];
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
n = sqrt(sigma2) * randn(1, length(t)); % AWGN
r = t + n; % Received signal after adding noise

disp('channel output after adding noise, r = t + n')
disp(r)

% Iterate over GF_0 and map each symbol to its corresponding channel output r
for i = 1:length(GF_0)
    sum=0;
    r_indices = (i-1)*b + (1:b);
    %disp(r_indices)
    for j=1:b
      if bitVectors(r_indices(j))==1
        sum=sum+r(r_indices(j));
      end
    end
    LLR_n_a(i) = sum;
end

disp('LLRs for each symbol L_n_a is ')
disp(LLR_n_a)

