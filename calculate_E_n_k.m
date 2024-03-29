if exist('calculate_E_n_k.m', 'file')
    disp('The file is in the current directory.');
else
    disp('The file is not in the current directory.');
end


function E_n_k = calculate_E_n_k(N, M, LLR_n_a, H, rec_vector, GF_0, alpha, Imax)
    % Number of symbols
    %N = length(rec_vector);
    %M = rows(H)
    s = mod((rec_vector*H'),4);
    disp('the syndrome vector s= ')
    disp(s);
    E_n_k = zeros(N, 1); % E_n_k is a vector with N elements
    % Determine N(m) and M(n) for each check m and symbol n
    N_m = cell(M, 1);
    M_n = cell(N, 1);
    symbol_set = cell (1,N);

    for m = 1:M
        N_m{m} = find(H(m,:) != 0); % set of symbols connected with mth checknode
        %disp('the set of symbol node index connected to mth checknode =')
        %disp(N_m{m})
    end
    for n = 1:N
        M_n{n} = find(H(:,n) != 0); %set of checks connected with nth symbolnode
    end

    % Calculate the metric E_n_a for each symbol n
    for k = 1:Imax
        for n = 1:N
            E_n_a = zeros(1, N);
            for a = 1:length(GF_0)
                sum_E = 0;
                for m = M_n{n}' % Transpose to ensure a column vector for iteration
                    % symbol index involved in mth check except nth symbol in rec_vector
                    sym_in_check_idx = setdiff(N_m{m}, n);
                    symbol_set = H(m,sym_in_check_idx);
                    [~, gf0_idx] = ismember(symbol_set, GF_0);
                    LLR_set = LLR_n_a(gf0_idx); % LLR for symbols except nth symbol in rec_vector
                    % w_n_m_a, the minimum absolute LLR for the symbols in check m
                    w_n_m_a = min(abs(LLR_set));
                    if s(m)==0
                      sum_E = sum_E - w_n_m_a; % S(m) is element from the Syndrome vector
                    else
                      sum_E = sum_E + w_n_m_a;
                    end

                end
                E_n_a(a) = sum_E - alpha * abs(LLR_n_a(a));
            end
            % Sum up the metrics for all elements of GF(q) to get E_n_k for the current symbol n
            E_n_k(n) = sum(E_n_a);
        end
            %disp(E_n_k);
            %Additional logic to update y and s for the next iteration would go here
    end
end

%calculation of LLR (LLR_n_a) for each symbol in GF field
%Excep the zero element

N = 3;    % Length of the codeword
K = 2;    % Length of the recieved messege
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
