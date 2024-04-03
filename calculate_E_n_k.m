%{
Title: Octave script for symbol flipping algorithm.
Reference paper: Bing Liu, Wei Tao, "Weighted Symbol-Flipping Decoding for Nonbinary
LDPC Codes", 2010, (https://ieeexplore.ieee.org/abstract/document/5480573)
Description: This script calculates the metric value E_n_k of symbol over GF_0(q) to
determine the flipped symbol, for each symbol and each GF_0(q) element.
GF_0 = vector of Galois field elements except zero
Notation and variables:
E_n_k = flipping metric
N = code length
K = message length
LLR_n_a
H = parity check matrix
rec_vector = received code vector of length N
GF_0 = Vector of Galois field elements except zero
alpha = weighting factor, not exactly defined in the paper, assumed to be 0.5 in this code
Imax = maximum number of iteration
s = syndrome vector
N_m = set of symbol indices connected with mth check node
M_n = set of checks connected with nth symbol-node
w_n_m_a = minimum absolute LLR for the symbols involved in check m except for the nth symbol in rec_vector
%}

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
        N_m{m} = find(H(m,:) != 0); % set of symbol indices connected with mth checknode
        %disp('the set of symbol node index connected to mth checknode =')
        %disp(N_m{m})
    end
    for n = 1:N
        M_n{n} = find(H(:,n) != 0); %set of checks connected with nth symbol-node
    end

    % Calculate the metric E_n_a for each symbol n
    for k = 1:Imax
        for n = 1:N
            E_n_a = zeros(1, N);
            for a = 1:length(GF_0)
                sum_E = 0;
                for m = M_n{n}' % Transpose to ensure a column vector for iteration
                    % symbol indices involved in mth check except nth symbol in rec_vector
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
            % Sum up the metrics for all elements of GF_0(q) to get E_n_k for the current symbol n
            E_n_k(n) = sum(E_n_a);
        end
            %Additional logic to update y and s for the next iteration would go here
    end
end
