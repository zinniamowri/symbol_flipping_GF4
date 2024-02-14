if exist('calculate_E_n_k.m', 'file')
    disp('The file is in the current directory.');
else
    disp('The file is not in the current directory.');
end


function E_n_k = calculate_E_n_k(s, L_n_a, H, rec_codeword, alpha, Imax)
    % Number of symbols
    N = length(rec_codeword);
    % Number of non-zero GF(q) elements (assuming GF(q) = GF(4) in this example)
    q = 4 - 1;

    E_n_k = zeros(N, 1); % E_n_k is a vector with N elements

    % Determine N(m) and M(n) for each check m and symbol n
    N_m = cell(rows(H), 1);
    M_n = cell(N, 1);
    for m = 1:rows(H)
        N_m{m} = find(H(m,:) != 0);
    end
    for n = 1:N
        M_n{n} = find(H(:,n) != 0);
    end

    % Calculate the metric E_n_a for each symbol n and each GF(q) element a
    for k = 1:Imax
        for n = 1:N
            E_n_a = zeros(1, q); % Temporary variable for the current symbol
            for a = 1:q
                % Initialize the sum for the current n and a
                sum_E = 0;
                for m = M_n{n}' % Transpose to ensure a column vector for iteration
                    % Find the symbols involved in the current check m
                    symbols_in_check = setdiff(N_m{m}, n);
                    % Calculate w_n_m_a, the minimum absolute LLR for the symbols in check m
                    w_n_m_a = min(abs(L_n_a(symbols_in_check, a)));
                    sum_E = sum_E + (2 * s(m) - 1) * w_n_m_a;
                end
                E_n_a(a) = sum_E - alpha * L_n_a(n, a);
            end
            % Sum up the metrics for all elements of GF(q) to get E_n_k for the current symbol n
            E_n_k(n) = sum(E_n_a);
        end
        % Here, E_n_k contains the sum E_n_a for all symbols for the current iteration k
        % Additional logic to update y and s for the next iteration would go here
        % ...
    end
end

s = [0, 1];  % Syndrome vector for the current iteration , need to replace with a funciton
L_n_a = [2, 2, 3; 2, 2, 3; 2, 2, 3];  % LLR values for each symbol and each GF(q) element, need to replace with a funciton
H = [1 0 1; 0 1 1];  % Parity-check matrix
rec_codeword = [2, 2, 3];  % Received symbol vector
alpha = 0.5;  % Weighting factor
Imax = 5;


E_n_k = calculate_E_n_k(s, L_n_a, H, rec_codeword, alpha, Imax);


disp('E_n_k:');
disp(E_n_k);


[E_max, symbol_to_flip] = max(E_n_k);

% Display the index of the symbol to flip
disp('Symbol to flip');
disp(symbol_to_flip);





