
gf4_elements = [0, 1, 2, 3];
H= [1,0,3;0,1,2];

% Generate all possible combinations for vectors of length 3
[A, B, C] = ndgrid(gf4_elements, gf4_elements, gf4_elements);

% Convert the ndgrid output to a list of vectors
all_combinations = [A(:), B(:), C(:)];

validCodewords = [];

% Initialize a matrix to store the syndromes
syndromes = zeros(size(all_combinations, 1), size(H, 1));

% Loop through all combinations to compute syndromes
for i = 1:size(all_combinations, 1)
    % Extract the current vector
    r = all_combinations(i, :);

    % Compute the syndrome for the current vector
    H_t = H';
    syndrome = mod(r * H_t, 4);

    if all(syndrome==0)
      validCodewords = [validCodewords;r];
    endif

end

disp(validCodewords);




