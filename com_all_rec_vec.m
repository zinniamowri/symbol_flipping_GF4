if exist('com_all_rec_vec.m', 'file')
    disp('The file is in the current directory.');
else
    disp('The file is not in the current directory.');
end


filename = 'received_vectors_and_syndromes.txt';

f_name = fopen(filename, 'w');

%Write a header
fprintf(f_name, 'v1 v2 v3 s1 s2\n');

gf4_elements = [0, 1, 2, 3];
H= [1,0,3;0,1,2];

% Generate all possible combinations for vectors of length 3
[A, B, C] = ndgrid(gf4_elements, gf4_elements, gf4_elements);

% Convert the ndgrid output to a list of vectors
all_combinations = [A(:), B(:), C(:)];

% Initialize a matrix to store the syndromes
syndromes = zeros(size(all_combinations, 1), size(H, 1));

% Loop through all combinations to compute syndromes
for i = 1:size(all_combinations, 1)
    % Extract the current received vector
    r = all_combinations(i, :);

    % Compute the syndrome for the current vector
    H_t = H';
    syndrome = mod(r * H_t, 4);

    if all(syndrome==0)
      disp('valid code words');
      disp(r);
    endif

    % Store the syndrome
    syndromes(i, :) = syndrome;


    %combine received vector and corresponding syndrome
    combined_data = [all_combinations, syndromes];


end

% Write the data
for i = 1:size(combined_data, 1)
    fprintf(f_name,'%d %d %d %d %d\n', combined_data(i, :));
end

% Close the file
fclose(f_name);




