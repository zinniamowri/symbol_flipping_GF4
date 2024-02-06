if exist('hammingDistance.m', 'file')
    disp('The file is in the current directory.');
else
    disp('The file is not in the current directory.');
end


function distance = hammingDistance(vector1, vector2)
    distance = sum(vector1 != vector2);
end


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

function closestCodeword = findClosestCodeword(receivedVector, validCodewords)
    minDistance = inf; % Initialize minimum distance with infinity
    closestCodeword = receivedVector; % Default to received vector in case no closer codeword is found

    for i = 1:size(validCodewords, 1)
        currentDistance = hammingDistance(receivedVector, validCodewords(i, :));
        if currentDistance < minDistance
            minDistance = currentDistance;
            closestCodeword = validCodewords(i, :);
        end
    end
end


% Received vector
receivedVector = [1 2 3];


disp('received vector')
disp(receivedVector)

% Find the closest valid codeword
closestCodeword = findClosestCodeword(receivedVector, validCodewords);

disp('Closest Codeword:');
disp(closestCodeword);


