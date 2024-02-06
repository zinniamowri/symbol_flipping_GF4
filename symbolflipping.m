if exist('symbolflipping.m', 'file')
    disp('The file is in the current directory.');
else
    disp('The file is not in the current directory.');
end


function output = symbolflipping(r,H)

    % Initialize symbol set
    gf4_elements = [0, 1, 2, 3];

    % Check if received vector is a valid codeword
    s = calculateSyndrome(r, H);

    if all(s==0)
      output=r;
      return
    endif

    % Attempt to correct the vector by flipping each symbol
    for i = 1:length(r)
        originalSymbol = r(i);
        for j=1:length(gf4_elements)
            symbol = gf4_elements(j)
                % Skip the original symbol
              if symbol == originalSymbol
                continue;
            endif

            % replace the symbol r(i) with GF elements
            r(i) = symbol;

            % Check if the syndrome is now zero

            s = calculateSyndrome(r, H);

            if all(s == 0)
                output = r;
                return
            endif

            % Revert to the original symbol before trying the next position
            r(i) = originalSymbol;
         endfor

           % If no correction was found, return the original vector

   endfor

   %output = r;
   disp('correct vector not found')

end


% parity check and syndrome claculation

function s = calculateSyndrome(r,H)

    H_t=H';
    s = mod(r * H_t, 4);

end

rec = [0, 0, 1];

disp('received vector')
disp (rec)

pcm= [1,0,3;0,1,2];

result = symbolflipping(rec,pcm);

disp('correct vector is')

disp(result)

