# Project Title
Experimenting Symbol Flipping Algorithm for LDPC code over $GF(4)$ field elements
- separate scripts are created for now to check whether each functionality is working or not
- while scaling up the field size and generalizing the algorithms, combining the separate scripts might help 

# Description

## com_all_rec_vec.m
- This script generates all possible received vector of a certain length according to the given parity check matrix.
- Calculates the syndrome for all possible received vector
- writes the received vector and its corresponding syndrome value in a text file

## validCodewords.m
- This script gives all the valid codewords by calculating syndrome for a given parity check matrix.
- This has been used to check the functionality
- same logic has been used in hammingDistance.m to calculate the valid codewords and to calculate the closest valid codewords to the received vectors

## symbolflipping.m

- This scripts first defines elements of GF field. For now, everything is hardcoded input. So for $GF(4)$ the elements are given $[0,1,2,3]$
- Also a parity check matrix is given to check the syndrome of a given received vector
- if the syndrome is not zero, it iterates through all the symbols of GF field and replace
- finally, it finds a valid codeword for the received erroneous vector and returns the corrected vector.
- this script is working in there is one position of error only
- there are some limitations as well, we need to locate the error position in some other way to avoid convergence to other codewords instead of the correct one.

## hammingDistance.m

- this script generates all the valid codewords for a given parity check matrix
- for a received vector it checks the hamming distance with all the valid codewords
- it returns the valid codeword which has a minimum hamming distance with the received codeword
- again, in this case, error location detection logic has to be implemented in a more precise way to avoid convergence to another valid codeword instead of the correct one.


