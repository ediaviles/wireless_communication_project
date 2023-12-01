



function [codeword] = encode_bits(bits, K, N, H)
    M = N - K;
    max_ones_row = 10;
    max_ones_col = 10;

    codeword = zeros(1, N);
    codeword(1:K) = bits;

    parity = zeros(1, M);
    
    for i = 1:M
        for j = 1:max_ones_row
            if ()
        end
    end
end