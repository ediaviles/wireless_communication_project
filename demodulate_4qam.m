function [signal_bits] = demodulate_4qam(modulated_symbols)
    M = 4;
    b = log2(M);
    d = 0.8;
    constellation = [
       -1 - j,
       -1 + j,
        1 - j,
        1 + j
    ] .* d/2;
    signal_bits = zeros(1, length(modulated_symbols)*b);
    for i = 1:length(modulated_symbols)
        symbol = modulated_symbols(i);
        distances = symbol - constellation;
        [~, min_idx] = min(abs(distances));
        guess = min_idx - 1;
        start_i = (i-1)*b + 1;
        end_i = i * b;
        bit_guess = de2bi(guess, b, "left-msb");
        signal_bits(start_i:end_i) = bit_guess;
    end
end