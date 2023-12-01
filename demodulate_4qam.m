function [signal_bits] = demodulate_4qam(modulated_symbols)
    M = 4;
    b = log2(M);
    constellation = [
       -1 - j,
       -1 + j,
        1 - j,
        1 + j
    ];
    signal_bits = zeros(1, length(modulated_symbols)*b);
    for i = 1:length(modulated_symbols)
        symbol = modulated_symbols(i);
        distances = symbol - constellation;
        [~, min_idx] = min(abs(distances));
        guess = min_idx - 1;
        bit_guess = de2bi(guess, b, "left-msb");
        signal_bits((i-1)*b + 1:(i*b)) = bit_guess;
    end
end