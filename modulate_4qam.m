function [modulated_symbols] = modulate_4qam(signal_bits)
    M = 4;
    b = log2(M);
    constellation = [
       -1 - j,
       -1 + j,
        1 - j,
        1 + j
    ];
    modulated_symbols = zeros(1, length(signal_bits)/b);
    for i = 1:length(modulated_symbols)
        bits = signal_bits((i - 1)*b + 1:i*b); % bits of length b
        bits
        decimal = bi2de(bits, b, 'left-msb');
        mapping = constellation(decimal + 1);
        modulated_symbols(i) = mapping;
    end
end