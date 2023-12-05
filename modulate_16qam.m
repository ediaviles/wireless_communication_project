function [modulated_symbols] = modulate_16qam(signal_bits)
    M = 16;
    b = log2(M);
    d = sqrt(2)/3/1.64;
    %d = 0.315;
    options = [1.5+1.5j, 0.5+1.5j, -1.5+1.5j, -0.5+1.5j
           1.5+0.5j, 0.5+0.5j, -1.5+0.5j, -0.5+0.5j
           1.5-1.5j, 0.5-1.5j, -1.5-1.5j, -0.5-1.5j
           1.5-0.5j, 0.5-0.5j, -1.5-0.5j, -0.5-0.5j] .* d;
    options = reshape(options, 1, []);

    modulated_symbols = zeros(1, length(signal_bits)/b);
    for i = 1:length(modulated_symbols)
        start_i = b*(i - 1) + 1;
        end_i = b*(i);
        grouping = signal_bits(start_i:end_i);
        
        % Convert binary grouping to decimal
        decimal_value = bi2de(grouping, 'left-msb');
        % QAM modulation
        modulated_symbols(i) = options(decimal_value+1);
    end
end