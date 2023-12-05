function [demodulated_signal] = demodulate_16qam(symbols)
    M = 16;
    b = log2(M);
    d = sqrt(2)/3/1.64;
    options = [1.5+1.5j, 0.5+1.5j, -1.5+1.5j, -0.5+1.5j
           1.5+0.5j, 0.5+0.5j, -1.5+0.5j, -0.5+0.5j
           1.5-1.5j, 0.5-1.5j, -1.5-1.5j, -0.5-1.5j
           1.5-0.5j, 0.5-0.5j, -1.5-0.5j, -0.5-0.5j] .* d;
    options = reshape(options, 1, []);
    demodulated_signal = zeros(1, length(symbols)*b);
    for i = 1:length(symbols)
        symbol = symbols(i);
        start_i = b*(i - 1) + 1;
        end_i = b*(i);
        

        [min_val, min_idx] = min(symbol - options);
        decimal = min_idx - 1;
        
        % Convert binary grouping to decimal
        binary_value = de2bi(decimal, b, 'left-msb');
        % QAM modulation
        demodulated_signal(start_i:end_i) = binary_value;
    end
end