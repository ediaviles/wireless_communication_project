function [llr, p1] = compute_llr(y, n_0)
    b = 2;
    d = 0.3;
    constellation_points = [
       -1 - j,
       -1 + j,
        1 - j,
        1 + j
    ] .* d;
    n_sym = 2^b;

    bit_sym_map = zeros(n_sym, b);
    for sym_index = 0 : n_sym - 1
        a = dec2bin(sym_index, b);
        for bit_index = 1 :b
            if a(b + 1 - bit_index) == '1'
                bit_sym_map(sym_index + 1, bit_index) = 1;
            end
        end
    end


    p0 = zeros(length(y) * b, 1);
    p1 = zeros(length(y) * b, 1);
    for y_index = 1 : length(y)
        for sym_index = 1  : 2^b
            if length(n_0) == 1
                p_sym = exp(-abs(y(y_index) - constellation_points(sym_index))^2/2/n_0);
            else
                p_sym = exp(-abs(y(y_index) - constellation_points(sym_index))^2/2/n_0(y_index));
            end
            for m_index = 1 : b
                if bit_sym_map(sym_index, m_index) == 0
                    p0((y_index-1) * b +  m_index) = p0((y_index-1) * b +  m_index) + p_sym;
                else
                    p1((y_index-1) * b +  m_index) = p1((y_index-1) * b +  m_index) + p_sym;
                end
            end
        end
    end
    llr = log(p0./p1);
    p1 = p1./(p0 + p1);
end