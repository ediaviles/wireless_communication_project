function decoded_bits = ldpc_decoder(received_bits, H, max_iterations)
    % received_bits: the received bit vector (likelihoods)
    % H: Parity-check matrix
    % max_iterations: Maximum number of iterations for the decoder

    % Initialize
    [m, n] = size(H); % Size of the parity-check matrix
    msg_from_var_to_chk = repmat(received_bits, m, 1); % Messages from variables to checks
    msg_from_chk_to_var = zeros(m, n); % Messages from checks to variables

    for iteration = 1:max_iterations
        % Check node update
        for i = 1:m
            for j = find(H(i, :))
                product = 1;
                for k = find(H(i, :))
                    if k ~= j
                        product = product * tanh(0.5 * msg_from_var_to_chk(k, i));
                    end
                end
                msg_from_chk_to_var(i, j) = 2 * atanh(product);
            end
        end

        % Variable node update
        for j = 1:n
            for i = find(H(:, j))
                sum = received_bits(j);
                for k = find(H(:, j))
                    if k ~= i
                        sum = sum + msg_from_chk_to_var(k, j);
                    end
                end
                msg_from_var_to_chk(j, i) = sum;
            end
        end

        % Decision
        for j = 1:n
            total_msg = received_bits(j);
            for i = find(H(:, j))
                total_msg = total_msg + msg_from_chk_to_var(i, j);
            end
            decoded_bits(j) = total_msg < 0;
        end

        % Check for convergence
        if all(mod(H * decoded_bits', 2) == 0)
            break;
        end
    end
end
