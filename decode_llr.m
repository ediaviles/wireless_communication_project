function decoded_bits = decode_llr(received_llrs, H, max_iterations)
    % received_llrs: the received LLRs (log-likelihood ratios)
    % H: Parity-check matrix
    % max_iterations: Maximum number of iterations for the decoder

    % Initialize 
    [m, n] = size(H); % Size of the parity-check matrix
    length(received_llrs)
    n
    if length(received_llrs) ~= n
        error('Length of received_llrs must match the number of columns in H');
    end
    msg_from_var_to_chk = repmat(received_llrs, m, 1); % Messages from variables to checks
    msg_from_chk_to_var = zeros(m, n); % Messages from checks to variables
    decoded_llrs = zeros(1, n); % Initialize decoded LLRs

    for iteration = 1:max_iterations
        % Check node update (using min-sum approximation for simplicity)
        for i = 1:m
            for j = find(H(i, :))
                temp = msg_from_var_to_chk(:, i);
                temp(j) = 0; % Exclude the current variable node
                msg_from_chk_to_var(i, j) = prod(sign(temp)) * min(abs(temp(temp ~= 0)));
            end
        end

        % Variable node update
        for j = 1:n
            for i = find(H(:, j))'
                msg_from_var_to_chk(j, i) = received_llrs(j) + sum(msg_from_chk_to_var(setdiff(find(H(:, j)), i), j));
            end
        end

        % Compute total LLRs for decision
        for j = 1:n
            decoded_llrs(j) = received_llrs(j) + sum(msg_from_chk_to_var(find(H(:, j)), j));
        end

        % Make a soft decision
        decoded_bits = decoded_llrs < 0;

        % Check for convergence
        if all(mod(H * double(decoded_bits)', 2) == 0)
            break;
        end
    end
end
