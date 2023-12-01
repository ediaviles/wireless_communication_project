function [trained_w] = LMS(zk, xk, delta) % zk is the entire signal, xk is the ideal signal
    L1 = -4;
    L2 = 4;
    trained_w = zeros(1, L2 - L1); % w is array of size L2 + abs(L1) + 1 -> Right now 21
    step_size = 0.001; % try normalizing the step size, might help (seems to be ok)
    % use zk and delta for start of pilot
    N = 1;
    order = length(trained_w);
    for j = 1:N
        for i = order+1:length(xk)
            u = zk(i + delta:-1:i + delta - order + 1); % received pilot
            v = trained_w * transpose(u); % whatever is the last vk we calculated becomes the next vk for the message
            e = v - xk(i);
            update = step_size * e * conj(u);
            trained_w = trained_w - reshape(update, 1, order);
        end
    end
end