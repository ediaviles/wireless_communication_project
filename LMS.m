function [trained_w] = LMS(trained_w, zk, xk, delta, gamma) % zk is the entire signal, xk is the ideal signal
    % use zk and delta for start of pilot
    N = 1;
    order = length(trained_w);
    received_pilot = zk(delta + 1: delta + length(xk));
    step_size = gamma / (received_pilot * conj(transpose(received_pilot)));
    for j = 1:N
        for i = order+1:length(xk)
            %u = zk(i + delta:-1:i + delta - order + 1); % received pilot
            u = received_pilot(i:-1:i - order + 1);
            v = trained_w * transpose(u); % whatever is the last vk we calculated becomes the next vk for the message
            e = v - xk(i);
            update = step_size * e * conj(u);
            trained_w = trained_w - update;
        end
    end
end

