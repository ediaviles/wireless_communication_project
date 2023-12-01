function [trained_w] = LMS(trained_w, zk, xk, delta, step_size) % zk is the entire signal, xk is the ideal signal
    % use zk and delta for start of pilot
    N = 10;
    order = length(trained_w);
    for j = 1:N
        for i = order+1:length(xk)
            u = zk(i + delta:-1:i + delta - order + 1); % received pilot
            v = trained_w * transpose(u); % whatever is the last vk we calculated becomes the next vk for the message
            e = v - xk(i);
            update = step_size * e * conj(u);
            trained_w = trained_w - update;
        end
    end
end

