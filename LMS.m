function [trained_w] = LMS(trained_w, zk, xk, delta, gamma) % zk is the entire signal, xk is the ideal signal
    % use zk and delta for start of pilot
    N = 5;
    filter_size = length(trained_w);
    step_size = gamma / (zk * conj(transpose(zk)));
    %step_size = 0.1;
    for j = 1:N
        for i = filter_size:length(xk)
            %u = zk(i + delta:-1:i + delta - order + 1); % received pilot
            u = zk(i:-1:i - filter_size + 1);
            v = trained_w * transpose(u); % whatever is the last vk we calculated becomes the next vk for the message
            e = v - xk(i);
            update = step_size * e * conj(u);
            trained_w = trained_w - update;
        end
    end
end