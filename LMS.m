function [trained_w] = LMS(zk, xk, delta) % zk is the entire signal, xk is the ideal signal
    L1 = -8;
    L2 = 8;
    trained_w = zeros(1, L2 - L1 + 1); % w is array of size L2 + abs(L1) + 1 -> Right now 21
    v = zeros(1, length(xk));
    e = zeros(1, length(xk));
    step_size = 5;
    % use zk and delta for start of pilot
    new_w = zeros(1, length(trained_w));
    for i = 1:length(xk)
        u = zk(i + delta + abs(L1):-1:i + delta - L2);
        v(i) = trained_w * u';
        e(i) = v(i) - xk(i);
        new_w =      - step_size * e(i) * conj(u);
        trained_w = new_w;
    end
end