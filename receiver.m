% sampling freq.:
Fs = 200 * 10^6;

% symbol period:
T = 1/20 * 10^-6;

% oversampling factor:
L = floor(Fs*T);

N = 51; % Length of filter in symbol periods.
Ns = floor(N*L); % Number of filter samples

y = receivedsignal;
y_symbols = (receivedsignal > 0) * 2 - 1; % turn to symbols
y_symbols = y_symbols * 0.3;

t = timing_sync_bits * 0.3;
t = upsample(t, L);
t = conv(t, fliplr(pt));

ps = pilot_sequence * 0.3;
ps = upsample(ps, L);
ps = conv(ps, fliplr(pt));

% Time sync
[corr, lags] = xcorr(y_symbols, t);
[~, timing_index] = max(abs(corr));
timing_offset = lags(timing_index);
delta = timing_offset + length(t); % determine the start of the first pilot

% Equalize based on number of chunks (n) and remove pilots from message
message = [];
for i = i:n
    % extract nth pilot
    pilot = y(delta + 1:delta + length(ps));

    % calculate one tap
    one_tap = (conj(ps)*pilot) / (conj(ps)*ps');

    % equalize only that chunk
    start_of_chunk = delta + length(ps);
    equalized_chunk = y(start_of_chunk + 1:start_of_chunk + chunk_size) / one_tap;

    % add equalized chunk to message
    message = [message, equalized_chunk];

    % update delta to start of next pilot sequence
    delta = start_of_chunk + chunk_size;
end

y = message;

y_I = real(y);
y_Q = imag(y);

% filter using matched filter
matched_filter = fliplr(pt);
z_I = conv(y_I, matched_filter);
z_Q = conv(y_Q, matched_filter);

% Sample
z_Ik = z_I(1:L:end);
z_Qk = z_Q(1:L:end);

z_k = z_Ik + j * z_Qk;

% demodulate (threshold)
z_demodulated = z_k > 0;


% BER
message = imread("shannon1440.bmp");
message_vec = reshape(message, 1, []);

bits = message_vec';
BER = mean(z_demodulated(1:length(bits)) ~= bits);
disp(['BER is ', num2str(BER)])

% Plot constellation
figure(1)
stem([1:length(bits)], message_vec,'bx')
hold on
stem([1:length(bits)], z_demodulated(1:length(bits)),'ro')
legend('x_k', 'z_k')
xlabel('discrete time k')
axis tight