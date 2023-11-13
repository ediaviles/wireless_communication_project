sigman = 0.2;

receivedsignal = transmitsignal + sigman/sqrt(2) * (randn(size(transmitsignal))+j*randn(size(transmitsignal)));
%T o test the effect of phase offset and delay, you could simulate such a channel as
padding = (randn(1,1000) > 0.5) * 2 - 1;
transmitsignalwithdelay = [zeros(1, 2147), transmitsignal, padding];
receivedsignal = exp(j*pi/6) * transmitsignalwithdelay + sigman/sqrt(2) * (randn(size(transmitsignalwithdelay))+j*randn(size(transmitsignalwithdelay)));

y = receivedsignal;
y_symbols = (receivedsignal > 0) * 2 - 1; % turn to symbols
y_symbols = y_symbols * 0.3;

t = timing_sync_bits * 0.3;
t = upsample(t, L);
t = conv(t, fliplr(pt));

ps = pilot_sequence * 0.3;
ps = upsample(ps, L);
ps = conv(ps, fliplr(pt));

f = fsync_sequence * 0.3;
f = upsample(f, L);
f = conv(f, fliplr(pt));


chunk_size = chunk * 0.3;
chunk_size = upsample(chunk_size, L);
chunk_size = conv(chunk_size, fliplr(pt));
chunk_size = length(chunk_size);

% Time sync
[corr, lags] = xcorr(y, t); % look at class slides which is based on analog signal
[~, timing_index] = max(abs(corr));
timing_offset = lags(timing_index);
delta = timing_offset + length(t); % determine the start of the first pilot

first_pilot = y(delta + 1: delta + length(ps)); % extract first pilot

% Frame sync
[fcorr, flags] = xcorr(y, f);
[~, frame_index] = max(abs(fcorr));
frame_offset = lags(frame_index);

start_first_chunk = frame_offset + length(f);

first_chunk = y(start_first_chunk + 1 : start_first_chunk + chunk_size);

filter_first_chunk = conv(first_chunk, fliplr(pt));
sample_first_chunk = filter_first_chunk(1:L:end);

one_tap = (conj(ps)*first_pilot') / (conj(ps)*ps');
sample_first_chunk = sample_first_chunk / one_tap;

delta = start_first_chunk + chunk_size;

% Equalize based on number of chunks (n) and remove pilots from message
equalized_message = [sample_first_chunk];
matched_filter = fliplr(pt);

for i = 1:1:n-1
    % extract nth pilot
    pilot = y(delta + 1:delta + length(ps));

    % calculate one tap
    one_tap = (conj(ps)*pilot') / (conj(ps)*ps');

    start_of_chunk = delta + length(ps);
    current_chunk = y(start_of_chunk + 1:start_of_chunk + chunk_size);

    % filter equalized chunk
    filtered_chunk = conv(current_chunk, matched_filter);

    % sample filtered chunk
    sampled_chunk = filtered_chunk(1:L:end);

    % Equalize chunk
    equalized_chunk = sampled_chunk / one_tap; % This should be applied after sampling

    % add equalized chunk to message
    equalized_message = [equalized_message, equalized_chunk];
    
    % update delta to start of next pilot sequence
    delta = start_of_chunk + chunk_size;
end

z_k = equalized_message;

% demodulate (threshold)
z_demodulated = z_k > 0;


% BER
message = imread("shannon1440.bmp");
message_vec = reshape(message, 1, []);

bits = message_vec;
BER = mean(z_demodulated(1:length(bits)) ~= bits);
disp(['BER is ', num2str(BER)])

figure(3);
plot(real(z_k(1:1440)), imag(z_k(1:1440)), 'rx');

figure(1);
plot(real(receivedsignal), imag(receivedsignal), 'rx');

