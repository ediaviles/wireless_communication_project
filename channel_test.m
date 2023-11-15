sigman = 0.2;

receivedsignal = transmitsignal + sigman/sqrt(2) * (randn(size(transmitsignal))+j*randn(size(transmitsignal)));
%T o test the effect of phase offset and delay, you could simulate such a channel as
padding = (randn(1,1000) > 0.5) * 2 - 1;
transmitsignalwithdelay = [zeros(1, 2147), transmitsignal, padding];
receivedsignal = exp(j*pi/6) * transmitsignalwithdelay + sigman/sqrt(2) * (randn(size(transmitsignalwithdelay))+j*randn(size(transmitsignalwithdelay)));

y = receivedsignal;

t = timing_sync_bits * 0.45;
%t = conj(t);

ps = pilot_sequence * 0.45;
%ps = conj(ps);

f = fsync_sequence * 0.45;
%f = conj(f);

chunk_size = chunk * 0.45;
chunk_size = length(chunk_size);

y_symbols = (y > 0) * 2 - 1;
y_symbols = y_symbols * 0.45;

%% Matched Filter
matched_filter = flipud(pt); % should be something else look at homework 3 and hw 4
y = conv(matched_filter, y);
y = y/L; % Divide by oversampling factor

y_downsampled = y(Ns:L:end);

%% Time sync
[corr, lags] = xcorr(y_downsampled, t); % look at class slides which is based on analog signal
[max_value, timing_index] = max(abs(corr));
timing_offset = lags(timing_index); % seems good
phase = angle(max_value);
delta = timing_offset + length(t); % determine the start of the first pilot

y_synced = y_downsampled(delta + 1:end);
%y_synced = y_synced * exp(-j*phase);

first_pilot = y_synced(1:length(ps)); % extract first pilot

%% Frame sync
[fcorr, flags] = xcorr(y_synced, f);
[~, frame_index] = max(abs(fcorr));
frame_offset = lags(frame_index);

start_first_chunk = frame_offset + length(f);
y_fsynced = y_synced(start_first_chunk + 1:end);

first_chunk = y_fsynced(1:chunk_size);

%% Equalizer
one_tap = (conj(ps)*first_pilot') / (conj(ps)*ps');

equalized_chunk = first_chunk / one_tap;

delta = start_first_chunk + chunk_size;

%% Equalize based on number of chunks (n) and remove pilots from message
equalized_message = [equalized_chunk];
matched_filter = fliplr(pt);

for i = 1:1:n-1
    % extract nth pilot
    pilot = y(delta + 1:delta + length(ps));

    % calculate one tap
    one_tap = (conj(ps)*first_pilot') / (conj(ps)*ps');

    start_of_chunk = delta + length(ps);
    current_chunk = y_fsynced(start_of_chunk + 1:start_of_chunk + chunk_size);

    % Equalize chunk
    equalized_chunk = current_chunk / one_tap; % This should be applied after sampling

    % add equalized chunk to message
    equalized_message = [equalized_message, equalized_chunk];
    
    % update delta to start of next pilot sequence
    delta = start_of_chunk + chunk_size;
end

z_k = equalized_message;

%% demodulate (threshold)
z_real = real(z_k);
z_demodulated = z_real > 0;


% BER
message = imread("shannon1440.bmp");
message_vec = reshape(message, 1, []);

bits = message_vec;

BER = mean(z_demodulated ~= bits);
disp(['BER is ', num2str(BER)])

figure(7)
subplot(2,1,1);
recovered_image = reshape(z_demodulated(1:length(bits)), [45, 32]);
imshow(recovered_image);
subplot(2,1,2);
imshow(message);

figure(11);
plot(real(z_k(1:1440)), imag(z_k(1:1440)), 'rx');

figure(12);
plot(real(receivedsignal), imag(receivedsignal));

t_synced = [1:length(y_synced)] / Fs * 10^6;

figure(13)
clf
subplot(2,1,1)
plot(t_synced, real(y_synced),'b')
hold on
plot(t_synced, imag(y_synced),'r')
legend('real','imag')
ylabel('yI(t)  and  yQ(t)')
xlabel('Time in microseconds')
subplot(2,1,2);
plot(([0:length(y_synced)-1]/length(y_synced)-0.5) * Fs / 10^6, abs(fftshift(fft(y_synced))))
xlabel('Frequency in MHz');
ylabel('abs(P(f))');

figure(14)
clf
subplot(2,1,1);
plot(z_demodulated(1:length(bits)), 'r');
subplot(2,1,2);
plot(bits, 'b')