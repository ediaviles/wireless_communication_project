sigman = 0.2;

%receivedsignal = transmitsignal + sigman/sqrt(2) * (randn(size(transmitsignal))+j*randn(size(transmitsignal)));
%T o test the effect of phase offset and delay, you could simulate such a channel as
%padding = (randn(1,1000) > 0.5) * 2 - 1;
%transmitsignalwithdelay = [zeros(1, 2147), transmitsignal];
%receivedsignal = exp(j*pi/6) * transmitsignalwithdelay + sigman/sqrt(2) * (randn(size(transmitsignalwithdelay))+j*randn(size(transmitsignalwithdelay)));
%receivedsignal = transmitsignal;
matched_filter = flipud(pt);

y = receivedsignal';

t = timing_sync_bits * 0.46;

ps = pilot_sequence * 0.46;

f = fsync_sequence * 0.46;

chunk_size = chunk * 0.46;
chunk_size = length(chunk_size);


%% Filter
y = conv(y, matched_filter);

%% Downsample
tau = mod(length(y), L);
y_synced2 = y(tau +1:end);
y_downsampled = y(1:L:end);

%% Time sync
[corr, lags] = xcorr(y_downsampled, t); % look at class slides which is based on analog signal
[max_value, timing_index] = max(abs(corr));
timing_offset = lags(timing_index);
phase = angle(max_value);
delta = timing_offset + length(t); % determine the start of the first pilot

y_synced = y_downsampled(delta + 1:end);
%y_synced = y_synced * exp(-j*phase);

first_pilot = y_downsampled(delta + 1:delta + length(ps)); % extract first pilot

%% Frame sync
[fcorr, flags] = xcorr(y_downsampled, f);
[~, frame_index] = max(abs(fcorr));
frame_offset = lags(frame_index);

start_first_chunk = frame_offset + length(f);
y_fsynced = y_downsampled(start_first_chunk+1:end);

first_chunk = y_fsynced(1:chunk_size);

%% Equalizer
one_tap = (first_pilot*conj(ps)') / (conj(ps)*ps');
equalizations = [one_tap];

%before_equalization = sample_first_chunk;

%first_chunk = first_chunk; % / one_tap;

delta = chunk_size;

%% Equalize based on number of chunks (n) and remove pilots from message
chunks = [first_chunk];

for i = 1:1:n-1
    % extract nth pilot
    pilot = y_fsynced(delta + 1:delta + length(ps));

    % calculate one tap
    one_tap = (pilot*conj(ps)') / (conj(ps)*ps');
    equalizations = [equalizations, one_tap];

    start_of_chunk = delta + length(ps);
    current_chunk = y_fsynced(start_of_chunk + 1:start_of_chunk + chunk_size);

    % Equalize chunk
    chunk = current_chunk; % / one_tap; % This should be applied after sampling

    % add equalized chunk to message
    chunks = [chunks, chunk];
    
    % update delta to start of next pilot sequence
    delta = start_of_chunk + chunk_size;
end

before_equalizations = chunks;
for i = 1:n
    chunks((i - 1)*chunk_size + 1:i*chunk_size) = chunks((i - 1)*chunk_size + 1:i*chunk_size) / equalizations(i);
end

z_k = chunks;

%% demodulate (threshold)
z_real = real(z_k);
z_demodulated = z_real > 0;


%% BER
message = imread("shannon1440.bmp");
message_vec = reshape(message, 1, []);

bits = message_vec;

BER = mean(z_demodulated ~= bits);
disp(['BER is ', num2str(BER)])


%% Plots

% transmited and received signals and pt
figure(1)
clf
subplot(2,1,1)
plot(t_transmitted, real(transmitsignal),'b')
hold on
plot(t_transmitted, imag(transmitsignal),'r')
legend('real','imag')
ylabel('xI(t)  and  xQ(t)')
xlabel('Time in microseconds')
subplot(2,1,2)
plot(t_received, real(receivedsignal),'b')
hold on
plot(t_received, imag(receivedsignal),'r')
zoom xon
legend('real','imag')
ylabel('yI(t)  and  yQ(t)')
xlabel('Time in microseconds')

figure(2)
clf
subplot(2,1,1)
plot(([0:length(transmitsignal)-1]/length(transmitsignal)-0.5) * Fs / 10^6, abs(fftshift(fft(transmitsignal))))
ylabel('abs(X(f))')
xlabel('Frequency in MHz')
subplot(2,1,2)
plot(([0:length(receivedsignal)-1]/length(receivedsignal)-0.5) * Fs / 10^6, abs(fftshift(fft(receivedsignal))))
ylabel('abs(Y(f))')
xlabel('Frequency in MHz')

figure(3)
clf
subplot(2,1,1);
plot(t_p, pt, 'b');
xlabel('Time in microseconds');
ylabel('Pulse function p(t)');
hold on;
subplot(2,1,2);
plot(([0:length(pt)-1]/length(pt)-0.5) * Fs / 10^6, abs(fftshift(fft(pt))))
xlabel('Frequency in MHz');
ylabel('abs(P(f))');


% recovered image
figure(4)
subplot(2,1,1);
recovered_image = reshape(z_demodulated(1:length(bits)), [45, 32]);
imshow(recovered_image);
subplot(2,1,2);
imshow(message);

% zk - before equalization
figure(5);
plot(real(before_equalizations), imag(before_equalizations), 'rx')

% vk - after equalization
figure(6);
plot(real(z_k(1:1440)), imag(z_k(1:1440)), 'rx');


t_synced = [1:length(y_synced)] / Fs * 10^6;

% y after time synchronization
figure(7)
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

% z demodulated
figure(8)
clf
subplot(2,1,1);
plot(z_demodulated(1:length(bits)), 'r');
subplot(2,1,2);
plot(bits, 'b')



