sigman = 0;

receivedsignal = transmitsignal + sigman/sqrt(2) * (randn(size(transmitsignal))+j*randn(size(transmitsignal)));
%T o test the effect of phase offset and delay, you could simulate such a channel as
%padding = (randn(1,1000) > 0.5) * 2 - 1;
transmitsignalwithdelay = [zeros(1, 2147), transmitsignal];
receivedsignal = exp(j*2*pi) * transmitsignalwithdelay + sigman/sqrt(2) * (randn(size(transmitsignalwithdelay))+j*randn(size(transmitsignalwithdelay)));


t_received = [1:length(receivedsignal)] / Fs * 10^6;

matched_filter = fliplr(pt);

y = receivedsignal;

t = timing_sync_bits; 
t_modulated = modulate_4qam(t);
t = t_modulated * d;


ps = pilot_sequence * d;
modulated_pilot = modulate_4qam(pilot_sequence);

f = fsync_sequence;
f_modulated = modulate_4qam(f);
f = f_modulated * d;
f = upsample(f, L);
f = conv(f, pt);
f = f(1:L:end);

chunk_size = length(message_vec)/n/b; % in symbols

%% Time sync
t = upsample(t, L);
t = conv(t, pt);

ideal = t;
[corr_id, lags_id] = xcorr(y, ideal);
[ideal_max_value, ideal_timing_index] = max(abs(corr_id));
ideal_timing_offset = lags_id(ideal_timing_index);

y_sync = y(ideal_timing_offset+1:end); % signal starts at the time sync bits
figure(111);
plot(real(y_sync),'b')
hold on
plot(imag(y_sync),'r')


%% Filter
y_filt = conv(y_sync, matched_filter);

%% DOWNSAMPLE
z_k = y_filt(1:L:end); % this starts at the time symbol seq


%% Start of first pilot
delta = length(t_modulated);

%y_synced = z_k(delta + 1:end); % after time sync
%% Frame sync
%[fcorr, flags] = xcorr(z_k, f);
%[~, frame_index] = max(abs(fcorr));
%frame_offset = flags(frame_index);
frame_offset = length(modulated_pilot) + delta; % I modified this

%% Extract first chunk
start_first_chunk = frame_offset + length(f_modulated); % maybe use different value?
y_fsynced = z_k(start_first_chunk+1:end);

first_chunk = y_fsynced(1:chunk_size);

chunks = [first_chunk']; % get all of our chunks
L2 = 4;
L1 = -4;
gamma = 1;
trained_w = zeros(1, L2 - L1);
trained_w = LMS(trained_w, z_k, modulated_pilot, delta, gamma);
filters = [trained_w']; % train their respective filters
delta = delta + length(modulated_pilot) + length(f_modulated) + chunk_size;
figure(69);
scatter(real(z_k), imag(z_k));
%% Train w for each chunk
for i = 1:n-1
    trained_w = LMS(trained_w, z_k, modulated_pilot, delta, gamma);
    filters = [filters, transpose(trained_w)]; % train their respective filters
    chunks = [chunks, transpose(z_k(delta + length(modulated_pilot) + 1:delta + length(modulated_pilot) + chunk_size))];
    delta = delta + length(modulated_pilot) + chunk_size;
end

before_equalizations = reshape(chunks, 1, []);
zk_equalized = []
%% Equalize
for i = 1:n
    w = transpose(filters(:,i));
    current_chunk = transpose(chunks(:,i));
    vk = filter(w, 1, current_chunk);
    figure(10000 + i);
    scatter(real(current_chunk), imag(current_chunk));
    %scatter(real(vk(L2: chunk_size + L2)), imag(vk(L2: chunk_size + L2)));
    zk_equalized = [zk_equalized, vk];
end


%% Soft decoding (TODO)

%% Demodulate
z_demodulated = demodulate_4qam(zk_equalized);

%% BER
message = imread("shannon1440.bmp");
%message = imread("shannon20520.bmp");

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
plot(real(zk_equalized), imag(zk_equalized), 'rx');

y_synced = y_sync;

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
%figure(8)
%clf
%subplot(2,1,1);
%plot(z_demodulated(1:length(bits)), 'r');
%subplot(2,1,2);
%plot(bits, 'b')

