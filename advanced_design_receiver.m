sigman = 0.2;

receivedsignal = transmitsignal + sigman/sqrt(2) * (randn(size(transmitsignal))+j*randn(size(transmitsignal)));
%T o test the effect of phase offset and delay, you could simulate such a channel as
%padding = (randn(1,1000) > 0.5) * 2 - 1;
transmitsignalwithdelay = [zeros(1, 2147), transmitsignal];
receivedsignal = exp(j*2*pi) * transmitsignalwithdelay + sigman/sqrt(2) * (randn(size(transmitsignalwithdelay))+j*randn(size(transmitsignalwithdelay)));
matched_filter = fliplr(pt);

r = transmitsignal;
%debug
g = conv(r, matched_filter);
g = g(1:L:end);
g = g(201:end);
g = g(1:20520/2);

g = demodulate_4qam(g);

figure(58)
subplot(2,1,1);
recovered_image = reshape(g, [171, 120]);
imshow(recovered_image);
subplot(2,1,2);
imshow(message);

t_received = [1:length(receivedsignal)] / Fs * 10^6;


y = receivedsignal;

t = timing_sync_bits; 
t_modulated = modulate_4qam(t);
t = t_modulated * d;


ps = pilot_sequence * d;
modulated_pilot = modulate_4qam(pilot_sequence);
p = modulated_pilot * d;
p = upsample(p, L);
p = conv(p, pt);

f = fsync_sequence;
f_modulated = modulate_4qam(f);
f = f_modulated * d;
f = upsample(f, L);
f = conv(f, pt);

%% Filter
zt = conv(y, matched_filter);

%% Time sync
t = upsample(t, L);
t = conv(t, pt);

ideal = t;
[corr_id, lags_id] = xcorr(zt, ideal);
[ideal_max_value, ideal_timing_index] = max(abs(corr_id));
ideal_timing_offset = lags_id(ideal_timing_index);

y_sync = y(ideal_timing_offset:end); % signal starts at the time sync bits
% figure(111);
% plot(real(y_sync),'b')
% hold on
% plot(imag(y_sync),'r')



%% Frame sync
[fcorr, flags] = xcorr(y_sync, f);
[~, frame_index] = max(abs(fcorr));
frame_offset = flags(frame_index);

y_fsync = y_sync(frame_offset + L:end);

ttt = y_sync(1:L:end);
ttt = ttt(length(t_modulated) + length(f_modulated) + 1:end);


[pcorr, plags] = xcorr(y_fsync, p);
[~, p_indexes] = maxk(abs(pcorr), n);
p_indexes = sort(p_indexes);

offset_indexes = zeros(1, n);
for i = 1:n
    offset_indexes(i) = plags(p_indexes(i));
end

y_psync = y_fsync(offset_indexes(1) + L:end);

%chunk_size = offset_indexes(2) - offset_indexes(1);
%chunk_size = chunk_size / L; % pilot + message
chunk_size = 20520/2 + 50;
%% DOWNSAMPLE
z_k = ttt; % this starts at the time symbol seq

z_k = z_k(1:chunk_size * n); % remove noise

test = z_k(1 + length(modulated_pilot):end);
test = demodulate_4qam(test);

figure(54)
subplot(2,1,1);
recovered_image = reshape(test, [171, 120]);
imshow(recovered_image);
subplot(2,1,2);
imshow(message);


%% Extract first chunk
%start_first_chunk = frame_offset + length(f_modulated); % maybe use different value?
%y_fsynced = z_k(start_first_chunk+1:end); 
%z_k = y_fsynced; % start of message
L2 = 3;
L1 = -3;

message_size = chunk_size - length(modulated_pilot);

chunks = zeros(message_size, n);
cs = zeros(1, message_size * n);
filters = zeros(L2-L1, n);

gamma = 0.00001;
trained_w = zeros(1, L2-L1);
delta = 1;
%% Train w for each chunk
for i = 1:n
    pilot = z_k(delta:delta + length(modulated_pilot) - 1);
    trained_w = LMS(trained_w, pilot, modulated_pilot, delta, gamma);
    filters(:,i) = trained_w; % train their respective filters
    c = z_k(delta + length(modulated_pilot):delta + chunk_size - 1);
    cs((i-1)*message_size + 1:i*message_size) = c;
    chunks(:,i) = transpose(c);
    delta = delta + chunk_size;
end


before_equalizations = reshape(chunks, 1, []);
zk_equalized = zeros(1, message_size * n);
% Equalize
for i = 1:n
    w = filters(:,i);
    current_chunk = chunks(:,i);
    %vk = filter(w, 1, current_chunk);
    vk = conv(current_chunk, w);
    %figure(10000 + i);
    %scatter(real(vk), imag(vk));
    zk_equalized((i - 1)*message_size + 1: i*message_size) = vk(1:message_size);
end

%% Soft decoding (TODO)

%% Demodulate
z_demodulated = demodulate_4qam(zk_equalized);

%% BER
%message = imread("shannon1440.bmp");
message = imread("shannon20520.bmp");

message_vec = reshape(message, 1, []);

bits = message_vec;

BER = mean(z_demodulated ~= bits);
disp(['BER is ', num2str(BER)])


%% Plots

%transmited and received signals and pt
% t_transmitted = [0:length(transmitsignal)-1] / Fs;
% t_received = [0:length(receivedsignal)-1] / Fs;
% figure(1)
% clf
% subplot(2,1,1)
% plot(t_transmitted, real(transmitsignal),'b')
% hold on
% plot(t_transmitted, imag(transmitsignal),'r')
% legend('real','imag')
% ylabel('xI(t)  and  xQ(t)')
% xlabel('Time in microseconds')
% subplot(2,1,2)
% plot(t_received, real(receivedsignal),'b')
% hold on
% plot(t_received, imag(receivedsignal),'r')
% zoom xon
% legend('real','imag')
% ylabel('yI(t)  and  yQ(t)')
% xlabel('Time in microseconds')

% figure(2)
% clf
% subplot(2,1,1)
% plot(([0:length(transmitsignal)-1]/length(transmitsignal)-0.5) * Fs / 10^6, abs(fftshift(fft(transmitsignal))))
% ylabel('abs(X(f))')
% xlabel('Frequency in MHz')
% subplot(2,1,2)
% plot(([0:length(receivedsignal)-1]/length(receivedsignal)-0.5) * Fs / 10^6, abs(fftshift(fft(receivedsignal))))
% ylabel('abs(Y(f))')
% xlabel('Frequency in MHz')

% figure(3)
% clf
% subplot(2,1,1);
% plot(t_p, pt, 'b');
% xlabel('Time in microseconds');
% ylabel('Pulse function p(t)');
% hold on;
% subplot(2,1,2);
% plot(([0:length(pt)-1]/length(pt)-0.5) * Fs / 10^6, abs(fftshift(fft(pt))))
% xlabel('Frequency in MHz');
% ylabel('abs(P(f))');


% recovered image
figure(4)
subplot(2,1,1);
recovered_image = reshape(z_demodulated, [171, 120]);
imshow(recovered_image);
subplot(2,1,2);
imshow(message);

%zk - before equalization
figure(5);
plot(real(before_equalizations), imag(before_equalizations), 'rx')

% vk - after equalization
figure(6);
plot(real(zk_equalized), imag(zk_equalized), 'rx');
% 
% y_synced = y_sync;
% 
% t_synced = [1:length(y_synced)] / Fs * 10^6;

% y after time synchronization
% figure(7)
% clf
% subplot(2,1,1)
% plot(t_synced, real(y_synced),'b')
% hold on
% plot(t_synced, imag(y_synced),'r')
% legend('real','imag')
% ylabel('yI(t)  and  yQ(t)')
% xlabel('Time in microseconds')
% subplot(2,1,2);
% plot(([0:length(y_synced)-1]/length(y_synced)-0.5) * Fs / 10^6, abs(fftshift(fft(y_synced))))
% xlabel('Frequency in MHz');
% ylabel('abs(P(f))');

% z demodulated
%figure(8)
%clf
%subplot(2,1,1);
%plot(z_demodulated(1:length(bits)), 'r');
%subplot(2,1,2);
%plot(bits, 'b')

