sigman = 0.1;
%receivedsignal = transmitsignal + sigman/sqrt(2) * (randn(size(transmitsignal))+j*randn(size(transmitsignal)));
%T o test the effect of phase offset and delay, you could simulate such a channel as
%padding = (randn(1,1000) > 0.5) * 2 - 1;
rng(8)
transmitsignalwithdelay = [zeros(1, 2147), transmitsignal];
receivedsignal = exp(j*pi/6) * transmitsignalwithdelay + sigman/sqrt(2) * (randn(size(transmitsignalwithdelay))+j*randn(size(transmitsignalwithdelay)));
%receivedpacket1 = receivedsignal;
%receivedpacket2 = receivedsignal;
%save("receivedpacket1.mat", 'receivedsignal');
%save("receivedpacket2.mat", 'receivedsignal');

%load receivedsignal.mat

matched_filter = fliplr(pt);

t_received = [1:length(receivedsignal)] / Fs * 10^6;


y = receivedsignal;

t = timing_sync_bits; 
t_modulated = modulate_4qam(t);
t = t_modulated * d;


ps = pilot_sequence * d;
modulated_pilot = modulate_4qam(pilot_sequence);
p = modulated_pilot;
p = upsample(p, L);
p = conv(p, pt);

% f = fsync_sequence;
% f_modulated = modulate_4qam(f);
% f = f_modulated * d;
% f = upsample(f, L);
% f = conv(f, pt);

%% Filter
zt = conv(y, matched_filter);

%% Time sync
t = upsample(t, L);
t = conv(t, pt);

%t = conv(t, matched_filter);
pr = modulate_4qam(preamble);
modulated_preamble = pr;
pr = upsample(pr, L);
pr = conv(pr, pt);

[corr_id, lags_id] = xcorr(zt, pr);
[ideal_max_value, ideal_timing_index] = maxk(abs(corr_id), 1);
ideal_timing_offset = lags_id(ideal_timing_index);

%phase = angle(corr_id(ideal_timing_offset));

y_sync = y(ideal_timing_offset:end); % signal starts at the time sync bits
%y_sync = y_sync * exp(-j* phase);

%% Frame sync




%% DOWNSAMPLE
z_k = y_sync(1:L:end);
z_k = z_k(length(modulated_preamble) + 1:end); %point to first pilot

%chunk_size = (length(message_vec)/n/b + length(modulated_pilot))/rate;
chunk_size = (length(message_vec)/n/b + length(modulated_pilot));

z_k = z_k(1:chunk_size * n); % remove noise


%% Extract first chunk
L2 = 3;
L1 = -4;

message_size = chunk_size - length(modulated_pilot);

chunks = zeros(message_size, n);
filters = zeros(L2-L1, n);

gamma = 10;
%trained_w = zeros(1, L2-L1);
training = 129;
%% Train w for each chunk
for j = 1:training
    delta = 1; 
    if j <= 2
        gamma = gamma / 1.5;
    else
        gamma = gamma / 1.2;
    end
    for i = 1:n
        pilot = z_k(delta:delta + length(modulated_pilot) - 1);
        trained_w = filters(:, i).';
        trained_w = LMS(trained_w, pilot, modulated_pilot, delta, gamma);
        filters(:,i) = trained_w; % train their respective filters
        c = z_k(delta + length(modulated_pilot):delta + chunk_size - 1);
        chunks(:,i) = transpose(c);
        delta = delta + chunk_size;
    end
end

before_equalizations = reshape(chunks, 1, []);
zk_equalized = zeros(1, message_size * n);
% Equalize
for i = 1:n
    w = filters(:,i);
    current_chunk = chunks(:,i);
    vk = conv(current_chunk, w);
    zk_equalized((i - 1)*message_size + 1: i*message_size) = vk(1:message_size);
end

%% Soft decoding (TODO)
max_iterations = 100;
%zk_llr = compute_llr(zk_equalized, 0.1);
%zk_decoded = decode_llr(zk_llr, H, max_iterations);

%% Demodulate
%z_demodulated = demodulate_4qam(zk_equalized);
% z_demodulated = qamdemod(zk_equalized, M);
% 
% guess = ones(1, length(z_demodulated) * b);
% for i = 1:length(z_demodulated)
%     start_i = b*(i - 1) + 1;
%     end_i = b*(i);
%     decimal = z_demodulated(i);
%     binary = de2bi(decimal, b, 'left-msb');
%     guess(start_i:end_i) = binary;
% end
%z_demodulated = demodulate_4qam(zk_decoded);
z_demodulated = demodulate_4qam(zk_equalized);
%% BER
%message = imread("shannon1440.bmp");
%message = imread("shannon20520.bmp");

message = imread("shannon13720.bmp");

message_vec = reshape(message, 1, []);

%bits = message_vec;
bits = message_vec(1:length(message_vec));

BER = mean(z_demodulated ~= bits);
disp(['BER is ', num2str(BER)])


%% Plots

%transmited and received signals and pt
t_received = [0:length(receivedsignal)-1] / Fs;
t_transmitted = [0:length(transmitsignal)-1] / Fs;
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
recovered_image = reshape(z_demodulated, [140, 98]);
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

