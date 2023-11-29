sigman = 0.2;

receivedsignal = transmitsignal + sigman/sqrt(2) * (randn(size(transmitsignal))+j*randn(size(transmitsignal)));
%T o test the effect of phase offset and delay, you could simulate such a channel as
%padding = (randn(1,1000) > 0.5) * 2 - 1;
transmitsignalwithdelay = [zeros(1, 2147), transmitsignal];
receivedsignal = exp(j*pi/6) * transmitsignalwithdelay + sigman/sqrt(2) * (randn(size(transmitsignalwithdelay))+j*randn(size(transmitsignalwithdelay)));


t_received = [1:length(receivedsignal)] / Fs * 10^6;

matched_filter = fliplr(pt);

y = receivedsignal;

t = timing_sync_bits;
% modulate t
t_modulated = ones(1, length(t)/b)
for i = 1:length(t_modulated)
    start_i = b*(i - 1) + 1
    end_i = b*(i)
    grouping = t(start_i:end_i);
    
    % Convert binary grouping to decimal
    decimal_value = bi2de(grouping, 'left-msb');

    % QAM modulation
    t_modulated(i) = qammod(decimal_value, M);
end
t = t_modulated * 0.3;


ps = pilot_sequence * 0.3;

% modulate f
f = fsync_sequence;
% modulate t
f_modulated = ones(1, length(f)/b)
for i = 1:length(f_modulated)
    start_i = b*(i - 1) + 1
    end_i = b*(i)
    grouping = f(start_i:end_i);
    
    % Convert binary grouping to decimal
    decimal_value = bi2de(grouping, 'left-msb');

    % QAM modulation
    f_modulated(i) = qammod(decimal_value, M);
end
f = f_modulated * 0.3;
f = upsample(f, L);
f = conv(f, pt);
f = f(1:L:end);

chunk_size = chunk * 0.3;
chunk_size = length(chunk_size);

%% Time sync
t = upsample(t, L);
t = conv(t, pt);

%ideal = xk(1001:2000);
ideal = t;
[corr_id, lags_id] = xcorr(y, ideal);
[ideal_max_value, ideal_timing_index] = max(abs(corr_id));
ideal_timing_offset = lags_id(ideal_timing_index);

%%
y_sync = y(ideal_timing_offset+1:end);
figure(111);
plot(real(y_sync),'b')
hold on
plot(imag(y_sync),'r')


%% Filter
y_filt = conv(y_sync, matched_filter);

%% DOWNSAMPLE
z_k = y_filt(1:L:end);


%%
delta = ideal_timing_offset + length(t)/L;
delta = delta / L;
first_pilot = z_k(delta + 1:delta + length(ps)/b);

y_synced = z_k(delta + 1:end);
%% 
[fcorr, flags] = xcorr(z_k, f);
[~, frame_index] = max(abs(fcorr));
frame_offset = flags(frame_index);

start_first_chunk = frame_offset + length(f);
y_fsynced = z_k(start_first_chunk+1:end);

first_chunk = y_fsynced(1:chunk_size);

%% Downsample
%tau = mod(length(y), L);
%y_synced2 = y(tau +1:end);
%y_downsampled = y(1:L:end);

%% Time sync
%[corr, lags] = xcorr(y_downsampled, t); % look at class slides which is based on analog signal
%[max_value, timing_index] = max(abs(corr));
%timing_offset = lags(timing_index);
%phase = angle(max_value);
%delta = timing_offset + length(t); % determine the start of the first pilot

%y_synced = y_downsampled(delta + 1:end);
%y_synced = y_synced * exp(-j*phase);

%first_pilot = y_downsampled(delta + 1:delta + length(ps)); % extract first pilot

%% Frame sync
%[fcorr, flags] = xcorr(y_downsampled, f);
%[~, frame_index] = max(abs(fcorr));
%frame_offset = lags(frame_index);

%start_first_chunk = frame_offset + length(f);
%y_fsynced = y_downsampled(start_first_chunk+1:end);

%first_chunk = y_fsynced(1:chunk_size);

%% Equalizer
%one_tap = (first_pilot*conj(modulated_pilot)') / (conj(modulated_pilot)*modulated_pilot');
%equalizations = [one_tap];

L1 = -8;
L2 = 8;
w = zeros(1, L2 - L1 + 1); % w is array of size L2 + abs(L1) + 1 -> Right now 21
u = first_pilot(end:-1:start);
v0 = conv(w, first_pilot);
e0 = v0 - modulated_pilot;
new_w = w;
step_size = 0.5;
for n = L2-L1+1 : N 
	u = inp(n:-1:n-sysorder+1) ;
    y(n)= w' * u;
    e(n) = d(n) - y(n) ;
% Start with big mu for speeding the convergence then slow down to reach the correct weights
    if n < 20
        mu=0.32;
    else
        mu=0.15;
    end
	w = w + mu * u * e(n) ;
end
%before_equalization = sample_first_chunk;

%first_chunk = first_chunk; % / one_tap;

delta = chunk_size;

%% Equalize based on number of chunks (n) and remove pilots from message
chunks = [first_chunk];

for i = 1:1:n-1
    % extract nth pilot
    pilot = y_fsynced(delta + 1:delta + length(ps)/b);

    % calculate one tap
    %one_tap = (pilot*conj(modulated_pilot)') / (conj(modulated_pilot)*modulated_pilot');
    %equalizations = [equalizations, one_tap];
    
    vk = conv(w, pilot);
    ek = vk - modulated_pilot;
    new_w = w;
    for j = 1:length(new_w)
        new_w(j) = w(j) - step_size * ek * conj(pilot(length(pilot) - j + 1));
    end
    w = new_w;

    start_of_chunk = delta + length(ps)/b;
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
z_demodulated = qamdemod(z_k, M);
zd = z_demodulated;
guess = ones(1, length(z_demodulated) * b);
for i = 1:length(z_demodulated)
    start_i = b*(i - 1) + 1
    end_i = b*(i)
    decimal = z_demodulated(i);
    binary = de2bi(decimal, b, 'left-msb');
    guess(start_i:end_i) = binary;
end
z_demodulated = guess;

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
plot(real(z_k(1:1440/b)), imag(z_k(1:1440/b)), 'rx');

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

