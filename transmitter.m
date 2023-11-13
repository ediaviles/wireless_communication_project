clear
clc

message = imread("shannon1440.bmp");

message_vec = reshape(message, 1, []);

% BPSK:

xk = 2 * message_vec - 1;

% sampling freq.:
Fs = 200 * 10^6;

% symbol period:
T = 1/20 * 10^-6;

% oversampling factor:
L = floor(Fs*T);

N = 51; % Length of filter in symbol periods.
Ns = floor(N*L); % Number of filter samples


pt = sinc([-floor(Ns/2):Ns-floor(Ns/2)-1]/L); pt = transpose(pt)/norm(pt)/sqrt(1/(L));

frequency_sync_bits = ones(1, 100);
rng(5);
timing_sync_bits = (randn(1,100) > 0.5) * 2 - 1;
pilot_sequence = (randn(1, 100) > 0.5) * 2 - 1;
fsync_sequence = (randn(1, 100) > 0.5) * 2 - 1;



preamble = [frequency_sync_bits, timing_sync_bits, pilot_sequence, fsync_sequence]; %use this for clarity

% x_k divide into n chunks -> pilot, n_1, pilot, n_2 ...
chunk = ones(1, 1440);
n = 1;
new_xk = [preamble];
chunk_size = floor(length(xk)/n);
new_xk = [new_xk, xk(1:chunk_size)]
for i = 1:1:n-1
    chunk = xk((i) * chunk_size + 1:(i) * chunk_size + chunk_size);
    new_xk = [new_xk, pilot_sequence, chunk];
end

xk = new_xk * 0.3;
xk = upsample(xk, L);
xk = conv(xk, pt);

transmitsignal = xk;

save('transmitsignal.mat','transmitsignal')

load receivedsignal.mat

% X axis values for plots %
t_transmitted = [1:length(transmitsignal)] / Fs * 10^6;
t_received = [1:length(receivedsignal)] / Fs * 10^6;


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
