clear
clc

message = imread("shannon1440.bmp");

message_vec = reshape(message, 1, []);

% BPSK:
xk = 2 * message_vec - 1; %multiply xk by some amplitude
xk = xk * 0.3;

% sampling freq.:
Fs = 200 * 10^6; 

% symbol period:
T = 1/20 * 10^-6; 

% oversampling factor:
L = floor(Fs*T);

N = 51; % Length of filter in symbol periods.
Ns = floor(N*L); % Number of filter samples


pt = sinc([-floor(Ns/2):Ns-floor(Ns/2)-1]/L); pt = transpose(pt)/norm(pt)/sqrt(1/(L)); % adjust sinc function according to bandwidth

xk_upsamp = upsample(xk, L);


xk_I = real(xk_upsamp);
xk_Q = imag(xk_upsamp);

x_I = conv(xk_I,  pt);
x_Q = conv(xk_Q,  pt);

x_base = x_I + 1j* x_Q;

transmitsignal = x_base;

figure(1)
clf
subplot(2,1,1)
plot(real(transmitsignal),'b')
hold on
plot(imag(transmitsignal),'r')
legend('real','imag')
ylabel('xI(t)  and  xQ(t)')
xlabel('Time in samples')


figure(2)
clf
subplot(2,1,1)
plot([0:length(transmitsignal)-1]/length(transmitsignal)-0.5, abs(fftshift(fft(transmitsignal))))
ylabel('abs(X(f))')
xlabel('Frequency in 1/samples')
