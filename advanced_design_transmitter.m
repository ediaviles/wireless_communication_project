clear
clc

%message = imread("shannon1440.bmp");
message = imread("shannon20520.bmp");

message_vec = reshape(message, 1, []);

% M-QAM variables
M = 4;
b = log2(M);
d = 0.3;

% sampling freq.:
Fs = 200 * 10^6;

% symbol period:
T = 1/20 * 10^-6; % 1/T is sampling rate without oversampling, since we oversample by L, then the sampling rate becomes 200 MHz from transmitter

% oversampling factor:
L = floor(Fs*T);

N = 51; % Length of filter in symbol periods.
Ns = floor(N*L); % Number of filter samples


pt = sinc([-floor(Ns/2):Ns-floor(Ns/2)-1]/L); pt = transpose(pt)/norm(pt)/sqrt(1/(L));

frequency_sync_bits = ones(1, 100);
rng(4);
timing_sync_bits = (randn(1,100) > 0.5);
pilot_sequence = (randn(1, 500) > 0.5);
fsync_sequence = (randn(1, 100) > 0.5);
preamble = [frequency_sync_bits, timing_sync_bits, fsync_sequence]; %use this for clarity


% x_k divide into n chunks -> pilot, n_1, pilot, n_2 ...
n = 6; % number of chunks
chunk_size = length(message_vec) / n;
xk = [preamble];
for i = 1:n
    xk = [xk,pilot_sequence,message_vec((i-1)*chunk_size + 1:i*chunk_size)];
end

%% LDPC Encoding (TODO)

%% Modulate
xk = modulate_4qam(xk);
xk = xk * d;

%% Upsample
xk = upsample(xk, L);

%% Filter
xk = conv(xk, pt);
transmitsignal = xk;

% te = conv(transmitsignal, fliplr(pt));
% te = te(1:L:end);
% te = te(201:200 + chunk_size/b);
% te = qamdemod(te, 4);
% new_te = zeros(1, length(te) * 2);
% for i = 1:length(te)
%     guess = de2bi(te(i), 2, 'right-msb');
%     new_te((i - 1)*2 + 1: i*2) = guess;
% end
%te = demodulate_4qam(te);

% figure(58)
% subplot(2,1,1);
% recovered_image = reshape(ttt, [171, 120]);
% imshow(recovered_image);
% subplot(2,1,2);
% imshow(message);

%%
save('transmitsignal.mat','transmitsignal')

load receivedsignal.mat

% X axis values for plots %
t_transmitted = [1:length(transmitsignal)] / Fs * 10^6;
t_received = [1:length(receivedsignal)] / Fs * 10^6;
t_p = [-floor(Ns/2):Ns-floor(Ns/2)-1]*T*10^6/(L);


