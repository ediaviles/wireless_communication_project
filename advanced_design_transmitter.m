clear
clc

message = imread("shannon1440.bmp");
%message = imread("shannon20520.bmp");

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
pilot_sequence = (randn(1, 100) > 0.5);
fsync_sequence = (randn(1, 100) > 0.5);
preamble = [frequency_sync_bits, timing_sync_bits, pilot_sequence, fsync_sequence]; %use this for clarity


% x_k divide into n chunks -> pilot, n_1, pilot, n_2 ...
n = 6; % number of chunks
chunks = reshape(message_vec, length(message_vec)/n, []);
xk = [preamble, chunks(:,1)'];
for i = 2:n
    xk = [xk,pilot_sequence,chunks(:,i)'];
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

%%
save('transmitsignal.mat','transmitsignal')

load receivedsignal.mat

% X axis values for plots %
t_transmitted = [1:length(transmitsignal)] / Fs * 10^6;
t_received = [1:length(receivedsignal)] / Fs * 10^6;
t_p = [-floor(Ns/2):Ns-floor(Ns/2)-1]*T*10^6/(L);


