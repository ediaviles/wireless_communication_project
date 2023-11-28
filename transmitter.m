clear
clc

message = imread("shannon1440.bmp");

message_vec = reshape(message, 1, []);

% M-QAM:
M = 4;
b = log2(M);

% Preallocate array for encoded values:
xk = zeros(1, length(message_vec)/b);

for i = 1:length(xk)
    start_i = b*(i - 1) + 1
    end_i = b*(i)
    grouping = message_vec(start_i:end_i);
    
    % Convert binary grouping to decimal
    decimal_value = bi2de(grouping, 'left-msb');

    % QAM modulation
    xk(i) = qammod(decimal_value, M);
end

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
pilot_sequence = (randn(1, 32) > 0.5);
% modulate pilot
modulated_pilot = ones(1, length(pilot_sequence)/b);
for i = 1:length(modulated_pilot)
    start_i = b*(i - 1) + 1
    end_i = b*(i)
    grouping = pilot_sequence(start_i:end_i);
    
    % Convert binary grouping to decimal
    decimal_value = bi2de(grouping, 'left-msb');

    % QAM modulation
    modulated_pilot(i) = qammod(decimal_value, M);
end

fsync_sequence = (randn(1, 100) > 0.5);
preamble = [frequency_sync_bits, timing_sync_bits, pilot_sequence, fsync_sequence]; %use this for clarity


% modulate our preamble
modulated_preamble = ones(1, length(preamble)/b);
for i = 1:length(modulated_preamble)
    start_i = b*(i - 1) + 1
    end_i = b*(i)
    grouping = preamble(start_i:end_i);
    
    % Convert binary grouping to decimal
    decimal_value = bi2de(grouping, 'left-msb');

    % QAM modulation
    modulated_preamble(i) = qammod(decimal_value, M);
end




% x_k divide into n chunks -> pilot, n_1, pilot, n_2 ...
chunk = ones(1, 1440);
n = 6; % number of chunks
new_xk = [modulated_preamble];
chunk_size = floor(length(xk)/n);
new_xk = [new_xk, xk(1:chunk_size)]
for i = 1:1:n-1
    chunk = xk((i) * chunk_size + 1:(i) * chunk_size + chunk_size);
    new_xk = [new_xk, modulated_pilot, chunk];
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
t_p = [-floor(Ns/2):Ns-floor(Ns/2)-1]*T*10^6/(L);


