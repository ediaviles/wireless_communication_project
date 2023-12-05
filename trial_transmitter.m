clear
clc

%message = imread("shannon1440.bmp");
%message = imread("shannon20520.bmp");
message = imread("shannon13720.bmp");


message_vec = reshape(message, 1, []);

%packet1 = message_vec(1:length(message_vec)/2);
%packet2 = message_vec(length(message_vec)/2 + 1:end);

%message_vec = packet1;

%message_vec = packet2;

% M-QAM variables
M = 4;
b = log2(M);
d = sqrt(2)/3/1.64;
%d = 0.315;

% sampling freq.:
Fs = 200 * 10^6;

% symbol period:
T = 1/20 * 10^-6; % 1/T is sampling rate without oversampling, since we oversample by L, then the sampling rate becomes 200 MHz from transmitter

% oversampling factor:
L = floor(Fs*T);

N = 51; % Length of filter in symbol periods.
Ns = floor(N*L); % Number of filter samples


pt = sinc([-floor(Ns/2):Ns-floor(Ns/2)-1]/L); pt = transpose(pt)/norm(pt)/sqrt(1/(L));

% 52,200,96,0,114

frequency_sync_bits = ones(1, 44);
rng(6);
timing_sync_bits = (randn(1, 100) > 0.5);
% 20 pilots ideal
pilot_sequence = (randn(1, 4 * 13) > 0.5);
fsync_sequence = (randn(1, 0) > 0.5);
preamble = [frequency_sync_bits, timing_sync_bits, fsync_sequence]; %use this for clarity


% x_k divide into n chunks -> pilot, n_1, pilot, n_2 ...
n = 14; % number of chunks
chunk_size = length(message_vec) / n;
xk = [preamble];
for i = 1:n
    xk = [xk,pilot_sequence,message_vec((i-1)*chunk_size + 1:i*chunk_size)];
end

%% LDPC Encoding (TODO)
% xk has length 14820  =  info_bits:
% xk has length 14512  =  info_bits:
info_bits = length(xk);
% rate:
r1 = 12;
r2 = 13;
rate = r1/r2;
% block_length;
block_length = info_bits / rate;
% num_constr:
num_constr = block_length - info_bits;
%%
% Gallager LDPC : did not work/not easy to convert into standard form
H = zeros(num_constr, block_length);
G = zeros(info_bits, block_length);
const = 2;
%w_c = (r2-r1)*const;
%w_r = (r2)*const;
%A_0 = zeros(num_constr/w_c, block_length);
% for i_row = 1 : num_constr/w_c
%     A_0(i_row, (i_row - 1)*w_r + 1:i_row*w_r) = 1;
% end
% for i_col = 1 : w_c
%     H((i_col-1)*num_constr/w_c + 1 : i_col * num_constr/w_c, :) =  A_0(:,randperm(size(A_0,2)));
% end

% Create H in standard form instead:
I_H = eye(num_constr);
H(:, info_bits+1:end) = I_H;

num_ones = 3*num_constr;
% Validate feasibility
if num_ones > num_constr * info_bits
    error('The number of ones exceeds the total number of elements in the matrix.');
end

% Initialize the matrix with zeros
P_H = sparse(num_constr, info_bits);

% Distribute the ones
for i = 1:num_ones
    placed = false;
    while ~placed
        row = randi(num_constr);
        col = randi(info_bits);
        if P_H(row, col) == 0
            P_H(row, col) = 1;
            placed = true;
        end
    end
end


% Convert H to G:
I_G = eye(info_bits);
P_G = P_H';
G(:, 1:info_bits) = I_G;
G(:, info_bits+1:end) = P_G;

encoded_xk = mod(xk * G, 2);


%encodercfg = ldpcEncoderConfig(H);

%codeword = ldpcEncode(xk, cfgLDPCEnc);

%% Modulate
xk_mod = modulate_4qam(encoded_xk);
%xk_mod = modulate_16qam(xk);
% xk_mod = zeros(1, length(xk)/b);
% for i = 1:length(xk_mod)
%     start_i = b*(i - 1) + 1;
%     end_i = b*(i);
%     grouping = xk(start_i:end_i);
% 
%     % Convert binary grouping to decimal
%     decimal_value = bi2de(grouping, 'left-msb');
%     % QAM modulation
%     xk_mod(i) = qammod(decimal_value, M);
% end


%% Upsample
xk_up = upsample(xk_mod, L);

%% Filter
xk_filt = conv(xk_up, pt);
transmitsignal = xk_filt;

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
%save('transmitpacket1.mat', 'transmitsignal')
%save('transmitpacket2.mat', 'transmitsignal')

%load receivedsignal.mat

% X axis values for plots %
% t_transmitted = [1:length(transmitsignal)] / Fs * 10^6;
% t_received = [1:length(receivedsignal)] / Fs * 10^6;
% t_p = [-floor(Ns/2):Ns-floor(Ns/2)-1]*T*10^6/(L);


