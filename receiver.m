% sampling freq.:
Fs = 200 * 10^6;

% symbol period:
T = 1/20 * 10^-6;

% oversampling factor:
L = floor(Fs*T);

N = 51; % Length of filter in symbol periods.
Ns = floor(N*L); % Number of filter samples

y = receivedsignal;
y_symbols = (receivedsignal > 0) * 2 - 1; % turn to symbols
y_symbols = y_symbols * 0.3;

t = timing_sync_bits * 0.3;
t = upsample(t, L);
t = conv(t, fliplr(pt));

f = fsync_bits * 0.3;
f = upsample(f, L);
f = conv(f, fliplr(pt));

ps = pilot_sequence * 0.3;
ps = upsample(ps, L);
ps = conv(ps, fliplr(pt));

% Time sync
%filter_timing_sync_bits = conv(pt, timing_sync_bits * 0.3);
%[corr, lags] = xcorr(y, filter_timing_sync_bits);
[corr, lags] = xcorr(y_symbols, t);
[~, timing_index] = max(abs(corr));
timing_offset = lags(timing_index);

% Eq
p = y_symbols(timing_offset + length(t) + 1:timing_offset + length(t) + length(ps)); % find the pilot sequence in the transmitted signal
%p = conv(p, fliplr(p));
one_tap = (conj(ps)*p) / (conj(ps)*ps');

y = y / one_tap;

% Frame sync
[fcorr, flags] = xcorr(y_symbols, f);
[~, f_index] = max(abs(fcorr));
f_offset = flags(f_index);

y = y(f_offset + length(f) + 1:end);


y_I = real(y);
y_Q = imag(y);

% filter using matched filter
matched_filter = fliplr(pt);
z_I = conv(y_I, matched_filter);
z_Q = conv(y_Q, matched_filter);

% Sample
z_Ik = z_I(1:L:end);
z_Qk = z_Q(1:L:end);

z_k = z_Ik + j * z_Qk;

% demodulate (threshold)
z_demodulated = z_k > 0;


% BER
message = imread("shannon1440.bmp");
message_vec = reshape(message, 1, []);

bits = transpose(message_vec);
BER = mean(z_demodulated(1:length(bits)) ~= bits);
disp(['BER is ', num2str(BER)])

% Plot constellation
figure(1)
stem([1:length(bits)], message_vec,'bx')
hold on
stem([1:length(bits)], z_demodulated(1:length(bits)),'ro')
legend('x_k', 'z_k')
xlabel('discrete time k')
axis tight