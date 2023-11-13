sigman = 0.2;

receivedsignal = transmitsignal + sigman/sqrt(2) * (randn(size(transmitsignal))+j*randn(size(transmitsignal)));
%T o test the effect of phase offset and delay, you could simulate such a channel as
padding = (randn(1,1000) > 0.5) * 2 - 1;
transmitsignalwithdelay = [zeros(1, 2147), transmitsignal, padding];
receivedsignal = exp(j*pi/6) * transmitsignalwithdelay + sigman/sqrt(2) * (randn(size(transmitsignalwithdelay))+j*randn(size(transmitsignalwithdelay)));

y = receivedsignal';
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
p = y(timing_offset + length(t) + 1:timing_offset + length(t) + length(ps)); % find the pilot sequence in the transmitted signal
%p = conv(p, fliplr(p));
%one_tap = (conj(ps)*p') / (conj(ps)*ps');
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
%z_k = -z_k;


% demodulate (threshold)
z_demodulated = z_k > 0;

figure(3);
plot(real(z_k(1:1440)), imag(z_k(1:1440)), 'rx');

figure(1);
plot(real(receivedsignal), imag(receivedsignal), 'rx');

figure(2);
plot(real(y(f_offset + length(f) + 1:end)));

% BER
message = imread("shannon1440.bmp");
message_vec = reshape(message, 1, []);

bits = message_vec';
figure(2);
plot(real(bits));
BER = mean(z_demodulated(1:length(bits)) ~= bits);
disp(['BER is ', num2str(BER)])