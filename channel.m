sigman = 0.2;

receivedsignal = transmitsignal + sigman/sqrt(2) * (randn(size(transmitsignal))+j*randn(size(transmitsignal)));
%T o test the effect of phase offset and delay, you could simulate such a channel as
padding = (randn(1,1000) > 0.5) * 2 - 1;
transmitsignalwithdelay = [zeros(2147,1); transmitsignal'; padding'];
receivedsignal = exp(j*pi/6) * transmitsignalwithdelay + sigman/sqrt(2) * (randn(size(transmitsignalwithdelay))+j*randn(size(transmitsignalwithdelay)));

y_base = receivedsignal;
y_symbols = (receivedsignal > 0) * 2 - 1; % turn to symbols

pt = sinc([-floor(Ns/2):Ns-floor(Ns/2)-1]/L); pt = transpose(pt)/norm(pt)/sqrt(1/(L)); %need to modify this


% Time Synchronization
[y_corr, lag] = xcorr(y_symbols, timing_sync_bits);
[~, max_index] = max(abs(y_corr));
timing_offset = lag(max_index);

% Equalization
p = y_symbols(timing_offset + length(timing_sync_bits) + 1:timing_offset + length(timing_sync_bits) + length(pilot_sequence)); % find the pilot sequence in the transmitted signal
%p = conv(fliplr(pt), p);
%ps = p(1:Lend);
one_tap = (conj(pilot_sequence)*p) / (conj(pilot_sequence)*pilot_sequence');

y_base = y_base / one_tap;

% Frame Synchronization
[y_fcorr, lags] = xcorr(y_symbols, fsync_bits);
[~, max_findex] = max(abs(y_fcorr));
framing_offset = lags(max_findex);

y = y_base(framing_offset + length(fsync_bits) + 1:end);

y_I = real(y);
y_Q = imag(y);

% filter using matched filter
matched_filter = pt(end:-1:1);
y_conv = conv(y, matched_filter);
z_I = conv(y_I, matched_filter);
z_Q = conv(y_Q, matched_filter);

% Sample
y_sample = y_conv(1:L:end);
z_Ik = z_I(1:L:end);
z_Qk = z_Q(1:L:end);

z_k = z_Ik + j * z_Qk;

% demodulate (threshold)
y_demo = y_sample > 0;
z_demodulated = z_k > 0;


figure(1)
plot(real(z_demodulated));

message = imread("shannon1440.bmp");
message_vec = reshape(message, 1, []);

bits = transpose(message_vec);
figure(2)
plot(real(bits));

BER = mean(z_demodulated(1:length(bits)) ~= bits);
disp(['BER is ', num2str(BER)])