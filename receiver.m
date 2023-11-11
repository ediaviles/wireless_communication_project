% sampling freq.:
Fs = 200 * 10^6;

% symbol period:
T = 1/20 * 10^-6;

% oversampling factor:
L = floor(Fs*T);

N = 51; % Length of filter in symbol periods.
Ns = floor(N*L); % Number of filter samples

y_base = receivedsignal;

pt = sinc([-floor(Ns/2):Ns-floor(Ns/2)-1]/L); pt = transpose(pt)/norm(pt)/sqrt(1/(L)); %need to modify this

% Time Synchronization
[y_corr, lag] = xcorr(y_base', timing_sync_bits);
[~, max_index] = max(abs(y_corr));
timing_offset = lag(max_index);

y_time_sync = y_base(timing_offset + 1:end); 

% Equalization
p = y_time_sync(length(timing_sync_bits) + 1:length(timing_sync_bits) + length(pilot_sequence)); % find the pilot sequence in the transmitted signal
one_tap = (conj(pilot_sequence)*p) / (conj(pilot_sequence)*pilot_sequence');

y_eq = y_time_sync / one_tap;

% Frame Synchronization
[y_fcorr, lags] = xcorr(y_eq', fsync_bits);
[~, max_findex] = max(abs(y_fcorr));
framing_offset = lags(max_findex);

y = y_eq(framing_offset + length(fsync_bits) + 1:end);

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