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

% Synchronization
[y_corr, lag] = xcorr(y_base, timing_sync_bits);
[~, max_index] = max(abs(y_corr));
timing_offset = lag(max_index);

index = timing_offset + length(timing_sync_bits) + length(pilot_sequence); % where our signal starts
y_sync = y_base(index + 1:end); % Check this

% Equalization
p = y_base(timing_offset + length(timing_sync_bits) + 1:index); % find the pilot sequence in the transmitted signal
one_tap = (pilot_sequence*p) / (pilot_sequence*pilot_sequence');

y = y_sync / one_tap;


y_I = real(y);
y_Q = imag(y);

% filter using matched filter
matched_filter = fliplr(pt);
z_I = conv(y_I, matched_filter);
z_Q = conv(y_Q, matched_filter);

% Sample
z_Ik = z_I(1:L:end);
z_Qk = z_Q(1:L:end);

z_k = z_Ik + 1i * z_Qk;

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