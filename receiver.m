% sampling freq.:
Fs = 200 * 10^6;

% symbol period:
T = 1/20 * 10^-6;

% oversampling factor:
L = floor(Fs*T);

N = 51; % Length of filter in symbol periods.
Ns = floor(N*L); % Number of filter samples

y_base = receivedsignal;

y_I = real(y_base);
y_Q = imag(y_base);

pt = sinc([-floor(Ns/2):Ns-floor(Ns/2)-1]/L); pt = transpose(pt)/norm(pt)/sqrt(1/(L)); %need to modify this

% filter using matched filter
matched_filter = fliplr(pt);
z_I = conv(y_I, matched_filter);
z_Q = conv(y_Q, matched_filter);

% Sample
z_Ik = z_I(1:L:length(z_I));
z_Qk = z_Q(1:L:length(z_Q));

% Synchronization


% Equalization


z_k = z_Ik + j * z_Qk;

% demodulate (threshold)
z_demodulated = z_k > 0;

% recover bits


% BER
message = imread("shannon1440.bmp");
message_vec = reshape(message, 1, []);

bits = transpose(message_vec);
BER = mean(z_demodulated(1:length(message_vec)) ~= bits);
disp(['BER is ', num2str(BER)])

% Plot constellation
figure(1)
stem([1:length(message_vec)], message_vec,'bx')
hold on
stem([1:length(z_demodulated)], z_demodulated,'ro')
legend('x_k', 'z_k')
xlabel('discrete time k')
axis tight