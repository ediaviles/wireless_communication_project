
%load('receivedpacket1.mat')
%load('receivedpacket2.mat')


packet1_demod = basic_receiver(receivedpacket1);
packet2_demod = basic_receiver(receivedpacket2);

entire_received = [packet1_demod, packet2_demod];

%% BER
%message = imread("shannon1440.bmp");
message = imread("shannon20520.bmp");

message_vec = reshape(message, 1, []);

bits = message_vec;

BER = mean(entire_received ~= bits);
disp(['BER is ', num2str(BER)])

figure(4)
subplot(2,1,1);
recovered_image = reshape(entire_received, [171, 120]);
imshow(recovered_image);
subplot(2,1,2);
imshow(message);