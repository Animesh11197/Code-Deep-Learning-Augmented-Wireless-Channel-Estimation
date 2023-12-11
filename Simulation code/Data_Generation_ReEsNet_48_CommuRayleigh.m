% Channel Regression
% Data generation

function [Xtraining_Array_RSRP, Ytraining_regression_double_RSRP, Xvalidation_RSRP, Yvalidation_regression_double_RSRP,Ytraining_Pilots,Yvalidation_Pilots] = Data_Generation_ReEsNet_48_CommuRayleigh(Training_set_ratio, SNR_Range, Num_of_frame_each_SNR)

Num_of_subcarriers = 72; %126
Num_of_FFT = Num_of_subcarriers + 1;
length_of_CP = 16;

Num_of_symbols = 12;
Num_of_pilot = 2;
Frame_size = Num_of_symbols + Num_of_pilot;

Pilot_interval = 7;
Pilot_starting_location = 1;
Pilot_location = [(1:3:Num_of_subcarriers)', (2:3:Num_of_subcarriers)'];
Pilot_value_user = 1 + 1j;

length_of_symbol = Num_of_FFT + length_of_CP;

MaxDopplerShift = 97;

Xtraining_Array_RSRP = zeros(size(Pilot_location, 1), size(Pilot_location, 2), 2, Training_set_ratio * Num_of_frame_each_SNR * size(SNR_Range, 2));
Ytraining_regression_double_RSRP = zeros(Num_of_subcarriers, Frame_size, 2, Training_set_ratio * Num_of_frame_each_SNR * size(SNR_Range, 2));
Xvalidation_RSRP = zeros(size(Pilot_location, 1), size(Pilot_location, 2), 2, Num_of_frame_each_SNR * size(SNR_Range, 2) - Training_set_ratio * Num_of_frame_each_SNR * size(SNR_Range, 2));
Yvalidation_regression_double_RSRP = zeros(Num_of_subcarriers, Frame_size, 2, Num_of_frame_each_SNR * size(SNR_Range, 2) - Training_set_ratio * Num_of_frame_each_SNR * size(SNR_Range, 2));

%%%
Ytraining_Pilots = zeros(size(Pilot_location, 1), size(Pilot_location, 2), 2, Training_set_ratio * Num_of_frame_each_SNR * size(SNR_Range, 2));
Yvalidation_Pilots = zeros(size(Pilot_location, 1), size(Pilot_location, 2), 2, Num_of_frame_each_SNR * size(SNR_Range, 2) - Training_set_ratio * Num_of_frame_each_SNR * size(SNR_Range, 2));
%%%

for SNR = SNR_Range
        
for Frame = 1 : Num_of_frame_each_SNR

%% Data generation

QPSK_signal = ones(Num_of_subcarriers, Num_of_symbols);

% Pilot inserted
[data_in_IFFT, data_location] = Pilot_Insert(Pilot_value_user, Pilot_starting_location, Pilot_interval, Pilot_location, Frame_size, Num_of_FFT, QPSK_signal);
[data_for_channel, ~] = Pilot_Insert(1, Pilot_starting_location, Pilot_interval, kron((1 : Num_of_subcarriers)', ones(1, Num_of_pilot)), Frame_size, Num_of_FFT, QPSK_signal);
data_for_channel(1, :) = 1;

% OFDM Transmitter
[Transmitted_signal, ~] = OFDM_Transmitter(data_in_IFFT, Num_of_FFT, length_of_CP);
[Transmitted_signal_for_channel, ~] = OFDM_Transmitter(data_for_channel, Num_of_FFT, length_of_CP);

%% Channel

% AWGN Channel
SNR_OFDM = SNR + 10 * log10((Num_of_subcarriers / Num_of_FFT));
awgnChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
awgnChan.SNR = SNR_OFDM;

% Multipath Rayleigh Fading Channel

% PathDelays = [0 50 120 200 230 500 1600 2300 5000] * 1e-9; % ETU
% AveragePathGains = [-1.0 -1.0 -1.0 0.0 0.0 0.0 -3.0 -5.0 -7.0]; % ETU
% PathDelays = [0 30 70 90 110 190 410] * 1e-9; % EPA
% AveragePathGains = [0 -1 -2 -3 -8 -17.2 -20.8]; % EPA
PathDelays = [0 30 150 310 370 710 1090 1730 2510] * 1e-9; % EVA
AveragePathGains = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9]; % EVA
rayleighchan = comm.RayleighChannel(...
    'SampleRate', 1065000, ...
    'PathDelays', PathDelays, ...
    'AveragePathGains', AveragePathGains, ...
    'MaximumDopplerShift', randi([0, MaxDopplerShift]), ...
    'PathGainsOutputPort', true);

chaninfo = info(rayleighchan); 
coeff = chaninfo.ChannelFilterCoefficients;
Np = length(rayleighchan.PathDelays);
state = zeros(size(coeff, 2) - 1, size(coeff, 1)); %initializing the delay filter state

[Multitap_Channel_Signal_user, Path_gain] = rayleighchan(Transmitted_signal);
fracdelaydata = zeros(size(Transmitted_signal, 1), Np);

for j = 1 : Np
    [fracdelaydata(:,j), state(:,j)] = filter(coeff(j, :), 1, Transmitted_signal_for_channel, state(:,j)); %fractional delay filter state is taken care of here.
end

SignalPower = mean(abs(Multitap_Channel_Signal_user) .^ 2);
Noise_Variance = SignalPower / (10 ^ (SNR_OFDM / 10));

Nvariance = sqrt(Noise_Variance / 2);
n = Nvariance * (randn(length(Transmitted_signal), 1) + 1j * randn(length(Transmitted_signal), 1)); % Noise generation

Multitap_Channel_Signal = Multitap_Channel_Signal_user + n;

Multitap_Channel_Signal_user_for_channel = sum(Path_gain .* fracdelaydata, 2);

%% OFDM Receiver
[Received_signal, H_Ref] = OFDM_Receiver(Multitap_Channel_Signal, Num_of_FFT, length_of_CP, length_of_symbol, Multitap_Channel_Signal_user_for_channel);

Pilot_location_symbols = Pilot_starting_location : Pilot_interval : Frame_size;
[Received_pilot, ~] = Pilot_extract(Received_signal(2:end, :), Pilot_location, Num_of_pilot, Pilot_location_symbols, data_location);
H_LS = Received_pilot / Pilot_value_user;

if Frame <= fix(Training_set_ratio * Num_of_frame_each_SNR)
    Training_index = Frame + fix(Training_set_ratio * Num_of_frame_each_SNR) * (find(SNR_Range == SNR) - 1);
    Xtraining_Array_RSRP(:, :, 1, Training_index) = real(H_LS);
    Xtraining_Array_RSRP(:, :, 2, Training_index) = imag(H_LS);
    Ytraining_regression_double_RSRP(:, :, 1, Training_index) = real(H_Ref(2:end, :));
    Ytraining_regression_double_RSRP(:, :, 2, Training_index) = imag(H_Ref(2:end, :));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [extracted_Ytraining, ~] = Pilot_extract(Ytraining_regression_double_RSRP(:, :, 1, Training_index), Pilot_location, Num_of_pilot, Pilot_location_symbols, data_location);
%     Ytraining_Pilots(:,:,1,Training_index)=extracted_Ytraining;
%     [extracted_Ytraining, ~] = Pilot_extract(Ytraining_regression_double_RSRP(:, :, 2, Training_index), Pilot_location, Num_of_pilot, Pilot_location_symbols, data_location);
%     Ytraining_Pilots(:,:,2,Training_index)=extracted_Ytraining;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Validation_index = Frame - Training_set_ratio * Num_of_frame_each_SNR + (find(SNR_Range == SNR) - 1) * (Num_of_frame_each_SNR - Training_set_ratio * Num_of_frame_each_SNR);
    Xvalidation_RSRP(:, :, 1, Validation_index) = real(H_LS);
    Xvalidation_RSRP(:, :, 2, Validation_index) = imag(H_LS);
    Yvalidation_regression_double_RSRP(:, :, 1, Validation_index) = real(H_Ref(2:end, :));
    Yvalidation_regression_double_RSRP(:, :, 2, Validation_index) = imag(H_Ref(2:end, :));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [extracted_Yval, ~] = Pilot_extract(Yvalidation_regression_double_RSRP(:, :, 1, Validation_index), Pilot_location, Num_of_pilot, Pilot_location_symbols, data_location);
%     Yvalidation_Pilots(:,:,1,Validation_index)=extracted_Yval;
%     [extracted_Yval, ~] = Pilot_extract(Yvalidation_regression_double_RSRP(:, :, 2, Validation_index), Pilot_location, Num_of_pilot, Pilot_location_symbols, data_location);
%     Yvalidation_Pilots(:,:,2,Validation_index)=extracted_Yval;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end

end

end
