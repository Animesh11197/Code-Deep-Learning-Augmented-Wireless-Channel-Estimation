% Deep Residual Learning Meets OFDM Channel Estimation
% SNR represents Es/N0, Es/N0 = Eb/N0 * log2(M)
% For SNR, the power of pilot is ignored when calculating the varience of
% noise and the pilots suffer from the noise on equal level with data
% Power is balanced since IFFT/FFT is applied
% InsertDCNull is true and the effect caused by DCNull is considered in the
% SNR_OFDM, which is adjusted to be SNR + 10 * log(Num_of_subcarriers_used
% / Num_of_FFT)
% h_channel is different for each frame and stored in H, or read from h_set
%clear
%clc
%close all;

SNR_Range = -5:5:25;
Num_of_frame_each_SNR = 500;

MSE_LS_over_SNR = zeros(length(SNR_Range), 1);
%
MSE_LS_over_SNR_Bilinear = zeros(length(SNR_Range), 1);
MSE_MMSE_over_SNR = zeros(length(SNR_Range), 1);
MSE_MMSE_over_SNR_2 = zeros(length(SNR_Range), 1);% testing RS at pilot positions in MMSE

MSE_DNN_over_SNR = zeros(length(SNR_Range), 1);
MSE_ReEsNet_over_SNR = zeros(length(SNR_Range), 1);
MSE_Improved_DNN_over_SNR = zeros(length(SNR_Range), 1);
MSE_Improved_DNN_over_SNR_2 = zeros(length(SNR_Range), 1);
MSE_Improved_DNN_over_SNR_3 = zeros(length(SNR_Range), 1);
MSE_LS_DNN_over_SNR=zeros(length(SNR_Range), 1);
MSE_LS_DNN_concat_over_SNR=zeros(length(SNR_Range), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Added by Asrar 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SER_LS = zeros(length(SNR_Range), 1);
SER_MMSE = zeros(length(SNR_Range), 1);
SER_DNN = zeros(length(SNR_Range), 1);
SER_ReEsNet = zeros(length(SNR_Range), 1);
SER_Improved_DNN = zeros(length(SNR_Range), 1);
SER_Improved_DNN_2 = zeros(length(SNR_Range), 1);
SER_Improved_DNN_3 = zeros(length(SNR_Range), 1);
SER_LS_DNN=zeros(length(SNR_Range), 1);
SER_LS_DNN_concat=zeros(length(SNR_Range), 1);



BER_LS = zeros(length(SNR_Range), 1);
BER_MMSE = zeros(length(SNR_Range), 1);
BER_DNN = zeros(length(SNR_Range), 1);
BER_ReEsNet = zeros(length(SNR_Range), 1);
BER_improved_DNN = zeros(length(SNR_Range), 1);
BER_improved_DNN_2 = zeros(length(SNR_Range), 1);
BER_improved_DNN_3 = zeros(length(SNR_Range), 1);
BER_LS_DNN=zeros(length(SNR_Range), 1);
BER_LS_DNN_concat=zeros(length(SNR_Range), 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Import Deep Neuron Network
% load('SL.mat'); % SL102_10filter SL
load('SL_ETU.mat');% for ETU channel model
% load('SL_EVA.mat');% for EVA channel model

load('Inter-ResNet_block3.mat')% inter-ResNet with 3 blocks
load('Inter_ResNet_block2.mat')% inter-ResNet with 2 blocks


% Import Deep Neuron Network
load('ReEsNet_50000.mat');

% Import Deep Neuron Network
load('ReEsNet.mat');

% Import LS_DNN_model_72_144
net=importKerasNetwork('doppler_97_model_72_144_bicubic.h5');
%net2=importKerasNetwork('new_doppler_97_model_concatenated.h5');
%net2=importKerasNetwork('newModel_1.0_frame2016.h5');
% net2_EPA=importKerasNetwork('newModel_2_frame2016_EPA_97.h5');% For EPA channel model
net2_ETU=importKerasNetwork('newModel_frame2016_ETU_97.h5');% For ETU channel model
% net2=importKerasNetwork('newModel_frame2016_EVA_97.h5');% For EVA channel model
% net2_EVA=importKerasNetwork('LSiDNN_EVA.h5');% For EVA channel model



% Rgg_value = zeros(73, 73, 4);

%%% DELETING PRE_EXISTING DATA FILES.!
%delete("test_input.dat");
delete("recieved_pilots.dat");
delete("Biliner_input.dat");
delete("test_golden_output.dat");
%delete("test_MSE_each_Frame.dat");
%delete("RS_user_24x2_real.dat");
%delete("RS_user_24x2_imag.dat");
%delete("RS_72x2_imag.dat");
%delete("RS_72x2_real.dat");
%delete("Received_pilot_LS_24x2_real.dat");
%delete("Received_pilot_LS_24x2_imag.dat");

% SNR_Range = 20;

for SNR = SNR_Range
 
    Num_of_symbols = 12;
    Num_of_pilot = 2;
    Frame_size = Num_of_symbols + Num_of_pilot;
    
    
    
    M = 4; % QPSK
    k = log2(M);
    
    Num_of_subcarriers = 72;
    Num_of_FFT = Num_of_subcarriers + 1;
    length_of_CP = 16;
    
    Pilot_interval = 7;
    Pilot_starting_location = 1;
    
    Pilot_location = [(1:3:Num_of_subcarriers)', (2:3:Num_of_subcarriers)'];
    Pilot_value_user = 1 + 1j;
    
    length_of_symbol = Num_of_FFT + length_of_CP;
    
    MaxDopplerShift = 97;
    
    Num_of_QPSK_symbols = Num_of_subcarriers * Num_of_symbols * Num_of_frame_each_SNR;
    Num_of_bits = Num_of_QPSK_symbols * k;
    
    LS_MSE_in_frame = zeros(Num_of_frame_each_SNR, 1);
    %
    LS_MSE_in_frame_Bilinear = zeros(Num_of_frame_each_SNR, 1);
    MMSE_MSE_in_frame = zeros(Num_of_frame_each_SNR, 1);
    MMSE_MSE_in_frame_2 = zeros(Num_of_frame_each_SNR, 1);% for testing RS @ pilot locations
    
    DNN_MSE_in_frame = zeros(Num_of_frame_each_SNR, 1);
    ReEsNet_MSE_in_frame = zeros(Num_of_frame_each_SNR, 1);
    Improved_DNN_MSE_in_frame = zeros(Num_of_frame_each_SNR, 1);
    Improved_DNN_MSE_in_frame_2 = zeros(Num_of_frame_each_SNR, 1);% for ResNet with 3 res blocks
    Improved_DNN_MSE_in_frame_3 = zeros(Num_of_frame_each_SNR, 1);% for ResNet with 2 res blocks
    LS_DNN_MSE_in_frame = zeros(Num_of_frame_each_SNR, 1);% LS_DNN
    LS_DNN_MSE_in_frame2 = zeros(Num_of_frame_each_SNR, 1);% LS_DNN_concat
    
    
    
    t_LS = zeros(Num_of_frame_each_SNR, 1);
    t_MMSE = zeros(Num_of_frame_each_SNR, 1);
    t_DNN = zeros(Num_of_frame_each_SNR, 1);
    t_Improved = zeros(Num_of_frame_each_SNR, 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Added by Asrar 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SER_LS(SNR_Range == SNR) = 0;
        SER_MMSE(SNR_Range == SNR) = 0;
        SER_DNN(SNR_Range == SNR) = 0;
        SER_ReEsNet(SNR_Range == SNR) = 0;
        SER_Improved_DNN(SNR_Range == SNR) = 0;
        SER_Improved_DNN_2(SNR_Range == SNR) = 0;
        SER_Improved_DNN_3(SNR_Range == SNR) = 0;
        SER_LS_DNN(SNR_Range == SNR) = 0;
        SER_LS_DNN_concat(SNR_Range == SNR) = 0;
    
        BER_LS(SNR_Range == SNR) = 0;
        BER_MMSE(SNR_Range == SNR) = 0;
        BER_DNN(SNR_Range == SNR) = 0;
        BER_ReEsNet(SNR_Range == SNR) = 0;
        BER_improved_DNN(SNR_Range == SNR) = 0;
        BER_improved_DNN_2(SNR_Range == SNR) = 0;
        BER_improved_DNN_3(SNR_Range == SNR) = 0;
        BER_LS_DNN(SNR_Range == SNR) = 0;
        BER_LS_DNN_concat(SNR_Range == SNR) = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%     Num_of_frame_each_SNR = 10000;
    for Frame = 1:Num_of_frame_each_SNR
    
        % Data generation
        N = Num_of_subcarriers * Num_of_symbols;
        data = randi([0 1], N, k);
        Data = reshape(data, [], 1);
        dataSym = bi2de(data);
        
        % QPSK modulator
        QPSK_symbol = QPSK_Modualtor(dataSym);
        QPSK_signal = reshape(QPSK_symbol, Num_of_subcarriers, Num_of_symbols);
        
        % Pilot inserted
        [data_in_IFFT, data_location] = Pilot_Insert(Pilot_value_user, Pilot_starting_location, Pilot_interval, Pilot_location, Frame_size, Num_of_FFT, QPSK_signal);
        [data_for_channel, ~] = Pilot_Insert(1, Pilot_starting_location, Pilot_interval, kron((1 : Num_of_subcarriers)', ones(1, Num_of_pilot)), Frame_size, Num_of_FFT, (ones(Num_of_subcarriers, Num_of_symbols)));
        data_for_channel(1, :) = 1;
        
        % OFDM Transmitter
        [Transmitted_signal, ~] = OFDM_Transmitter(data_in_IFFT, Num_of_FFT, length_of_CP);
        [Transmitted_signal_for_channel, ~] = OFDM_Transmitter(data_for_channel, Num_of_FFT, length_of_CP);
        
        %% Channel
        
        % AWGN Channel
        SNR_OFDM = SNR + 10 * log10((Num_of_subcarriers / Num_of_FFT));
        %awgnChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
        %awgnChan.SNR = SNR_OFDM;
        
        %% Multipath Rayleigh Fading Channel
        
        PathDelays = [0 50 120 200 230 500 1600 2300 5000] * 1e-9; % ETU
        AveragePathGains = [-1.0 -1.0 -1.0 0.0 0.0 0.0 -3.0 -5.0 -7.0]; % ETU
%         PathDelays = [0 30 70 90 110 190 410] * 1e-9; % EPA
%         AveragePathGains = [0 -1 -2 -3 -8 -17.2 -20.8]; % EPA
%         PathDelays = [0 30 150 310 370 710 1090 1730 2510] * 1e-9; % EVA
%         AveragePathGains = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9]; % EVA
        
        %PathDelays = 1e-9.*[0, 1, 2, 100, 101, 200, 201,202,300,301,302];%VTV-EX
        %AveragePathGains = [0,0,0,-6.3,-6.3,-25.1,-25.1, -25.1,-22.7,-22.7,-22.7];%VTV-EX
        %PathDelays = 1e-9.*[0, 1, 2, 100, 101, 102, 200, 201, 300, 301, 400, 401];%RTV-EX
        %AveragePathGains = [0, 0, 0, -9.3, -9.3, -9.3, -20.3, -20.3, -21.3, -21.3, -28.8,-28.8];%RTV-EX 

        rayleighchan = comm.RayleighChannel(...
            'SampleRate', 1065000, ...
            'PathDelays', PathDelays, ...
            'AveragePathGains', AveragePathGains, ...
            'MaximumDopplerShift', randi([0, MaxDopplerShift]), ...
            'PathGainsOutputPort', true);
        
        chaninfo = info(rayleighchan); 
        coeff = chaninfo.ChannelFilterCoefficients;
        Np = length(rayleighchan.PathDelays);
        state = zeros(size(coeff, 2) - 1, size(coeff, 1)); % initializing the delay filter state
        
        [Multitap_Channel_Signal_user, Path_gain] = rayleighchan(Transmitted_signal);
        fracdelaydata = zeros(size(Transmitted_signal, 1), Np);
        
        for j = 1 : Np
            [fracdelaydata(:,j), state(:,j)] = filter(coeff(j, :), 1, Transmitted_signal_for_channel, state(:,j)); % fractional delay filter state is taken care of here.
        end
        
        SignalPower = mean(abs(Multitap_Channel_Signal_user) .^ 2);
        Noise_Variance = SignalPower / (10 ^ (SNR_OFDM / 10));
        
        Nvariance = sqrt(Noise_Variance / 2);
        n = Nvariance * (randn(length(Transmitted_signal), 1) + 1j * randn(length(Transmitted_signal), 1)); % Noise generation
        
        Multitap_Channel_Signal = Multitap_Channel_Signal_user + n;
        
        Multitap_Channel_Signal_user_for_channel = sum(Path_gain .* fracdelaydata, 2);
        
        %% OFDM Receiver
        [Unrecovered_signal, RS_User] = OFDM_Receiver(Multitap_Channel_Signal, Num_of_FFT, length_of_CP, length_of_symbol, Multitap_Channel_Signal_user);
        [~, RS] = OFDM_Receiver(Multitap_Channel_Signal_user_for_channel, Num_of_FFT, length_of_CP, length_of_symbol, Multitap_Channel_Signal_user_for_channel);
        Pilot_location_symbols = Pilot_starting_location : Pilot_interval : Frame_size;
        
        [Received_pilot, ~] = Pilot_extract(RS_User(2:end, :), Pilot_location, Num_of_pilot, Pilot_location_symbols, data_location);
        H_Ref = Received_pilot ./ Pilot_value_user;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Added by Asrar
                data_loc = [2 3 4 5 6 7 9 10 11 12 13 14];
                data_sym = Unrecovered_signal(2:end,data_loc);
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %% Channel estimation
        
        % Perfect knowledge on Channel
        
        %% LS
        [Received_pilot_LS, ~] = Pilot_extract(Unrecovered_signal(2:end, :), Pilot_location, Num_of_pilot, Pilot_location_symbols, data_location);
        
        H_LS = Received_pilot_LS / Pilot_value_user;
        
        H_LS_frame = imresize(H_LS, [Num_of_subcarriers, Frame_size]);
        H_LS_frame_BiLinear = imresize(H_LS, [Num_of_subcarriers, Frame_size],"bilinear");
        
        MSE_LS_frame = mean(abs(H_LS_frame - RS(2:end, :)).^2, 'all');
        MSE_LS_frame_Bilinear = mean(abs(H_LS_frame_BiLinear - RS(2:end, :)).^2, 'all');
        
recieved_pilots(:, :, 1) = real(Received_pilot_LS);
recieved_pilots(:, :, 2) = imag(Received_pilot_LS);

extracted_pilots(:, :, 1) = real(H_LS);
extracted_pilots(:, :, 2) = imag(H_LS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dlmwrite("recieved_pilots.dat",recieved_pilots,"-append",'delimiter','\n');
%%%%%%%% Writing DNN inputs for all frames into a File %%%%%%%%%%%%
dlmwrite("Biliner_input.dat",extracted_pilots,"-append",'delimiter','\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %            Added by Asrar to get dataset
        %%%%%%%%%%-----------------------------------------------%%%%%%%%%%
        %                       CIR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %extract pilot symbols
        out_pilot_sym(:,1) =  RS(2:end,1);
        out_pilot_sym(:,2) =  RS(2:end,8);
        %extraxt pilots
        out_pilot_sc(:,1) = out_pilot_sym(Pilot_location(:,1),1);
        out_pilot_sc(:,2) = out_pilot_sym(Pilot_location(:,2),2);
        
        out_new(1:24,:)=real(out_pilot_sc);
        out_new(25:2*24,:)=imag(out_pilot_sc);

        out_new_vector(1:2*24,:) =  out_new(:,1);
        out_new_vector(2*24+1:4*24,:) =  out_new(:,2);   %96x1 vector : pil1 real;imag;pil2 real;imag

        %CIR of entire frame 
        cir_frame_comp(1:72,:)=real(RS(2:end,:));
        cir_frame_comp(73:2*72,:)=imag(RS(2:end,:));

        cir_frame = reshape(cir_frame_comp,2016,1);        %2016x1 vector : sym1 real;imag; sym2 real;imag; ... ; sym14 real;imag


        HLS_pilot(1:24,:)=real(H_LS);
        HLS_pilot(25:2*24,:)=imag(H_LS);

        HLS_pilot_vector(1:2*24,:) =  HLS_pilot(:,1);
        HLS_pilot_vector(2*24+1:4*24,:) =  HLS_pilot(:,2);   %96x1 vector : pil1 real;imag;pil2 real;imag

        X(:,Frame) = HLS_pilot_vector;
%         Y(:,Frame) = out_new_vector;
        Y(:,Frame) = cir_frame;

        %%%%%%%%%%-----------------------------------------------%%%%%%%%%%
        %%                       LSDNN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        H_LS_DNN_feature2 = (predict(net2_ETU, HLS_pilot_vector'))';
        lsdnn_real = reshape(H_LS_DNN_feature2,144,14);

        H_LS_DNN2 = complex(lsdnn_real(1:72,:),lsdnn_real(73:144,:));

        MSE_LS_DNN_frame2 = mean(abs(H_LS_DNN2 - RS(2:end, :)).^2, 'all');


        % ------------------------- END ----------------------------------% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        recieved_pilots(:, :, 1) = real(Received_pilot_LS);
        recieved_pilots(:, :, 2) = imag(Received_pilot_LS);
        
        extracted_pilots(:, :, 1) = real(H_LS);
        extracted_pilots(:, :, 2) = imag(H_LS);
        
        
        %% MMSE
        
        % linear MMSE
        
        H_MMSE = zeros(Num_of_subcarriers, Num_of_pilot);
        H_MMSE_2 = zeros(Num_of_subcarriers, Num_of_pilot); %// MMSE using RS(2:end,i) only at piot positions
        
        for i = 1:size(Pilot_location, 2) 
        j=Pilot_location_symbols;
        H_pilot = H_Ref(:, i);
        Rhh = H_pilot * H_pilot';
        H_MMSE(:, i) = (RS(2:end, i) * (H_pilot')) * pinv(Rhh + (1 / 10^(SNR / 10)) * eye(size(Rhh, 1))) * H_LS(:, i);
        H_MMSE_2(:, i) = (RS(2:j:end, i) * (H_pilot')) * pinv(Rhh + (1 / 10^(SNR / 10)) * eye(size(Rhh, 1))) * H_LS(:, i);
        
      
        end
        
        H_MMSE_frame = imresize(H_MMSE, [Num_of_subcarriers, Frame_size]);
        H_MMSE_frame_2 = imresize(H_MMSE_2, [Num_of_subcarriers, Frame_size],'bilinear');
        
        MSE_MMSE_frame = mean(abs(H_MMSE_frame - RS(2:end, :)).^2, 'all');
        MSE_MMSE_frame_2 = mean(abs(H_MMSE_frame_2 - RS(2:end, :)).^2, 'all');
        
        %% ReEsNet A
        Res_feature_signal(:, :, 1) = real(Received_pilot_LS / Pilot_value_user);
        Res_feature_signal(:, :, 2) = imag(Received_pilot_LS / Pilot_value_user);
        
        H_DNN_feature = predict(ReEsNet_A, Res_feature_signal);
        
        H_DNN1 = H_DNN_feature(:, :, 1) + 1j * H_DNN_feature(:, :, 2);
        MSE_DNN_frame = mean(abs(H_DNN1 - RS(2:end, :)).^2, 'all');
        
        %% ReEsNet B
        Res_feature_signal(:, :, 1) = real(Received_pilot_LS / Pilot_value_user);
        Res_feature_signal(:, :, 2) = imag(Received_pilot_LS / Pilot_value_user);
        
        H_DNN_feature = predict(ReEsNet_B, Res_feature_signal);
        
        H_DNN2 = H_DNN_feature(:, :, 1) + 1j * H_DNN_feature(:, :, 2);
        MSE_DNN_frame_ReEsNet = mean(abs(H_DNN2 - RS(2:end, :)).^2, 'all');
        
        %% Improved Deep learning
        Res_feature_signal(:, :, 1) = real(Received_pilot_LS / Pilot_value_user);
        Res_feature_signal(:, :, 2) = imag(Received_pilot_LS / Pilot_value_user);
        
 
        
        H_DNN_feature = predict(SL_ETU, Res_feature_signal);%Improved_DNN_Trained
        H_DNN_feature2 = predict(Inter_ResNet_block3, Res_feature_signal);%Improved_DNN_Trained with 3 blocks
        H_DNN_feature3 = predict(Inter_ResNet_block2, Res_feature_signal);%Improved_DNN_Trained with 2 blocks
        
        
        
        H_DNN = H_DNN_feature(:, :, 1) + 1j * H_DNN_feature(:, :, 2);
        H_DNN_2 = H_DNN_feature2(:, :, 1) + 1j * H_DNN_feature2(:, :, 2);
        H_DNN_3 = H_DNN_feature3(:, :, 1) + 1j * H_DNN_feature3(:, :, 2);
        
        
        
        MSE_Improved_DNN_frame = mean(abs(H_DNN - RS(2:end, :)).^2, 'all');
        MSE_Improved_DNN_frame_2 = mean(abs(H_DNN_2 - RS(2:end, :)).^2, 'all');
        MSE_Improved_DNN_frame_3 = mean(abs(H_DNN_3 - RS(2:end, :)).^2, 'all');
        
%%%%%%%% Writing Golden Outputs for all frames into a File %%%%%%%%%%%%
out(:,:,1)=real(RS(2:end, :));
out(:,:,2)=imag(RS(2:end, :));
dlmwrite("test_golden_output.dat",out,"-append",'delimiter','\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% LS_DNN:
        DNN_input=cat(1,imresize(Res_feature_signal(:,:,1),[72,14]),imresize(Res_feature_signal(:,:,2),[72,14]));
        H_LS_DNN_feature=(predict(net,DNN_input'))';
        
        H_LS_DNN=H_LS_DNN_feature(1:72,:)+1j*H_LS_DNN_feature(73:144,:);
        
        MSE_LS_DNN_frame = mean(abs(H_LS_DNN - RS(2:end, :)).^2, 'all');
        
        % Concat input LS_DNN:
        %DNN_input_concat=cat(1,DNN_input,[ DNN_input(:,1) DNN_input(:,1:end-1)]);
        %H_LS_DNN_feature2=(predict(net2,DNN_input_concat'))';
        
        %H_LS_DNN2=H_LS_DNN_feature2(1:72,:)+1j*H_LS_DNN_feature2(73:144,:);
        
        %MSE_LS_DNN_frame2 = mean(abs(H_LS_DNN2 - RS(2:end, :)).^2, 'all');
        
%         % Pilot input LS_DNN:
%         LS_DNN_input=reshape(Res_feature_signal,96,[]);
%         H_LS_DNN_feature2 = (predict(net2, LS_DNN_input'))';
%         H_LS_DNN2_=reshape(H_LS_DNN_feature2,24,2,2);
%         H_LS_DNN2=zeros(72,14,2);
%         for i=1:2
%             H_LS_DNN2(:,:,i)=imresize(H_LS_DNN2_(:,:,i),[72 14]);
%         end
%         %H_LS_DNN2=imresize(H_LS_DNN2_,[72 14],2)
%         MSE_LS_DNN_frame2 = mean(abs(H_LS_DNN2 - RS(2:end, :)).^2, 'all');
        
        %% LS MSE calculation in each frame
        LS_MSE_in_frame(Frame, 1) = MSE_LS_frame;
        %
        LS_MSE_in_frame_Bilinear(Frame, 1) = MSE_LS_frame_Bilinear;
        
        %% MMSE MSE calculation in each frame
        MMSE_MSE_in_frame(Frame, 1) = MSE_MMSE_frame;
        %
        MMSE_MSE_in_frame_2(Frame, 1) = MSE_MMSE_frame_2;
        
        %% DNN MSE calculation in each frame
        DNN_MSE_in_frame(Frame, 1) = MSE_DNN_frame;
        
        %% DNN MSE calculation in each frame
        ReEsNet_MSE_in_frame(Frame, 1) = MSE_DNN_frame_ReEsNet;
        
        %% DNN MSE calculation in each frame
        Improved_DNN_MSE_in_frame(Frame, 1) = MSE_Improved_DNN_frame;
        %
        Improved_DNN_MSE_in_frame_2(Frame, 1) = MSE_Improved_DNN_frame_2;
        %
        Improved_DNN_MSE_in_frame_3(Frame, 1) = MSE_Improved_DNN_frame_3;
        
        %% LS_DNN MSE calculation in each frame
        LS_DNN_MSE_in_frame(Frame, 1) = MSE_LS_DNN_frame;
        %
        LS_DNN_MSE_in_frame2(Frame, 1) = MSE_LS_DNN_frame2;
        
        
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Added by Asrar
        %equalization
                LS_eq_dataSym = data_sym./H_LS_frame(:,data_location);
                MMSE_eq_dataSym = data_sym./H_MMSE_frame(:,data_location);
                DNN_eq_dataSym = data_sym./H_DNN1(:,data_location);
                ReEsNet_eq_dataSym = data_sym./H_DNN2(:,data_location);
                improved_DNN_eq_dataSym = data_sym./H_DNN(:,data_location);
                improved_DNN_2_eq_dataSym = data_sym./H_DNN_2(:,data_location);
                improved_DNN_3_eq_dataSym = data_sym./H_DNN_3(:,data_location);
                LS_DNN_eq_dataSym = data_sym./H_LS_DNN(:,data_location);
                LS_DNN2_eq_dataSym = data_sym./H_LS_DNN2(:,data_location);
        
        
        
                %demodulation
                qpskdemod = comm.QPSKDemodulator;
        
                LS_eq_dataSym_col = reshape(LS_eq_dataSym, [], 1);
                MMSE_eq_dataSym_col = reshape(MMSE_eq_dataSym, [], 1);
                DNN_eq_dataSym_col = reshape(DNN_eq_dataSym, [], 1);
                ReEsNet_eq_dataSym_col = reshape(ReEsNet_eq_dataSym, [], 1);
                improved_DNN_eq_dataSym_col = reshape(improved_DNN_eq_dataSym, [], 1);
                improved_DNN_2_eq_dataSym_col = reshape(improved_DNN_2_eq_dataSym, [], 1);
                improved_DNN_3_eq_dataSym_col = reshape(improved_DNN_3_eq_dataSym, [], 1);
                LS_DNN_eq_dataSym_col = reshape(LS_DNN_eq_dataSym, [], 1);
                LS_DNN2_eq_dataSym_col = reshape(LS_DNN2_eq_dataSym, [], 1);
        
        
                LS_demodSig = qpskdemod(LS_eq_dataSym_col);
                release(qpskdemod);
                MMSE_demodSig = qpskdemod(MMSE_eq_dataSym_col);
                release(qpskdemod);
                DNN_demodSig = qpskdemod(DNN_eq_dataSym_col);
                release(qpskdemod);
                ReEsNet_demodSig = qpskdemod(ReEsNet_eq_dataSym_col);
                release(qpskdemod);
                improved_DNN_demodSig = qpskdemod(improved_DNN_eq_dataSym_col);
                release(qpskdemod);
                improved_DNN_2_demodSig = qpskdemod(improved_DNN_2_eq_dataSym_col);
                release(qpskdemod);
                improved_DNN_3_demodSig = qpskdemod(improved_DNN_3_eq_dataSym_col);
                release(qpskdemod);
                LS_DNN_demodSig = qpskdemod(LS_DNN_eq_dataSym_col);
                release(qpskdemod);
                LS_DNN2_demodSig = qpskdemod(LS_DNN2_eq_dataSym_col);
                release(qpskdemod);
        
           
                SER_LS(SNR_Range == SNR) = SER_LS(SNR_Range == SNR) + length(find(dataSym ~= LS_demodSig))/(Num_of_subcarriers * Num_of_symbols);
                SER_MMSE(SNR_Range == SNR) =SER_MMSE(SNR_Range == SNR) + length(find(dataSym ~= MMSE_demodSig))/(Num_of_subcarriers * Num_of_symbols);
                SER_DNN(SNR_Range == SNR) = SER_DNN(SNR_Range == SNR) + length(find(dataSym ~= DNN_demodSig))/(Num_of_subcarriers * Num_of_symbols);
                SER_ReEsNet(SNR_Range == SNR) = SER_ReEsNet(SNR_Range == SNR) + length(find(dataSym ~= ReEsNet_demodSig))/(Num_of_subcarriers * Num_of_symbols);
                SER_Improved_DNN(SNR_Range == SNR) = SER_Improved_DNN(SNR_Range == SNR) + length(find(dataSym ~= improved_DNN_demodSig))/(Num_of_subcarriers * Num_of_symbols);
                SER_Improved_DNN_2(SNR_Range == SNR) = SER_Improved_DNN_2(SNR_Range == SNR) + length(find(dataSym ~= improved_DNN_2_demodSig))/(Num_of_subcarriers * Num_of_symbols);
                SER_Improved_DNN_3(SNR_Range == SNR) = SER_Improved_DNN_3(SNR_Range == SNR) + length(find(dataSym ~= improved_DNN_3_demodSig))/(Num_of_subcarriers * Num_of_symbols);
                SER_LS_DNN(SNR_Range == SNR) = SER_LS_DNN(SNR_Range == SNR) + length(find(dataSym ~= LS_DNN_demodSig))/(Num_of_subcarriers * Num_of_symbols);
                SER_LS_DNN_concat(SNR_Range == SNR) = SER_LS_DNN_concat(SNR_Range == SNR) + length(find(dataSym ~= LS_DNN2_demodSig))/(Num_of_subcarriers * Num_of_symbols);
        
        
        
                LS_Rx_Bits = de2bi(LS_demodSig);
                LS_Rx = reshape(LS_Rx_Bits, [], 1);
                MMSE_Rx_Bits = de2bi(MMSE_demodSig);
                MMSE_Rx = reshape(MMSE_Rx_Bits, [], 1);
                DNN_Rx_Bits = de2bi(DNN_demodSig);
                DNN_Rx = reshape(DNN_Rx_Bits, [], 1);
                ReEsNet_Rx_Bits = de2bi(ReEsNet_demodSig);
                ReEsNet_Rx = reshape(ReEsNet_Rx_Bits, [], 1);
                impDNN_Rx_Bits = de2bi(improved_DNN_demodSig);
                impDNN_Rx = reshape(impDNN_Rx_Bits, [], 1);
                impDNN_2_Rx_Bits = de2bi(improved_DNN_2_demodSig);
                impDNN_2_Rx = reshape(impDNN_2_Rx_Bits, [], 1);
                impDNN_3_Rx_Bits = de2bi(improved_DNN_3_demodSig);
                impDNN_3_Rx = reshape(impDNN_3_Rx_Bits, [], 1);
                LS_DNN_Rx_Bits = de2bi(LS_DNN_demodSig);
                LS_DNN_Rx = reshape(LS_DNN_Rx_Bits, [], 1);
                LS_DNN2_Rx_Bits = de2bi(LS_DNN2_demodSig);
                LS_DNN2_Rx = reshape(LS_DNN2_Rx_Bits, [], 1);
        
                
                BER_LS(SNR_Range == SNR) = BER_LS(SNR_Range == SNR) + length(find(Data ~= LS_Rx))/(Num_of_subcarriers * Num_of_symbols*k);
                BER_MMSE(SNR_Range == SNR) =BER_MMSE(SNR_Range == SNR) + length(find(Data ~= MMSE_Rx))/(Num_of_subcarriers * Num_of_symbols*k);
                BER_DNN(SNR_Range == SNR) = BER_DNN(SNR_Range == SNR) + length(find(Data ~= DNN_Rx))/(Num_of_subcarriers * Num_of_symbols*k);
                BER_ReEsNet(SNR_Range == SNR) = BER_ReEsNet(SNR_Range == SNR) + length(find(Data ~= ReEsNet_Rx))/(Num_of_subcarriers * Num_of_symbols*k);
                BER_improved_DNN(SNR_Range == SNR) = BER_improved_DNN(SNR_Range == SNR) + length(find(Data ~= impDNN_Rx))/(Num_of_subcarriers * Num_of_symbols*k);
                BER_improved_DNN_2(SNR_Range == SNR) = BER_improved_DNN_2(SNR_Range == SNR) + length(find(Data ~= impDNN_2_Rx))/(Num_of_subcarriers * Num_of_symbols*k);
                BER_improved_DNN_3(SNR_Range == SNR) = BER_improved_DNN_3(SNR_Range == SNR) + length(find(Data ~= impDNN_3_Rx))/(Num_of_subcarriers * Num_of_symbols*k);
                BER_LS_DNN(SNR_Range == SNR) = BER_LS_DNN(SNR_Range == SNR) + length(find(Data ~= LS_DNN_Rx))/(Num_of_subcarriers * Num_of_symbols*k);
                BER_LS_DNN_concat(SNR_Range == SNR) = BER_LS_DNN_concat(SNR_Range == SNR) + length(find(Data ~= LS_DNN2_Rx))/(Num_of_subcarriers * Num_of_symbols*k);
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    %% BER calculation
        SER_LS(SNR_Range == SNR) = SER_LS(SNR_Range == SNR)/ Num_of_frame_each_SNR;
        SER_MMSE(SNR_Range == SNR) =SER_MMSE(SNR_Range == SNR) / Num_of_frame_each_SNR;
        SER_DNN(SNR_Range == SNR) = SER_DNN(SNR_Range == SNR) / Num_of_frame_each_SNR;
        SER_ReEsNet(SNR_Range == SNR) = SER_ReEsNet(SNR_Range == SNR) / Num_of_frame_each_SNR;
        SER_Improved_DNN(SNR_Range == SNR) = SER_Improved_DNN(SNR_Range == SNR) / Num_of_frame_each_SNR;
        SER_Improved_DNN_2(SNR_Range == SNR) = SER_Improved_DNN_2(SNR_Range == SNR) / Num_of_frame_each_SNR;
        SER_Improved_DNN_3(SNR_Range == SNR) = SER_Improved_DNN_3(SNR_Range == SNR) / Num_of_frame_each_SNR;
        SER_LS_DNN(SNR_Range == SNR) = SER_LS_DNN(SNR_Range == SNR)/ Num_of_frame_each_SNR;
        SER_LS_DNN_concat(SNR_Range == SNR) = SER_LS_DNN_concat(SNR_Range == SNR)/ Num_of_frame_each_SNR;
    
    
        BER_LS(SNR_Range == SNR) = BER_LS(SNR_Range == SNR)/ Num_of_frame_each_SNR;
        BER_MMSE(SNR_Range == SNR) =BER_MMSE(SNR_Range == SNR) / Num_of_frame_each_SNR;
        BER_DNN(SNR_Range == SNR) = BER_DNN(SNR_Range == SNR) / Num_of_frame_each_SNR;
        BER_ReEsNet(SNR_Range == SNR) = BER_ReEsNet(SNR_Range == SNR) / Num_of_frame_each_SNR;
        BER_improved_DNN(SNR_Range == SNR) = BER_improved_DNN(SNR_Range == SNR) / Num_of_frame_each_SNR;
        BER_improved_DNN_2(SNR_Range == SNR) = BER_improved_DNN_2(SNR_Range == SNR) / Num_of_frame_each_SNR;
        BER_improved_DNN_3(SNR_Range == SNR) = BER_improved_DNN_3(SNR_Range == SNR) / Num_of_frame_each_SNR;
        BER_LS_DNN(SNR_Range == SNR) = BER_LS_DNN(SNR_Range == SNR)/ Num_of_frame_each_SNR;
        BER_LS_DNN_concat(SNR_Range == SNR) = BER_LS_DNN_concat(SNR_Range == SNR)/ Num_of_frame_each_SNR;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    MSE_LS_over_SNR(SNR_Range == SNR, 1) = sum(LS_MSE_in_frame, 1) / Num_of_frame_each_SNR;
    %
    MSE_LS_over_SNR_Bilinear(SNR_Range == SNR, 1) = sum(LS_MSE_in_frame_Bilinear, 1) / Num_of_frame_each_SNR;
    
    MSE_MMSE_over_SNR(SNR_Range == SNR, 1) = sum(MMSE_MSE_in_frame, 1) / Num_of_frame_each_SNR;
    %
    MSE_MMSE_over_SNR_2(SNR_Range == SNR, 1) = sum(MMSE_MSE_in_frame_2, 1) / Num_of_frame_each_SNR;
    
    
    MSE_DNN_over_SNR(SNR_Range == SNR, 1) = sum(DNN_MSE_in_frame, 1) / Num_of_frame_each_SNR;
    
    MSE_ReEsNet_over_SNR(SNR_Range == SNR, 1) = sum(ReEsNet_MSE_in_frame, 1) / Num_of_frame_each_SNR;
    
    MSE_Improved_DNN_over_SNR(SNR_Range == SNR, 1) = sum(Improved_DNN_MSE_in_frame, 1) / Num_of_frame_each_SNR;
    MSE_Improved_DNN_over_SNR_2(SNR_Range == SNR, 1) = sum(Improved_DNN_MSE_in_frame_2, 1) / Num_of_frame_each_SNR;
    MSE_Improved_DNN_over_SNR_3(SNR_Range == SNR, 1) = sum(Improved_DNN_MSE_in_frame_3, 1) / Num_of_frame_each_SNR;
    
    %LS_DNN
    MSE_LS_DNN_over_SNR(SNR_Range == SNR, 1) = sum(LS_DNN_MSE_in_frame, 1) / Num_of_frame_each_SNR;
    
    MSE_LS_DNN_concat_over_SNR(SNR_Range == SNR, 1) = sum(LS_DNN_MSE_in_frame2, 1) / Num_of_frame_each_SNR;

end

% Preamble_Error_Correction_Dataset.('X') =  X;
% Preamble_Error_Correction_Dataset.('Y') =  Y ;
% save('dataset_new',  'Preamble_Error_Correction_Dataset','-v7.3');

%////////////////////////////PLOT MSE//////////////////////////////////////////
figure(1);

semilogy(SNR_Range,MSE_LS_over_SNR ,'k--','LineWidth', 1.5);
hold on
semilogy(SNR_Range,MSE_MMSE_over_SNR_2,'r-', 'LineWidth', 1.5);
hold on
semilogy(SNR_Range,MSE_ReEsNet_over_SNR,'c-s', 'LineWidth', 1);
%hold on
%semilogy(SNR_Range,MSE_Improved_DNN_over_SNR_3,'m--', 'LineWidth', 1.5);
%hold on
%semilogy(SNR_Range,MSE_Improved_DNN_over_SNR_2,'--', 'LineWidth', 1.5);
hold on
semilogy(SNR_Range,MSE_Improved_DNN_over_SNR,'m-p', 'LineWidth', 1.5);
hold on
%semilogy(SNR_Range,MSE_LS_DNN_over_SNR,'g-+', 'LineWidth', 1);
%hold on
semilogy(SNR_Range,MSE_LS_DNN_concat_over_SNR,'b-h', 'LineWidth', 1.5);

ylim([0.0001, 1])

legend('MSE LS over SNR', ...
    'MMSE over SNR', ...
    'ReEsNet over SNR',...%'IResNet Block 2',...'IResNet Block 3',...
    'IResNet over SNR',...%'LS DNN',...
    'LS DNN ','location','southwest');
xlabel('SNR in dB');
ylabel('MSE');
grid on;
hold off;

%////////////////////////////PLOT SER//////////////////////////////////////////
figure(2);
semilogy(SNR_Range,SER_LS ,'k--','LineWidth', 1);
hold on
semilogy(SNR_Range,SER_MMSE,'r-', 'LineWidth', 1.5);
hold on
%semilogy(SNR_Range,SER_DNN,'g-+', 'LineWidth', 1.5);
%hold on
semilogy(SNR_Range,SER_ReEsNet,'c-s', 'LineWidth', 1.5);
hold on
semilogy(SNR_Range,SER_Improved_DNN,'m-p', 'LineWidth', 1.5);
% hold on
% semilogy(SNR_Range,SER_Improved_DNN_2,'--', 'LineWidth', 1.5);
% hold on
% semilogy(SNR_Range,SER_Improved_DNN_3,'m--', 'LineWidth', 1.5);
hold on
%semilogy(SNR_Range,SER_LS_DNN,'g-+', 'LineWidth', 1.5);
%hold on
semilogy(SNR_Range,SER_LS_DNN_concat,'b-h', 'LineWidth', 1.5);
ylim([0.0001, 1])

legend('SER LS over SNR', ...
    'SER MMSE over SNR', ...
    'SER ReEsNet over SNR',...
    'SER Improved DNN over SNR',...%'SER LS DNN',...'SER Inter-ResNet Block 3',...'SER Inter-ResNet Block 2',...
    'SER LS DNN ','location','southwest');
xlabel('SNR in dB');
ylabel('SER');
grid on;
hold off;

%////////////////////////////PLOT BER//////////////////////////////////////////
figure(3);
semilogy(SNR_Range,BER_LS ,'k--','LineWidth', 1.5);
hold on
semilogy(SNR_Range,BER_MMSE,'r-', 'LineWidth', 1.5);
hold on
%semilogy(SNR_Range,BER_DNN,'g-+', 'LineWidth', 1.5);
%hold on
semilogy(SNR_Range,BER_ReEsNet,'c-s', 'LineWidth', 1.5);
% hold on
% semilogy(SNR_Range,BER_improved_DNN_3,'m--', 'LineWidth', 1.5);
% hold on
% semilogy(SNR_Range,BER_improved_DNN_2,'--', 'LineWidth', 1.5);
hold on
semilogy(SNR_Range,BER_improved_DNN,'m-p', 'LineWidth', 1.5);
hold on
%semilogy(SNR_Range,BER_LS_DNN,'g-+', 'LineWidth', 1);
%hold on
semilogy(SNR_Range,BER_LS_DNN_concat,'b-h', 'LineWidth', 1.5);
ylim([0.0001, 1])

legend('BER LS over SNR', ...
    'BER MMSE over SNR', ...
    'BER ReEsNet over SNR',...   'BER IResNet Block 2',...'BER IResNet Block 3',...
    'BER IResNet',...%'BER LS DNN',...
    'BER LS DNN ','location','southwest');
xlabel('SNR in dB');
ylabel('BER');
grid on;
hold off;