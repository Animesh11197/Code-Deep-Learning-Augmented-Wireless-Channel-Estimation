Data_Samples = 5000;
Training_Set_Rate = 0.995;
ValidationFrequency = 300;
SNR = 0:5:25;

[XTrain_RSRP, YTrain_RSRP, XValidation_RSRP, YValidation_RSRP,YTrain_Pilots_RSRP,YValidation_Pilots_RSRP] = Data_Generation_ReEsNet_48_CommuRayleigh(Training_Set_Rate, SNR, Data_Samples);

Input_Layer_Size = size(XTrain_RSRP, [1, 2, 3]);

%Y_train=reshape(cat(1,YTrain_RSRP(:,:,1,:),YTrain_RSRP(:,:,2,:)),2016,[]);
%X_train=reshape(cat(1,XTrain_RSRP(:,:,1,:),XTrain_RSRP(:,:,2,:)),96,[]);
%Y_val=reshape(cat(1,YValidation_RSRP(:,:,1,:),YValidation_RSRP(:,:,2,:)),2016,[]);
%X_val=reshape(cat(1,XValidation_RSRP(:,:,1,:),XValidation_RSRP(:,:,2,:)),96,[]);

% % reshaped_Xtrain=reshape(XTrain_RSRP,96,[]);
% % reshaped_Xval=reshape(XValidation_RSRP,96,[]);
% % 
% % reshaped_Ytrain_Pilots=reshape(YTrain_Pilots_RSRP,96,[]);
% % reshaped_Yval_Pilots=reshape(YValidation_Pilots_RSRP,96,[]);
% % 
% % reshaped_Ytrain_full_Frame=reshape(YTrain_RSRP,72*14*2,[]);
% % reshaped_Yval_Full_Frame=reshape(YValidation_RSRP,72*14*2,[]);

SL102
%Improved_ReEsNet
%Improved_Neural_Network
%ReEsNet

% Option settings
Options = trainingOptions('adam', ...
    'MaxEpochs',100, ...
    'MiniBatchSize',128, ...
    'InitialLearnRate',1e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',20, ...
    'ValidationData',{XValidation_RSRP,YValidation_RSRP}, ...
    'ValidationFrequency',ValidationFrequency, ...
    'Shuffle','every-epoch', ...
    'Verbose',1, ...
    'L2Regularization',0.0000000001, ...
    'ExecutionEnvironment','auto', ...%'parallel'
    'Plots','training-progress');

% Train Network
[SL_EVA, info] = trainNetwork(XTrain_RSRP, YTrain_RSRP, lgraph, Options);

% Network Pruning
%CDF_Layerweights
