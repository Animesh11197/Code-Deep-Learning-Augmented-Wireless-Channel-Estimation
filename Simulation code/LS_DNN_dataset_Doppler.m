
% for SNR=0:5:20 
    Data_Samples = 5000;
    Training_Set_Rate = 0.995;
    ValidationFrequency = 300;
    SNR = 0:5:20;
    MaxDopplerShift=97;
    
    [XTrain_RSRP, YTrain_RSRP, XValidation_RSRP, YValidation_RSRP] = Data_Generation_LS_DNN(Training_Set_Rate, SNR, Data_Samples,MaxDopplerShift);
    
    new_Xtrain=reshape(cat(1,XTrain_RSRP(:,:,1,:),XTrain_RSRP(:,:,2,:)),96,[]);
    new_Ytrain=reshape(cat(1,YTrain_RSRP(:,:,1,:),YTrain_RSRP(:,:,2,:)),72*14*2,[]);
    new_Xval=reshape(cat(1,XValidation_RSRP(:,:,1,:),XValidation_RSRP(:,:,2,:)),96,[]);
    new_Yval=reshape(cat(1,YValidation_RSRP(:,:,1,:),YValidation_RSRP(:,:,2,:)),72*14*2,[]);


    
    save(['./LS_DNN_dataset/' 'ETU_new_LS_DNN_SNR_' '0to20' '_Doppler_97.mat'],'new_Xtrain','new_Xval','new_Ytrain','new_Yval');
% end