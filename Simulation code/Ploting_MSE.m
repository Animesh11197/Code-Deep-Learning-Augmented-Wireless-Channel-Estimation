compiled_MSE=[MSE_LS_over_SNR MSE_MMSE_over_SNR MSE_DNN_over_SNR MSE_ReEsNet_over_SNR MSE_Improved_DNN_over_SNR];
figure;
semilogy(SNR_Range,MSE_LS_over_SNR,'b');
hold on
semilogy(SNR_Range,MSE_MMSE_over_SNR,'r');
hold on
semilogy(SNR_Range,MSE_LS_DNN_concat_over_SNR,'g');
hold on
semilogy(SNR_Range,MSE_ReEsNet_over_SNR,'*--');
hold on
semilogy(SNR_Range,MSE_Improved_DNN_over_SNR_2,'y');
hold on
semilogy(SNR_Range,MSE_Improved_DNN_over_SNR,'c');
hold on
semilogy(SNR_Range,MSE_Improved_DNN_over_SNR_3,'k');
hold on
semilogy(SNR_Range,MSE_LS_DNN_over_SNR,'b--');

legend('LS', ...
    'MMSE', ...
    'LS DNN test', ...
    'ResNet B', ...
    'Interpolated ResNet with 3 blocks', ...
    'Interpolated ResNet', ...
    'Interpolation ReEsNet with 2 blocks', ...
    'LS DNN');
xlabel('SNR in dB');
ylabel('MSE');
title('MSE 48 pilot');
grid on;
hold off;
