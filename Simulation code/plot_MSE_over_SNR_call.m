 figure;
semilogy(SNR_Range,SW_over_SNR,'c*');
hold on
semilogy(SNR_Range,half_precision_over_SNR,'-r');
hold on
semilogy(SNR_Range,fixed_12_5,'*-g');
hold on
semilogy(SNR_Range,fixed_8_5,'*-y');
hold on
semilogy(SNR_Range,fixed_16_4,'b--o');
legend('SW over SNR', ...
'half precision over SNR', ...
'fixed 12,5', ...
'fixed 8,5', ...
'fixed 16,4');
xlabel('SNR in dB');
ylabel('MSE');
title('Interpolated-ResNet Performance:');
grid on;
hold off