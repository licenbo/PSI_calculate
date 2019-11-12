function [phase_synchrony_indices] = calculatePSI(f_min,f_max,data_channel1,data_channel2,fs)
%% test parameter
% load S001R01.mat;
% data = eeg(1:64,:);
% lpass=0.5;
% hpass=50;
% filterorder=3;
% fs=160;
% filtercutoff = [2*lpass/fs 2*hpass/fs]; 
% [f_b, f_a] = butter(filterorder,filtercutoff);
% filteddata=[];
% for k=1:size(data,1)
%     filteddata(k,:)=filtfilt(f_b, f_a,data(k,:));
% end
% data_to_calculatePSI = filteddata(:,3201:3360);
% % C3导联为第9行数据 C4导联为第13行数据
% data_channel1 = data_to_calculatePSI(9,:);
% data_channel2 = data_to_calculatePSI(13,:);
% f_max = 20; f_min = 1;
%%
morlet_wavelet = [];
fi_delta=[];
realvalue_fi_delta=[];
imagevalue_fi_delta=[];
indices_no_mean_real=[];
phase_synchrony_indices_no_mean=[];
fi_Channel1=[];
fi_Channel2=[];
frequency = f_max-f_min+1;
for f = f_min:f_max
    for t = 1:fs
        delta_t = (7/(2*pi*f))*1000;
        delte_f = f/7;
        morlet_wavelet(f,t) = ((delta_t*((pi)^(1/2)))^(-1/2))*(exp(-1*(t^2)/(2*((delta_t)^2))))*(exp(1i*2*pi*f*t));
    end
end
fi_Channel1_conv2 = conv2(morlet_wavelet,data_channel1,'same');
fi_Channel1_abs = abs(conv2(morlet_wavelet,data_channel1,'same'));
fi_Channel2_conv2 = conv2(morlet_wavelet,data_channel2,'same');
fi_Channel2_abs = abs(conv2(morlet_wavelet,data_channel2,'same'));
for x = 1:frequency
    for y = 1:fs
        fi_Channel1(x,y) =  fi_Channel1_conv2(x,y)/fi_Channel1_abs(x,y);
        fi_Channel2(x,y) =  fi_Channel2_conv2(x,y)/fi_Channel2_abs(x,y);
    end
end
for j = 1:frequency
    for k =1:fs
    fi_delta(j,k) = fi_Channel1(j,k)/fi_Channel2(j,k);
    end
end
for a =1:frequency
    for b=1:fs
        realvalue_fi_delta(a,b) = real(fi_delta(a,b));
        imagevalue_fi_delta(a,b) = imag(fi_delta(a,b));
    end
end
realvalue_fi_delta_sum = sum(realvalue_fi_delta,2);
imagevalue_fi_delta_sum = sum(imagevalue_fi_delta,2);
for m = 1:frequency
    indices_no_mean_real(m,1) = (realvalue_fi_delta_sum(m,1))^2;
    indices_no_mean_imag(m,1) = (imagevalue_fi_delta_sum(m,1))^2;
end
for n = 1:frequency
     phase_synchrony_indices_no_mean(n,1) = (sqrt(indices_no_mean_real(n,1)+indices_no_mean_imag(n,1)))/fs;
end
phase_synchrony_indices = mean(phase_synchrony_indices_no_mean,1);
end