eegsample = 1000; %取樣
samplerate = 500; %downsample頻率
downsample = eegsample/samplerate; %每兩個sample點只取一個
window = 700; %每次一個chunk長度
phase_error_window = 524; %用來預測的長度
edge = 41;
forwardsample = 256; %預測點
sample = 0;
p = 32; %ar參數
fpass=[8,13];%passband

phase_idx = 1;
idx = downsample;
%以下是用來畫圖的
count0 = 0;
count1 = 1;
count2 = 1;

x0 = (1:phase_error_window-edge+forwardsample)*(1/samplerate);
x1 = (phase_error_window-edge:phase_error_window-edge+forwardsample)*(1/samplerate);
x2 = (1:window)*(1/samplerate);

x3 = 1:15;
%
allvec = nan(1,1000000);%用來存lsl的
allts = nan(1,1000000);%用來存lsl的
phase = nan(1000,forwardsample);%存預測資料
%error_forward_phase = nan(1000,phase_error_foward);
op_phase = nan(1000,window);%存真實資料
real_phase = nan(1000,phase_error_window-edge);
true_phase = nan(1000,window);%存真實相位
%error_real_phase = nan(1000,phase_error_window-edge);
% error_phase = nan(1000,phase_error_window-edge+forwardsample);

concate_phase = nan(1000,window-edge+forwardsample);
disp('Loading the library...');
lib = lsl_loadlib();
disp('Resolving an EEG stream...');
result = {}; 

while isempty(result)
    result = lsl_resolve_byprop(lib,'type','EEG1');
end

disp('Opening an inlet...');
inlet = lsl_inlet(result{1});


disp('Now receiving data...');

while 1
    [vec,ts] = inlet.pull_sample();
  
        
    if downsample == idx %downsample
        idx = 1; %當idx=2才會存 所以重新設1
        sample = sample + 1;
        allvec(1,sample) = vec(1);
        allts(1,sample) = ts;
        if mod(sample,window) == 0 %資料數夠一個window
            count0 = count0 + 1;
            disp('sample: ')
            fprintf('%.2f\t\n',sample);
            chunk = allvec(1,sample-window+1:sample);
            chunk = bandpass(chunk,fpass,samplerate);
            phase_error_chunk = allvec(1,sample-window+1:sample-window+phase_error_window);%用來預測的chunk 長度比window短
            phase_error_chunk = bandpass(phase_error_chunk,fpass,samplerate);
            error_coeffs = aryule(phase_error_chunk(edge+1:end-edge), p); 
            error_coeffs = -error_coeffs(end:-1:2);
            error_nextvalues = zeros(1,p+forwardsample);
            error_nextvalues(1:p) = phase_error_chunk(end-p-edge+1:end-edge);
            for i = 1:forwardsample
                error_nextvalues(p+i) = error_coeffs*error_nextvalues(i:p+i-1)'; %預測出的資料
            end
%             phase(phase_idx,:) = angle(hilbert(error_nextvalues(p+1:end)));
%             real_phase(phase_idx,:) = angle(hilbert(phase_error_chunk(1:end-edge)));
%             true_phase(phase_idx,:) = angle(hilbert(chunk));
            phase(phase_idx,:) = error_nextvalues(p+1:end);%預測出的資料
            real_phase(phase_idx,:) = phase_error_chunk(1:end-edge);%用來預測的資料
            true_phase(phase_idx,:) = chunk;%原始資料
            
            phase_idx = phase_idx + 1;
                
        end
            
    
        
        else
            idx = idx + 1;
    end

    if sample == 6300
        n_trial=sample/window;
        break
    end
    
end

concate_phase = [real_phase,phase];

figure(10);
for i = 1:9
    subplot(3,3,i);
    hold on;
%     plot(x1,concate_phase(i,:),'r');
    plot(x0,concate_phase(i,:),'r');
    plot(x2,true_phase(i,:),'b');
    xline(numel(real_phase(i,:))/eegsample,'k--',{'now'});
    
    xlabel('time(sec)');
    ylabel('amplitude');
    grid on;
%     y = concate_phase(i,phase_error_window-edge:phase_error_window-edge+forwardsample);
%     p1 = find(abs(concate_phase(i,phase_error_window-edge:phase_error_window-edge+forwardsample)-pi) < 0.15);
%     plot(x1(p1),y(p1),'o','color','g')
%     
%     p2 = find(abs(concate_phase(i,phase_error_window-edge:phase_error_window-edge+forwardsample)-0) < 0.15);
%     plot(x1(p2),y(p2),'o','color','r')
    
    hold off;
end

% for x = 1:1
%     figure(x)
%     plot(x2,op_phase(x,:),'b');
%     hold;
% end
%% plv heatmap fft
now_time=0.8;
max_future_time=0.5;
past_windows=0.02:0.02:now_time;%in secs
future_windows=0.02:0.02:max_future_time;%in secs
PLVmap=nan(numel(past_windows),numel(future_windows));
idx1=1;
idx2=1;
for past_window=past_windows
    idx2=1;
    for future_window=future_windows
        past_samples=round(past_window*samplerate);
        future_samples=round(future_window*samplerate);
        pastEEG=true_phase(1:n_trial,now_time*samplerate-past_samples+1:now_time*samplerate);
        real_futureEEG=true_phase(1:n_trial,now_time*samplerate+1:...
            now_time*samplerate+future_samples);
        %find the estimate of parameters(peak frequency and current phase)
        L=4096;
        pastFFT=fft(pastEEG,L,2);
        [~,fidx] = max(abs(pastFFT),[],2);
        peak_freq=samplerate*(fidx-1)/L;
        now_phase=angle(hilbert(pastEEG'))';
        now_phase=now_phase(:,end);
        %make prediction using cosine wave
        predicted_futureEEG=nan(n_trial,future_samples);
        for i=1:n_trial
            predicted_futureEEG(i,:)=cos(2*pi*peak_freq(i)*(1:future_samples)+now_phase(i));
        end
        fprintf('predicted_futureEEG:');
        print_arr(size(predicted_futureEEG));
        fprintf('real_futureEEG:');
        print_arr(size(real_futureEEG));
        
        filtSpec.order = 50;
        filtSpec.range = fpass; %Hz
        plv=0;
        for trial=1:n_trial
            plv=plv+ my_eegPLV(predicted_futureEEG(trial,:),real_futureEEG(trial,:), samplerate, filtSpec);
        end
        plv=plv/n_trial;
        PLVmap(idx1,idx2)=plv;
        idx2=idx2+1;
    end
    idx1=idx1+1;
end
xvalues = future_windows;
yvalues = past_windows;
figure;
h = heatmap(xvalues,yvalues,PLVmap);

h.Title = 'PLV in alpha band';
h.XLabel = 'future window(sec)';
h.YLabel = 'past window(sec)';
%% plv heatmap hilbert
now_time=0.8;
max_future_time=0.5;
ar_windows=0.02:0.02:0.3;%in secs;ar_order/samplerate
future_windows=0.02:0.02:max_future_time;%in secs
PLVmap=nan(numel(ar_windows),numel(future_windows));
idx1=1;
idx2=1;
for past_window=ar_windows
    idx2=1;
    for future_window=future_windows
        p=round(past_window*samplerate);%ar order
        future_samples=round(future_window*samplerate);
        pastEEG=true_phase(1:n_trial,1:now_time*samplerate);
        real_futureEEG=true_phase(1:n_trial,now_time*samplerate+1:...
            now_time*samplerate+future_samples);
        %ar prediction
        error_coeffs = aryule(pastEEG', p); 
        error_coeffs = -error_coeffs(1:n_trial,end:-1:2);
        predicted_futureEEG = zeros(n_trial,p+future_samples);
        predicted_futureEEG(1:n_trial,1:p) = pastEEG(1:n_trial,end-p+1:end);
        for i = 1:future_samples
            %預測出的資料%(n_trial x p)(n_trial x p)
            predicted_futureEEG(1:n_trial,p+i) =sum( error_coeffs.*...
                predicted_futureEEG(1:n_trial,i:p+i-1),2); 
        end
        predicted_futureEEG=predicted_futureEEG(:,p+1:end);
        fprintf('predicted_futureEEG:');
        print_arr(size(predicted_futureEEG));
        fprintf('real_futureEEG:');
        print_arr(size(real_futureEEG));
        
        filtSpec.order = 50;
        filtSpec.range = fpass; %Hz
        plv=0;
        for trial=1:n_trial
            plv=plv+ my_eegPLV(predicted_futureEEG(trial,:),real_futureEEG(trial,:), samplerate, filtSpec);
        end
        plv=plv/n_trial;
        PLVmap(idx1,idx2)=plv;
        idx2=idx2+1;
    end
    idx1=idx1+1;
end
xvalues = future_windows;
yvalues = ar_windows;
figure;
h = heatmap(xvalues,yvalues,PLVmap);

h.Title = 'PLV in alpha band';
h.XLabel = 'future window(sec)';
h.YLabel = 'ar window(sec)';     
        
