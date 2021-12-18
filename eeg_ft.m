eegsample = 1000;
adjusted_desired_phase = 0;
phase_tolereance = 0.1;
samplerate = 500;
downsample = eegsample/samplerate;
window = 1000;
ps_window = 900;
N = 0:window-ps_window-1;
edge = 41;
forwardsample = 64;
sample = 0;
fs = 500;
n = 500;
%hilsample = 0;
p = 32;
phase_idx = 1;
idx = downsample;
x0 = (1:window)*(1/eegsample);
x1 = (ps_window:window)*(1/eegsample);
x2 = (1:window)*(1/eegsample);
x3 = (1:ps_window)*(1/eegsample);
allvec = nan(1,1000000);
allts = nan(1,1000000);
phase = nan(1000,forwardsample);
op_phase = nan(1000,window);
real_phase = nan(1000,ps_window);
predict_phase = nan(1000,window-ps_window);
concate_phase = nan(1000,window);
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
    if downsample == idx
        idx = 1;
        sample = sample + 1;
        allvec(1,sample) = vec(1);
        allts(1,sample) = ts;
        if mod(sample,window) == 0
            disp('sample: ')
            fprintf('%.2f\t\n',sample);
            chunk = allvec(1,sample-window+1:sample);
            phase_error_chunk = allvec(1,sample-window+1:sample-window+ps_window);
            Xf = fft(phase_error_chunk,4096);
            true_phase = angle(hilbert(phase_error_chunk));
            phase_now = true_phase(ps_window);
            [~,mad] = max(abs(Xf));
            
            
            
            f_est = (mad-1)*fs/length(Xf);
            y = cos(2*pi*f_est*N/fs+phase_now);
            
%             phase_est = angle(Xf(mad)); 
%             phase = mod(2*pi*f_est*(phase_error_window-1)/fs+phase_est,2*pi);
%             phase = wrapToPi(phase);
%             op_phase(phase_idx,:) = angle(hilbert(chunk));
%             real_phase(phase_idx,:) = angle(hilbert(phase_error_chunk));
%             predict_phase(phase_idx,:) = angle(hilbert(y));
            op_phase(phase_idx,:) = chunk;
            real_phase(phase_idx,:) = phase_error_chunk;
            predict_phase(phase_idx,:) = y;
            concate_phase(phase_idx,:) = [real_phase(phase_idx,:),predict_phase(phase_idx,:)];
            phase_idx = phase_idx + 1;
        end
        
    else
        idx = idx + 1;
    end

    if sample == 9000
        break
    end
    
end

% concate_phase = [real_phase,phase];



figure(10);
for i = 1:9
    subplot(3,3,i);
    hold on;
%     plot(x1,concate_phase(i,:),'r');
    plot(x0,concate_phase(i,:),'r');
    plot(x2,op_phase(i,:),'b');
    xline(numel(real_phase(i,:))/eegsample,'k--',{'now'});
    grid on;
%     t = concate_phase(i,ps_window:window);
%     p1 = find(abs(concate_phase(i,ps_window:window)-pi) < 0.15);
%     plot(x1(p1),t(p1),'o','color','g')
%     
%     p2 = find(abs(concate_phase(i,ps_window:window)-0) < 0.15);
%     plot(x1(p2),t(p2),'o','color','r')
    xlabel('time(sec)');
    ylabel('amplitude');
    hold off;
end






