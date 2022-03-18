clc;clear all;close all
%% SIGNAL ACQUISITION ==
rec1 = importdata('rec_1m.mat');
unfiltecg = rec1(1,:);

%% == TUNABLE PARAMETERS ==
IMPROVE=0.5; %% CAMBIANDO ESTA DETECTA PICOS EN SEÑALES MAS RUIDOSAS
sq = 0;

%% Signal chacracteristics
Fs = 512;
L = length(unfiltecg);
% Time vector
time = 0:(1/Fs):(L/Fs)-(1/Fs); 
unfiltecg = unfiltecg - mean(unfiltecg);


%% 
figure('Name','Raw')
subplot(211)
plot(time,unfiltecg);
title('Raw Signal');
xlabel('Time(s)')
ylabel('Amplitude')
% Frequency domain
sp= fft(unfiltecg);norm = abs(sp)/L;half= norm(1:(L/2)+1);
% Frequency vector
frequency = Fs.*(0:(L/2))/L;
subplot(212)
plot(frequency,half)
title ('Original spectre')
xlabel('Frequency (Hz)')
ylabel('Magnitude of FT')

%% High Pass Filter

[b,a] = cheby2(4,20,0.6/(Fs/2),'high'); %% cheby2(order,stopband ripple,stopband-edge frequency)
figure 
freqz(b,a); 
HPecg = filtfilt(b,a,unfiltecg);
figure ('Name','Highpass Filtered')
subplot(211)
plot(time,HPecg)
title('High Pass filtered signal')
xlabel('Time(s)')
ylabel('Amplitude')

sp2= fft(HPecg);norm2 = abs(sp2)/L;half2= norm2(1:(L/2)+1);frequency2 = Fs.*(0:(L/2))/L;
subplot(212)
plot(frequency,half2);
xlabel('Frequency (Hz)')
ylabel('Magnitude of FT')

%% Average Filter

avgFiltB = ones(1,8)/8; %% Window design of 8-points
MAecg = filter(avgFiltB,1,HPecg);
figure ('Name','Average Filter')
plot(time,MAecg)
title ('Average filter')
xlabel('Time(s)')
ylabel('Amplitude')

%% Notch filter

BandStop = designfilt('bandstopiir', ...       % Response type
       'PassbandFrequency1',40, ...    % Frequency constraints
       'StopbandFrequency1',50, ...
       'StopbandFrequency2',60, ...
       'PassbandFrequency2',70, ...
       'PassbandRipple1',1, ...         % Magnitude constraints
       'StopbandAttenuation',55, ...
       'PassbandRipple2',1, ...
       'DesignMethod','ellip', ...      % Design method
       'MatchExactly','both', ...       % Design method options
       'SampleRate',Fs); 
Bsdata = filter(BandStop,MAecg);

figure('Name','Notch filter');
subplot(211)
plot(time,Bsdata);
xlabel('Time(s)')
ylabel('Amplitude')
%in frequency
Bs1 = fft(Bsdata);Ba1 = abs(Bs1)/L;BaS1 = Ba1(1:(L/2)+1);f = Fs.*(0:(L/2))/L;
subplot(212)
plot(f,BaS1);
xlabel('Frequency (Hz)')
ylabel('Magnitude of FT')

%% Bandpass 5 to 15

BP = designfilt('bandpassiir', ...       % Response type
       'StopbandFrequency1',2, ...    % Frequency constraints
       'PassbandFrequency1',5, ...
       'PassbandFrequency2',15, ...
       'StopbandFrequency2',18, ...
       'StopbandAttenuation1',40, ...   % Magnitude constraints
       'PassbandRipple',1, ...
       'StopbandAttenuation2',60, ...
       'DesignMethod','ellip', ...      % Design method
       'MatchExactly','passband', ...   % Design method options
       'SampleRate',Fs); 
BPdata = filter(BP,Bsdata);

figure('Name','Bandpass 5-15 Hz');
subplot(211)
%% Squaring
if sq== 1
    BPdata = BPdata.^2;
end
%%
plot(time,BPdata);
xlabel('Time(s)')
ylabel('Amplitude')
%in frequency
BP1 = fft(BPdata);BPa1 = abs(BP1)/L;BPS1 = BPa1(1:(L/2)+1);f = Fs.*(0:(L/2))/L;
subplot(212)
plot(f,BPS1);
xlabel('Frequency (Hz)')
ylabel('Magnitude of FT')


%% SIGNAL DENOISING END

%% First peak finding

figure("Name",'Findpeaks function')
findpeaks(BPdata,'MinPeakHeight',0.05*100,'MinPeakDistance',round(0.196*Fs));
xlabel('N. samples')
ylabel('Amplitude')
[A B]=findpeaks(BPdata,'MinPeakHeight',0.05*100,'MinPeakDistance',round(0.196*Fs));
title('Findpeaks function')

B = [min(B) B];
A = [min(A) A];

%% == Indixes discrimination ==
maxPeak = 0;
lastPeak=0;
lastRR=0;
RRinterval =0;
finepeaks = [];
finepeakspos=[];
RR = [];
maxPEAk=[];
RRthrs=[];
lastthrs=0;


for i=1: length(B)-1
    if(((A(i+1) - A(i))>= 1)) && ((B(i+1) - lastPeak)>=round(0.6*RRinterval)) && (A(i+1)>= IMPROVE*maxPeak)
        
%         finepeaks(i+1) = B(i+1);

        finepeaks(end+1) = B(i+1);
        finepeakspos(end+1) = (i+1);
        
        RRinterval = B(i+1)-lastPeak;
        RR(end+1) = RRinterval;

        if i>1

        if RRinterval >= round(1.5*RR(end-1))
            aa = (1+finepeakspos(end-1))
            bb = (-1+finepeakspos(end))
            [sa sb] = max(A(aa:bb));

            if A(-1+aa+sb)>= IMPROVE*maxPeak

                finepeaks(end-1) = B(-1+aa+sb)
                finepeaks(end+1)= B(i+1)

                finepeakspos(end-1) = (-1+aa+sb)
                finepeakspos(end+1)= (i+1)

                RRinterval = B(-1+aa+sb)-lastPeak;
                RR(end-1) = RRinterval;
                RR(end+1) = RRinterval;

            end

            

        end
        end

        
        lastthrs = lastRR + RRinterval;
        lastRR = RRinterval;
        


        lastPeak = B(i+1);

        maxPeak =  BPdata(B(i+1)); 
        maxPEAk(end+1) = BPdata(B(i+1));

    end
end
% finepeaks = finepeaks
maxPEAk = maxPEAk.* 0.7;


%%      == Peaks in raw signal ==

for i = 1:length(finepeaks)-2
    [rawamp(i),rawindex(i)]=max(unfiltecg(finepeaks(i)-(round(RR(i)*0.5)):finepeaks(i)));
   
    rawindex(i) = (rawindex(i)-1)+(finepeaks(i)-(round(RR(i)*0.5)));
    
    if i <= length(finepeaks)-1
    shotindex(i) = rawindex(i)+round(0.7*RR(i));
    shotamp(i) = unfiltecg(shotindex(i));
    end
end

    


figure ("Name",'Raw signal w. R n shoot')
h1=stem(rawindex(2:end),rawamp(2:end),'b');
hold on
h2=plot(unfiltecg,'k');
hold on 
h3=xline(shotindex(2:end),'r');
hold on
h4=plot(shotindex(2:end),shotamp(2:end),'*r');
legend([h3],{'Shoot'})
title('Raw signal w. R n shoot')
xlabel('Number of Samples')
ylabel('Amplitude')



