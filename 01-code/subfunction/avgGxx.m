%This function computes the the averages autpower of given time domain. The
%program accepts different windows (hann, flat, and rect), different levels
%of overlap, and different correction factors. The program outputs the
%averaged autopower (Gxx),all of the linear spectra used (Gx), and an
%frequency vector (Fx).
%
% Example: 
% [Gx,Gxx,Fx] = fftPlotData('hann',50,'ACF',10,3200,3200,timeData)
% semilogy(Fx,Gxx)
function [Gxx,Gx,Fx] = avgGxx(wind,PercOverlap,CF,numAvgs,Fs,N,timeData)
dF=Fs/N; %Compute deltaF
PercOverlap=PercOverlap/100;
%Create array with columns = number of averages
reshapedO=zeros(N,numAvgs); %Initialize arrays to make for loop more memory efficient.
for ii=1:numAvgs
    reshapedO(:,ii)=timeData((((1-PercOverlap)*(ii-1)*N)+1):((((1-PercOverlap)*(ii-1)*N+N)))); %Create matrix of data to average
end
%reshapedO=reshape(voltageO,[N,numAvgs]); no overlap easy way

%Pick window
if strcmp(wind,'hann')
    windowToUse=window(@hann,N);
end
if strcmp(wind,'flat')
    windowToUse=window(@flattopwin,N);
end
if strcmp(wind,'rect')
    windowToUse=window(@rectwin,N);
end
windowToUse=repmat(windowToUse,1,numAvgs); %Make a matrix of windows that is numAvgs long to window each time history
windowedO=reshapedO.*windowToUse; %Window in time domain

%FFT and add amplitude correction to get linear spectra
if strcmp(CF,'ACF')
    Gx=(1/mean(windowToUse(:,1))*fft(windowedO)/N); %First part is the amplitude correction factor times fft
end
if strcmp(CF,'ECF')
    Gx=(1/rms(windowToUse(:,1))*fft(windowedO)/N); %First part is the amplitude correction factor times fft
end

%only look at positive values
Gx=Gx(1:floor(N/2)+1,:);
%Only multiply frequencies above zero by2.
Gx(2:end,:)=Gx(2:end,:).*2;

%Do averaging on autopowers and crosspowers
for ii=1:numAvgs
    if ii==1
        Gxx=conj(Gx(:,ii)).*Gx(:,ii);
    else
        Gxx=conj(Gx(:,ii)).*Gx(:,ii)+Gxx;
    end
end
Gxx=Gxx/ii;

Fx = (0:1:(length(Gxx)-1))*dF; %Create x axis labels
end