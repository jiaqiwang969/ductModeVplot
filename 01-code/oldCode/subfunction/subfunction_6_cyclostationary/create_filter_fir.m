function b=create_filter_fir(n,Wn,ftype,win)

if nargin==2
     ftype='low';
     win='hamming';
 elseif nargin==3
     win='hamming';
 elseif nargin<2||nargin>4
     error('at least two inputs and at most four inputs.');
 end

if ismember(ftype,{'low','high','stop','bandpass'})
else
    error('the string represesting the type of the filter should be low, high,stop or bandpass.');
end

if strcmp(win,'kaiser')
    para=60;
    beta=0.1102*(para-8.7);
    window=kaiser(n+1,beta);
elseif strcmp(win,'chebwin');
    para=40;
    window=chebwin(n+1,para);
else
    window=eval(strcat(win,'(n+1)'));
end

b=fir1(n,Wn,ftype,window);

if nargout==0
    freqz(b,1,512,2);
end

