function h=create_filter_win(N,M,sign,type,passrad,para);

if nargin<6
    para=[];
end
n=1:N;

if isempty(para)
    filter_str=[sign '(' num2str(N) ')'];
else
    filter_str=[sign '(' num2str(N) ',' num2str(para) ')'];
end
w=eval(filter_str);                                                 

switch type
    case 'low'
        hd=sin(passrad*(n-ceil(N/2)))./(pi*(n-ceil(N/2)));
        hd(ceil(N/2))=passrad/pi;
    case 'high'
        hd=(sin(pi*(n-ceil(N/2)))-sin(passrad*(n-ceil(N/2))))./(pi*(n-ceil(N/2)));
        hd(ceil(N/2))=(pi-passrad)/pi;
    case 'strip_pass'
        hd=(sin(passrad(1)*(n-ceil(N/2)))-sin(passrad(2)*(n-ceil(N/2))))./(pi*(n-ceil(N/2)));
        hd(ceil(N/2))=(passrad(1)-passrad(2))/pi;
    case 'strip_stop'
        hd=(sin(passrad(1)*(n-ceil(N/2)))+sin(pi*(n-ceil(N/2)))-sin(passrad(2)*(n-ceil(N/2))))./(pi*(n-ceil(N/2)));
        hd(ceil(N/2))=(passrad(1)+pi-passrad(2))/pi;
end

h=hd.*rot90(w);

if nargout==0
    [mag,rad]=freqz(h,1,M);
    figure(gcf+1);
    plot(rad/2/pi,20*log10(abs(mag)));
    ylabel('dB');
    xlabel('*fs[Hz]');
    grid on;
    title(sign);
    figure(gcf+1);
    if fix(N/2)*2==N
        h=[h 0];
    end
    plot([-1*fix(N/2):fix(N/2)],h,'bo-');
    title('time domain');
end