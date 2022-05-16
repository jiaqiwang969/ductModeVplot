function y=Envelop(x)


y=abs(hilbert(x));

if nargout==0
    figure( gcf + 1 );
    plot(y);
    title('Envelop Wave');
    xlabel('Time [s]');
    ylabel('Amplitude');
end