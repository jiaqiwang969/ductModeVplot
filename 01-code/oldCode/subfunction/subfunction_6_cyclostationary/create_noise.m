function noise=create_noise(tyrand,N,W);

if (W>1/1.28)
    error('cutoff frequency is 1/1.28 ');
end
Wdiff=diff(W);

if isempty(Wdiff)
    Wdiff=W;
    Wp=W-Wdiff/20;
    Ws=W+Wdiff/20;
else
    Wp=[W(1)+Wdiff/20,W(2)-Wdiff/20];
    Ws=[W(1)-Wdiff/20,W(2)+Wdiff/20];
end


noise=eval(strcat(tyrand,'(2^nextpow2(N),1)-0.5'));
[b,a]=create_filter('ellip',Wp,Ws,3,40);
noise=filter(b,a,noise);
noise=noise/std(noise);
noise=noise(length(noise)-(N-1:-1:0));
noise=noise-mean(noise);


if nargout==0
    draw_x(noise);
    showfft(noise);
end
