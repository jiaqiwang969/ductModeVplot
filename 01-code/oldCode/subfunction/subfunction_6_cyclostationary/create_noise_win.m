function Noisesignal=create_noise_win(type,amp,M,passrad,N)

if nargin==0
    type='whole';
end

randn('state',0);
Noisesignal=amp*randn(1,N); 
switch type 
    case 'whole'
        Noisesignal=Noisesignal; 
    case 'low'
        w1=boxcar(M);
        n=1:M;
        hd=sin(passrad*(n-ceil(M/2)))./(pi*(n-ceil(M/2)));
        hd(41)=passrad/pi;
        h1=hd.*rot90(w1);
        Noisesignal=filter_win(Noisesignal,h1);
    case 'high'
        w1=boxcar(M);
        n=1:M; 
        hd=(sin(pi*(n-ceil(M/2)))-sin(passrad*(n-ceil(M/2))))./(pi*(n-ceil(M/2)));
        hd(ceil(M/2))=(pi-passrad)/pi;
        h1=hd.*rot90(w1);
        Noisesignal=filter_win(Noisesignal,h1);
    case 'strip_pass'
        w1=boxcar(M);
        n=1:M; 
        hd=(sin(passrad(1)*(n-ceil(M/2)))-sin(passrad(2)*(n-ceil(M/2))))./(pi*(n-ceil(M/2)));
        hd(ceil(M/2))=(passrad(1)-passrad(2))/pi;
        h1=hd.*rot90(w1);
        Noisesignal=filter_win(Noisesignal,h1);
    case 'strip_stop'
        w1=boxcar(M);
        n=1:M; 
        hd=(sin(passrad(1)*(n-ceil(M/2)))+sin(pi*(n-ceil(M/2)))-sin(passrad(2)*(n-ceil(M/2))))./(pi*(n-ceil(M/2)));
        hd(ceil(M/2))=(passrad(1)+pi-passrad(2))/pi;
        h1=hd.*rot90(w1);
        Noisesignal=filter_win(Noisesignal,h1);
end