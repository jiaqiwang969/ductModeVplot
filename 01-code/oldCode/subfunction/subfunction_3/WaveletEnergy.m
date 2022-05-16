function A=WaveletEnergy(x,N,Type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%x为输入的一维信号，n为小波包的层数
    
    T = WPDEC(x,N,Type);  %小波包分解
    M=2^N;
    for j=0:(M-1)
        B=WPCOEF(T,[N,j]);
        A(j+1)=sum(B.^2)/length(B);
    end
    A=A/sum(A);
end