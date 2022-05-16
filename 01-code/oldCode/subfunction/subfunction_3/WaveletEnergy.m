function A=WaveletEnergy(x,N,Type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%xΪ�����һά�źţ�nΪС�����Ĳ���
    
    T = WPDEC(x,N,Type);  %С�����ֽ�
    M=2^N;
    for j=0:(M-1)
        B=WPCOEF(T,[N,j]);
        A(j+1)=sum(B.^2)/length(B);
    end
    A=A/sum(A);
end