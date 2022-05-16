function Smax=part_integral_s(Suse,cutsign)

if nargin==1
    cutsign=0;
end

if cutsign==0
    Smax=max(Suse);
else
    Smax=max(Suse');
end

if nargout==0  
    [r,c]=size(Smax);
    if floor(c/2)*2==c
        Smax_cl=floor(c/2)+1:c;
    else
        Smax_cl=floor(c/2)+2:c;
    end
    dispSmax=Smax(Smax_cl);
    xl=((0:floor(c/2)-1))/c;
    figure(gcf+1);
    plot(xl,abs(dispSmax));
    xlabel('*2pi ÆµÂÊ');
    
end