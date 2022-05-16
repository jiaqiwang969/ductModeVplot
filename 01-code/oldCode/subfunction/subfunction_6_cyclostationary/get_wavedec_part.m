
function s=get_wavedec_part(x,scale,wname,which_part)

if nargin<4
    which_part='d';
end

[col row]=size(x);
if col>row
    x=x';    
end

[c,l]=wavedec(x,scale,wname);
a=appcoef(c,l,wname);

if which_part=='d'   
    cnew=zeros(1,length(a));
    for i=scale:-1:1
        tmpd=detcoef(c,l,i);
        cnew=[cnew tmpd];
    end
else
    cnew=a;
    for i=scale:-1:1
        tmpd=detcoef(c,l,i);
        zerod=zeros(1,length(tmpd));
        cnew=[cnew zerod];
    end
end
    
s=waverec(cnew,l,wname);

if nargout==0
    figure(gcf+1);
    plot(s);
    figure(gcf+1);
    plot(x);
end