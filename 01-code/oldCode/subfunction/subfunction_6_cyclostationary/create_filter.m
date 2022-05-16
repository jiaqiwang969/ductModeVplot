function [b,a,n]=create_filter(filter_class,Wp,Ws,Rp,Rs)

if nargin~=5
    error('five inputs');
end

if strcmp(filter_class,'butt')
    tysign={'butter','buttord'};
    paramstr='(n,Wn)';
elseif strcmp(filter_class,'cheb1')
    tysign={'cheby1','cheb1ord'};
    paramstr='(n,Rp,Wn)';
elseif strcmp(filter_class,'cheb2')
    tysign={'cheby2','cheb2ord'};
    paramstr='(n,Rs,Wn)';
elseif strcmp(filter_class,'ellip')
    tysign={'ellip','ellipord'};
    paramstr='(n,Rp,Rs,Wn)';
else
    tysign={'butter','buttord'};
    paramstr='(n,Wn)';
end

[n,Wn]=eval(strcat(tysign{2},'(Wp,Ws,Rp,Rs);'));
[b,a]=eval(strcat(tysign{1},paramstr));

if nargout==0
    figure;
    set(gcf,'color','white');
    impz(b,a);
    figure;
    set(gcf,'color','white');
    freqz(b,a);
end



