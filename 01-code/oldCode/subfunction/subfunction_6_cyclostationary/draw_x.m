function draw_x(x,df,plot_sign,showstr,len,range)


if nargin<2
    df=1/length(x);
    plot_sign=[1 1];
    showstr=char('','Time (s)','Acceleration (m/s^2)');
    len=length(x)*df;
    range=[0,0];
elseif nargin<3
    plot_sign=[1 1];
    showstr=char('','Time (s)','Acceleration (m/s^2)');
    len=length(x)*df;
    range=[0,0];
elseif nargin<4
    showstr=char('','Time (s)','Acceleration (m/s^2)');
    len=length(x)*df;
    range=[0,0];
elseif nargin<5
    len=length(x)*df;
    range=[0,0];
elseif nargin<6
    range=[0,0];
end

if length(plot_sign)~=2
    error('plot_sign must have two elements.');
elseif plot_sign(2)>plot_sign(1)
    error('subplot wrong');
end

if len>length(x)*df
    len=length(x)*df;
end
[strrow,strcol]=size(showstr);
if strrow~=3
    error('you should give tiltle,xlabel and ylabel at the same time.');
end

showL=ceil(len/df)+1;
if showL>length(x)
    showL=length(x);
end
x=x(1:showL);

xaxismax=df*length(x);

if plot_sign(1)>1
    if plot_sign(2)==1
        figure;


        subplot(plot_sign(1),1,plot_sign(2));
        
    else
        
        subplot(plot_sign(1),1,plot_sign(2));
    end
else
    figure;

end
plot((0:showL-1)*df,x,'k');
set(gcf,'color','white');
set(gca,'color','none');
set(gca,'fontsize',12);
set(get(gca,'xlabel'),'string','\fontsize{12}');
set(get(gca,'ylabel'),'string','\fontsize{12}');
title(deblank(showstr(1,:)));
xlabel(deblank(showstr(2,:)));
ylabel(deblank(showstr(3,:)));

if sum(abs(range))==0
    range=[2*min(x),2*max(x)];
    axis([0,len,range]);
else
    axis([0,len,range]);
end

