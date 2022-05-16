%FigSetting_12000rpm
%µ¼ÈëÒÑ¾­¼ÆËãºÃµÄcpsdÍ¼Ïñ
clc
clear
close all
[fname,location]=uigetfile({'*.fig';'*.*'},'fig²ÎÊıÎÄ¼ş¶ÁÈ¡','MultiSelect','on');%MultiSelectµ¥Ñ¡
% % //========================================  
if isstr(fname)
   fname=cellstr(fname);
end
temp1=strsplit(location,'\');
save_directory=[temp1{end-2},'-',temp1{end-1}];%ÆµÆ×Í¼´æ´¢ÎÄ¼ş¼Ğ
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('ÎÄ¼ş¼Ğ´æÔÚ£¡');
end


for i_file=57:length(fname)
    h=openfig(fullfile(location,char(fname(i_file))),'reuse');ylim([1,30]);
    D1=get(gca,'Children');%get the handle of the line object
    XData1=get(D1,'XData');%get the x data
    YData1=get(D1,'YData');%get the y data;
    ZData1=get(D1,'CData');%get the z data;
    t1=strsplit(char(fname(i_file)),'-');t2=strsplit(t1{5},'rpm');rotor_speed=str2num(t2{1});   
%     if rotor_speed<11000
%         continue;
%     end
    df=YData1(2)-YData1(1);
    [x0,y0]=find(ZData1==max(max(ZData1(round(28.5/df):round(29.5/df),:))),1);Zdata=ZData1(x0,:);Zdata(y0)=0;[y0_1]=find(Zdata==max(Zdata));  
    [x1,y1]=find(ZData1==max(max(ZData1(round(10/df):round(20/df),:))),1);Zdata=ZData1(x1,:);Zdata(y1)=0;[y1_1]=find(Zdata==max(Zdata));
    [x2,y2]=find(ZData1==max(max(ZData1(round(10/df):x1-20,:))),1);Zdata=ZData1(x2,:);Zdata(y2)=0;[y2_1]=find(Zdata==max(Zdata));
    [x3,y3]=find(ZData1==max(max(ZData1(x1+20:round(20/df),:))),1);Zdata=ZData1(x3,:);Zdata(y3)=0;[y3_1]=find(Zdata==max(Zdata));
    o_RI0=round(x0*df*10)/10;m_RI0=y0-17;o_RI1=round(x1*df*10)/10;m_RI1=y1-17;o_RI2=round(x2*df*10)/10;m_RI2=y2-17;o_RI3=round(x3*df*10)/10;m_RI3=y3-17;
    o_RI0_1=round(x0*df*10)/10;m_RI0_1=y0_1-17;o_RI1_1=round(x1*df*10)/10;m_RI1_1=y1_1-17;o_RI2_1=round(x2*df*10)/10;m_RI2_1=y2_1-17;o_RI3_1=round(x3*df*10)/10;m_RI3_1=y3_1-17;
    
    hold on
    text(m_RI0,o_RI0,{[num2str(m_RI0),',',num2str(round(ZData1(x0,y0)))];[num2str(o_RI0)]});text(m_RI0_1,o_RI0_1,[num2str(m_RI0_1),',',num2str(round(ZData1(x0,y0_1)))]);
    text(m_RI1,o_RI1,{[num2str(m_RI1),',',num2str(round(ZData1(x1,y1)))];[num2str(o_RI1)]});text(m_RI1_1,o_RI1_1,[num2str(m_RI1_1),',',num2str(round(ZData1(x1,y1_1)))]);
    text(m_RI2,o_RI2,{[num2str(m_RI2),',',num2str(round(ZData1(x2,y2)))];[num2str(o_RI2)]});text(m_RI2_1,o_RI2_1,[num2str(m_RI2_1),',',num2str(round(ZData1(x2,y2_1)))]);
    text(m_RI3,o_RI3,{[num2str(m_RI3),',',num2str(round(ZData1(x3,y3)))];[num2str(o_RI3)]});text(m_RI3_1,o_RI3_1,[num2str(m_RI3_1),',',num2str(round(ZData1(x3,y3_1)))]);



    
    h1=figure;
    bar(XData1,ZData1(x1,:));hold on
    ylim([70 125]);
    title({[char(fname(i_file))];['f=',num2str(rotor_speed/60*o_RI1),',','Order=',num2str(o_RI1),',m=',num2str(m_RI1)]},'FontSize',14)
    h2=figure;
    bar(XData1,ZData1(x2,:));hold on
    ylim([70 125]);
    title({[char(fname(i_file))];['f=',num2str(rotor_speed/60*o_RI2),',','Order=',num2str(o_RI2),',m=',num2str(m_RI2)]},'FontSize',14)
    h3=figure;
    bar(XData1,ZData1(x3,:));hold on
    ylim([70 125]);
    title({[char(fname(i_file))];['f=',num2str(rotor_speed/60*o_RI3),',','Order=',num2str(o_RI3),',m=',num2str(m_RI3)]},'FontSize',14)

   
%    switch(t1{1})
%    case '±ï' 
%      mkdir([save_directory,'\','1-±ï']) 
%      saveas(h,[save_directory,'\','1-±ï','\',char(fname(i_file))]);
%      saveas(h,[save_directory,'\','1-±ï','\',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h1,[save_directory,'\','1-±ï','\','R1-',char(fname(i_file))]);
%      saveas(h1,[save_directory,'\','1-±ï','\','R1-',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h2,[save_directory,'\','1-±ï','\','R2-',char(fname(i_file))]);
%      saveas(h2,[save_directory,'\','1-±ï','\','R2-',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h3,[save_directory,'\','1-±ï','\','R3-',char(fname(i_file))]);
%      saveas(h3,[save_directory,'\','1-±ï','\','R3-',strrep(char(fname(i_file)),'.fig','.png')]);
%    case '±ï2' 
%      mkdir([save_directory,'\','2-±ï2']) 
%      saveas(h,[save_directory,'\','2-±ï2','\',char(fname(i_file))]);
%      saveas(h,[save_directory,'\','2-±ï2','\',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h1,[save_directory,'\','2-±ï2','\','R1-',char(fname(i_file))]);
%      saveas(h1,[save_directory,'\','2-±ï2','\','R1-',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h2,[save_directory,'\','2-±ï2','\','R2-',char(fname(i_file))]);
%      saveas(h2,[save_directory,'\','2-±ï2','\','R2-',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h3,[save_directory,'\','2-±ï2','\','R3-',char(fname(i_file))]);
%      saveas(h3,[save_directory,'\','2-±ï2','\','R3-',strrep(char(fname(i_file)),'.fig','.png')]);
%    case '±ï3' 
%      mkdir([save_directory,'\','3-±ï3']) 
%      saveas(h,[save_directory,'\','3-±ï3','\',char(fname(i_file))]);
%      saveas(h,[save_directory,'\','3-±ï3','\',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h1,[save_directory,'\','3-±ï3','\','R1-',char(fname(i_file))]);
%      saveas(h1,[save_directory,'\','3-±ï3','\','R1-',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h2,[save_directory,'\','3-±ï3','\','R2-',char(fname(i_file))]);
%      saveas(h2,[save_directory,'\','3-±ï3','\','R2-',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h3,[save_directory,'\','3-±ï3','\','R3-',char(fname(i_file))]);
%      saveas(h3,[save_directory,'\','3-±ï3','\','R3-',strrep(char(fname(i_file)),'.fig','.png')]);
%    case '±ï60'
%      mkdir([save_directory,'\','4-±ï60-50']) 
%      saveas(h,[save_directory,'\','4-±ï60-50','\',char(fname(i_file))]);
%      saveas(h,[save_directory,'\','4-±ï60-50','\',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h1,[save_directory,'\','4-±ï60-50','\','R1-',char(fname(i_file))]);
%      saveas(h1,[save_directory,'\','4-±ï60-50','\','R1-',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h2,[save_directory,'\','4-±ï60-50','\','R2-',char(fname(i_file))]);
%      saveas(h2,[save_directory,'\','4-±ï60-50','\','R2-',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h3,[save_directory,'\','4-±ï60-50','\','R3-',char(fname(i_file))]);
%      saveas(h3,[save_directory,'\','4-±ï60-50','\','R3-',strrep(char(fname(i_file)),'.fig','.png')]);
%    case '±ï50' 
%      mkdir([save_directory,'\','5-±ï50-40']) 
%      saveas(h,[save_directory,'\','5-±ï50-40','\',char(fname(i_file))]);
%      saveas(h,[save_directory,'\','5-±ï50-40','\',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h1,[save_directory,'\','5-±ï50-40','\','R1-',char(fname(i_file))]);
%      saveas(h1,[save_directory,'\','5-±ï50-40','\','R1-',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h2,[save_directory,'\','5-±ï50-40','\','R2-',char(fname(i_file))]);
%      saveas(h2,[save_directory,'\','5-±ï50-40','\','R2-',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h3,[save_directory,'\','5-±ï50-40','\','R3-',char(fname(i_file))]);
%      saveas(h3,[save_directory,'\','5-±ï50-40','\','R3-',strrep(char(fname(i_file)),'.fig','.png')]);
%    case '±ï40' 
%      mkdir([save_directory,'\','6-±ï40-0']) 
%      saveas(h,[save_directory,'\','6-±ï40-0','\',char(fname(i_file))]);
%      saveas(h,[save_directory,'\','6-±ï40-0','\',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h1,[save_directory,'\','6-±ï40-0','\','R1-',char(fname(i_file))]);
%      saveas(h1,[save_directory,'\','6-±ï40-0','\','R1-',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h2,[save_directory,'\','6-±ï40-0','\','R2-',char(fname(i_file))]);
%      saveas(h2,[save_directory,'\','6-±ï40-0','\','R2-',strrep(char(fname(i_file)),'.fig','.png')]);
%      saveas(h3,[save_directory,'\','6-±ï40-0','\','R3-',char(fname(i_file))]);
%      saveas(h3,[save_directory,'\','6-±ï40-0','\','R3-',strrep(char(fname(i_file)),'.fig','.png')]);
%    end
     saveas(h,[save_directory,'\',char(fname(i_file))]);
     saveas(h,[save_directory,'\',strrep(char(fname(i_file)),'.fig','.png')]);
     saveas(h1,[save_directory,'\','R1-',char(fname(i_file))]);
     saveas(h1,[save_directory,'\','\','R1-',strrep(char(fname(i_file)),'.fig','.png')]);
     saveas(h2,[save_directory,'\','R2-',char(fname(i_file))]);
     saveas(h2,[save_directory,'\','R2-',strrep(char(fname(i_file)),'.fig','.png')]);
     saveas(h3,[save_directory,'\','R3-',char(fname(i_file))]);
     saveas(h3,[save_directory,'\','R3-',strrep(char(fname(i_file)),'.fig','.png')]);
    

     close all
end