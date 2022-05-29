function [xbins,N]=matrix_sort(data)
    
M1=real(data);
n_max=max(M1(:,2));
b=find(M1(:,2)==0);
size_raw=size(b,1);
size_col=n_max+1;
N=zeros(size_raw,size_col);
A=zeros(size_raw,3*(n_max+1));
k=1;
for i=0:n_max
    a=find(M1(:,2)==i);
    for j=1:size(a)
        aa(j,:)=M1(a(j),:);        
    end
    A(1:size(aa,1),k:k+2)=sortrows(aa);
    k=k+3;
    clear aa 
end
N1=A(:,1:3);
N(:,1)=N1(:,3);
for i=2:n_max+1
    c=A(:,(i-1)*3+1:(i-1)*3+3);
    n_star=find(N1(:,1)==c(1,1));
    n_end=find(N1(:,1)==-c(1,1));
    N(n_star:n_end,i)=c(1:(n_end-n_star+1),3);   
    clear c
end
xbins=N1(:,1);
end
% M1=real(data);
% size_raw=2*max(abs(real(M1(:,1))))+1;
% size_col=max(abs(real(M1(:,2))))+1;
% N=zeros(size_raw,size_col);
% a1=find(M1(:,2)==0);
% b1=find(M1(:,2)==1);
% bb1=isempty(b1);
% c1=find(M1(:,2)==2);
% cc1=isempty(c1);
% d1=find(M1(:,2)==3);
% dd1=isempty(d1);
% e1=find(M1(:,2)==4);
% ee1=isempty(e1);
% if bb1==1
%     for i=1:size(a1)
%         aa1(i,:)=M1(a1(i),:);
%     end
%     aa1(i,:)=M1(a1(i),:);
%     A1=sortrows(aa1);
%     N(1:size(aa1),1)=A1(:,3);
% elseif bb1==0 && cc1==1
%     for i=1:size(a1)
%         aa1(i,:)=M1(a1(i),:);
%     end
%     aa1(i,:)=M1(a1(i),:);
%     A1=sortrows(aa1);
%     clear i 
%     for i=1:size(b1)
%         aa2(i,:)=M1(b1(i),:);
%     end
%     A2=sortrows(aa2);
%     n_star2=find(A1(:,1)==A2(1,1));
%     n_end2=find(A1(:,1)==A2(end,1));
%     N(1:size(aa1),1)=A1(:,3);
%     N(n_star2:n_end2,2)=A2(:,3);
% else
%     for i=1:size(a1)
%         aa1(i,:)=M1(a1(i),:);
%     end
%     aa1(i,:)=M1(a1(i),:);
%     A1=sortrows(aa1);
%     clear i 
%     for i=1:size(b1)
%         aa2(i,:)=M1(b1(i),:);
%     end
%     A2=sortrows(aa2);
%     for i=1:size(c1)
%         aa3(i,:)=M1(c1(i),:);
%     end
%     clear i
%     A3=sortrows(aa3);
%     n_star2=find(A1(:,1)==A2(1,1));
%     n_end2=find(A1(:,1)==A2(end,1));
%     n_star3=find(A1(:,1)==A3(1,1));
%     n_end3=find(A1(:,1)==A3(end,1));
%     N(1:size(aa1),1)=A1(:,3);
%     N(n_star2:n_end2,2)=A2(:,3);
%     N(n_star3:n_end3,3)=A3(:,3);
% end
% xbins=A1(:,1);




