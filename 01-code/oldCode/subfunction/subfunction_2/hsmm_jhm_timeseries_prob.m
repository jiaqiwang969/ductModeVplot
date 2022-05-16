function [state,state2,stateseq,lkh]=hsmm_jhm_timeseries_prob(feature1,PAI,A,P,mu,Sigma,mixmat)
%IterationNo迭代次数，0时仅计算最优状态序列和似然概率，大于0时重估参数且计算最优状态序列和似然概率
[O,T,nex]=size(feature1);%时间序列特征的维数
D=size(P,2);
cov_type='full';
%初始状态概率向量、状态转移矩阵，平均初始化
[O,Q,M]=size(mu);
%中间计算量初始化
ALPHA=zeros(Q,D);%前向变量，alpha_t(i,d)=P(q_t=i,tao_t=d|o_1_t-1)
b=zeros(Q,T); %b_*_i(o_t)=b_i(o_t)/P(o_t|o_1_t-1)
S=zeros(Q,T);%S_t(i)=P(tao_t=1,q_t+1=i|o_1_t)
E=zeros(Q,T);%E_t(i)=P(q_t=i,tao_t=1|o_1_t)
BETA=ones(Q,D);
Ex=ones(Q,D);
Sx=ones(Q,D);
GAMMA=zeros(Q,1); 
Qest=zeros(T,1);
data=num2cell(feature1,[1 2]);
lkh=0;  
stateseq=[];
for on=1:nex	 % for each observation sequence
    obs=data{on};	 % the n'th observation sequence
    [B,B2]=mixgauss_prob(obs,mu,Sigma,mixmat);
    % B(i,t) = Pr(y(t) | Q(t)=i) 
    % B2(i,k,t) = Pr(y(t) | Q(t)=i, M(t)=k) 
    %    starttime=clock;
    %++++++++++++++++++     Forward     +++++++++++++++++
    %---------------    Initialization    ---------------
    ALPHA(:)=0; ALPHA=repmat(PAI,1,D).*P;		%Equation (13)t=1时alpha
    r=(B(:,1)'*sum(ALPHA,2));			%Equation (3)t=1时r_t_-1
    b(:,1)=B(:,1)./r;				%Equation (2)t=1时b^*_i(o_t)
     
    E(:)=0; E(:,1)=b(:,1).*ALPHA(:,1);		%Equation (5)
    S(:)=0; S(:,1)=A'*E(:,1);			%Equation (6)
    lkh=lkh+log(r);
    %---------------    Induction    ---------------
    for t=2:T
        ALPHA=[repmat(S(:,t-1),1,D-1).*P(:,1:D-1)+repmat(b(:,t-1),1,D-1).*ALPHA(:,2:D),S(:,t-1).*P(:,D)];		%Equation (12)
        r=(B(:,t)'*sum(ALPHA,2));		%Equation (3)
%         if r==0
%            pause 
%         end
        b(:,t)=B(:,t)./r;			%Equation (2)
        E(:,t)=b(:,t).*ALPHA(:,1);		%Equation (5)
        S(:,t)=A'*E(:,t);				%Equation (6)
        lkh=lkh+log(r);

    end
        
    %++++++++ To check if the likelihood is increased ++++++++
    %         if ir>1
    %             %    clock-starttime
    %             if (lkh-lkh1)/T<0.001
    %                 break
    %             end
    %         end
    %         lkh1=lkh;
    %++++++++ Backward and Parameter Restimation ++++++++
    %---------------    Initialization    ---------------
    GAMMA=b(:,T).*sum(ALPHA,2);
    [X,Qest(T)]=max(GAMMA);
    BETA=repmat(b(:,T),1,D);				%Equation (7)
    Ex=sum(P.*BETA,2);					%Equation (8)
    Sx=A*Ex;						%Equation (9)
    %---------------    Induction    ---------------
    for t=(T-1):-1:1
        %% for estimate of state at time t
        GAMMA=GAMMA+E(:,t).*Sx-S(:,t).*Ex;
        GAMMA(GAMMA<0)=0;           % eleminate errors due to inaccurace of the computation.
        [X,Qest(t)]=max(GAMMA);
        BETA=repmat(b(:,t),1,D).*[Sx,BETA(:,1:D-1)];	%Equation (14)
        Ex=sum(P.*BETA,2);					%Equation (8)
        Sx=A*Ex;						%Equation (9)
    end
    stateseq=[stateseq Qest];
end								% End for multiple observation sequences
for i=1:Q
    num(i)=length(find(stateseq==i));
end
state=find(num==max(num));
state2=Qest(T);