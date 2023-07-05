function [F,A]=my_mathmode(P,W,H,R,w_max,h_max,changdu)
juxing_size=zeros(changdu,2);%表示小矩形的尺寸.M(i,1)、M(i,2)分别表示序号为 i 的小矩形的宽和高.
D=zeros(changdu,2);%排样方案数组.D(i,1)、D(i,2)分别表示第i个矩形的序号以及r(i).
rest=zeros(changdu,4);%剩余矩形数组.记录每个矩形的左下角坐标（x,y）、宽和高.
A=zeros(changdu,4);%记录数组（记录每个矩形件在样板上的位置）.
%初始值均为零。（zeros）
juxing_size=[W',H'];%小矩形件的尺寸.
D=[P',R'];     %选择排样方案.

w=zeros(changdu,1);h=zeros(changdu,1);%矩形件的宽和高
N=1;  %N是剩余矩形的个数.
rest(1,:)=[0,0,w_max,h_max];%样板的尺寸设为：宽15，高60.
HH=0;% 占用高度，用来求样板利用率.
for i=1:changdu %矩形件i的宽和高
    w(i)=juxing_size(D(i,1),1+D(i,2));h(i)=juxing_size(D(i,1),2-D(i,2));
end
re=zeros(30,4);
for i=1:changdu  %放置changdu个矩形件
    %i=3;
    n=1;j=0;
    while n<=N  %求能包含小矩形件的剩余矩形 re，j表示个数.
        if rest(n,3)>=w(i)&&rest(n,4)>=h(i)
            j=j+1;re(j,:)=rest(n,:);
        end
        n=n+1;
    end
    k=2;
    while k<=j %求用到的剩余矩形（放在re的第一行），根据BL条件.
        if re(k,2)<re(1,2)||(re(k,2)==re(1,2)&&re(k,1)<re(1,1))
            re(1,:)=re(k,:);
        end
        k=k+1;
    end
    A(i,1)=re(1,1); A(i,2)=re(1,2);A(i,3)=w(i);A(i,4)=h(i);%第i个矩形件的位置信息
    if HH<(A(i,2)+A(i,4))  % 占用高度的更新.
        HH=(A(i,2)+A(i,4));
    end
    %剩余矩形数组的处理更新.N A(i,)
    n=1;Now=N;
    while n<=Now
        m=0;
        k=N+1;         %
        if (A(i,2)>rest(n,2))&&(A(i,2)<(rest(n,2)+rest(n,4)))&&(A(i,1)<rest(n,1)+rest(n,3))&&(A(i,1)+A(i,3)>rest(n,1))
            N=N+1;
            rest(k,1)=rest(n,1);rest(k,2)=rest(n,2);
            rest(k,3)=rest(n,3);rest(k,4)=A(i,2)-rest(n,2);
            k=k+1;m=1;
        end
        if rest(n,1)+rest(n,3)>A(i,1)+A(i,3)&&A(i,1)+A(i,3)>rest(n,1)&&(A(i,2)<rest(n,2)+rest(n,4))&&(A(i,2)+A(i,4)>rest(n,2))
            N=N+1;
            rest(k,1)=A(i,1)+A(i,3);rest(k,2)=rest(n,2);
            rest(k,3)=rest(n,1)+rest(n,3)-A(i,1)-A(i,3);
            rest(k,4)=rest(n,4);
            k=k+1;m=1;
        end
        if A(i,1)>rest(n,1)&&(A(i,1)<rest(n,1)+rest(n,3))&&(A(i,2)<rest(n,2)+rest(n,4))&&(A(i,2)+A(i,4)>rest(n,2))
            N=N+1;
            rest(k,1)=rest(n,1);rest(k,2)=rest(n,2);
            rest(k,4)=rest(n,4);rest(k,3)=A(i,1)-rest(n,1);
            k=k+1;m=1;
        end
        if rest(n,2)+rest(n,4)>A(i,2)+A(i,4)&&A(i,2)+A(i,4)>rest(n,2)&&(A(i,1)<rest(n,1)+rest(n,3))&&(A(i,1)+A(i,3)>rest(n,1))
            N=N+1;
            rest(k,2)=A(i,2)+A(i,4);rest(k,1)=rest(n,1);
            rest(k,4)=rest(n,2)+rest(n,4)-A(i,2)-A(i,4);
            rest(k,3)=rest(n,3);
            k=k+1;m=1;
        end
        if m==0
            n=n+1;
        else
            rest(n,:)=[];N=N-1;Now=Now-1;rest(50,:)=0;%删除要注意！
        end
    end  %求出剩余矩形数组，下面要对其进行整理.

    n=1;m=0;
    while n<=N
        k=i+1;
        while k<=changdu
            if (rest(n,3)>=w(k)&&rest(n,4)>=h(k))
                m=1;break;
            else
                k=k+1;%

            end

        end
        if m==1
            n=n+1;
        else
            N=N-1;rest(n,:)=[];rest(50,:)=0;
            %删去面积为零的或已无法排下所剩的任何一个矩形件的剩余矩形
        end
    end
    n=1;
    while n<=N
        k=1;
        while k<N
            if rest(n,1)>=rest(k,1)&&rest(n,2)>=rest(k,2)&&rest(n,1)+rest(n,3)<=rest(k,1)+rest(k,3)&&rest(n,2)+rest(n,4)<=rest(k,2)+rest(k,4)
                if n~=k
                    N=N-1;m=1;
                    rest(n,:)=[];rest(50,:)=0;%删去被包含的剩余矩形
                    break;
                else
                    k=k+1;m=0;
                end
            else
                k=k+1;m=0;
            end
        end
        k=k+1;m=0;
        if m==0
            n=n+1;
        end
    end  %
end
F=sum(W.*H)/(w_max*h_max)*100;
end

