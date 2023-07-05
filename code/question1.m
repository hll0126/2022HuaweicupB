clear;close all;clc
% for cishu = 1:10
    original_data = readmatrix(['E:\Documents\Research_Topic\比赛\2022华为杯\' ...
        '2022年中国研究生数学建模竞赛试题\2022年B题\子问题1-数据集A\dataA4.csv']); % 读取数据
    [M,N] = size(original_data); % 求取数据规格
    original_data2 = original_data(~isnan(original_data)); % 删除material和order两列
    data = reshape(original_data2,[M,4]);
    W1 = data(:,4)';
    H1 = data(:,3)';
    w_max=1220; % 宽
    h_max=2440; % 长/高
    result = [];
    P_f = [];
    %% 找出符合面积约束的方形组件
    changdu = 0;
    tic
    for ii = 1:M
        area1 = cumsum(W1(:,changdu+1:end)'.*H1(:,changdu+1:end)');
        weizhi = find(area1<w_max*h_max);
        changdu2 = length(weizhi);
        changdu = changdu + changdu2;
        area1 = [];
        W = W1(:,changdu-changdu2+1:changdu);
        H = H1(:,changdu-changdu2+1:changdu);

        R1 = ones(changdu2,1);
        R2 = round(rand(1,changdu2));
        R3 = round(rand(1,changdu2));
        R=R3;
        P=[];
        yfw=[];

        %% 初始化参数
        m = 50;         % 蚂蚁数量，视情况而定，坐标点多的话可以适当增加蚂蚁数量
        a= 1;           % 信息素重要程度因子
        b= 5;           % 启发函数重要程度因子
        r = 0.1;        % 信息素挥发因子
        Q = 100;          % 常数
        n=length(W);
        qfhs = ones(n,n)./n;    % 启发函数
        xxsjz = ones(n,n);       % 信息素矩阵初始化
        ljjl = zeros(m,n);       % 路径记录表矩阵初始化
        ddcs_max = 100;          % 最大迭代次数
        Lujin_best = zeros(ddcs_max,n);      % 各代最佳路径
        L_best = zeros(ddcs_max,1);     % 各代最佳路径的长度
        L_ave = zeros(ddcs_max,1);      % 各代路径的平均长度
        for ddcs=1:ddcs_max%在ddcs小于ddcs_max前，一直循环
            %% 随机产生各个蚂蚁的起点
            start = zeros(m,1);
            for i = 1:m
                temp = randperm(changdu2);%功能是随机打乱一个数字序列,也就是现将坐标点排号再打乱，相当于将蚂蚁随机分布在各个地点
                start(i) = temp(1);
            end
            ljjl(:,1) = start;
            %% 构建解空间
            zuobiao_index = 1:n;
            % 逐个蚂蚁路径选择
            for i = 1:m
                % 逐个地点路径选择
                for j = 2:n
                    yfw = ljjl(i,1:(j - 1));           % 已访问的地点集合(禁忌表)
                    allow_index = ~ismember(zuobiao_index,yfw);%ismember用于判断矩阵某个元素是否存在，用法详见后文函数讲解
                    allow = zuobiao_index(allow_index);  % 待访问的城市集合
                    P = allow;
                    % 计算节点转移概率
                    for k = 1:length(allow)
                        P(k) = xxsjz(yfw(end),allow(k))^a * qfhs(yfw(end),allow(k))^b;
                    end
                    P = P/sum(P);
                    % 选择下一节点
                    Plj = cumsum(P);
                    yidong_index = find(Plj >= rand);
                    yidong = allow(yidong_index(1));
                    ljjl(i,j) = yidong;
                end
            end
            % 计算各个蚂蚁的路径距离
            for i=1:m
                L(i,1)=my_mathmode(ljjl(i,:),W,H,R,w_max,h_max,changdu2);
            end
            % 计算最短路径距离及平均距离
            if ddcs == 1
                [max_L,max_index] = max(L);
                L_best(ddcs) = max_L;
                L_ave(ddcs) = mean(L);
                Lujin_best(ddcs,:) = ljjl(max_index,:);
            else
                [max_L,max_index] = max(L);
                L_best(ddcs) = max(L_best(ddcs - 1),max_L);
                L_ave(ddcs) = mean(L);
                if L_best(ddcs) == max_L
                    Lujin_best(ddcs,:) = ljjl(max_index,:);
                else
                    Lujin_best(ddcs,:) = Lujin_best((ddcs-1),:);
                end
            end
            %% 更新信息素
            S = zeros(n,n);
            % 逐个蚂蚁计算
            for i = 1:m
                % 逐个节点计算
                for j = 1:(n - 1)
                    S(ljjl(i,j),ljjl(i,j+1)) = S(ljjl(i,j),ljjl(i,j+1)) + Q/L(i);
                end
                S(ljjl(i,n),ljjl(i,1)) = S(ljjl(i,n),ljjl(i,1)) + Q/L(i);
            end
            xxsjz = (1-r) * xxsjz + S;
            % 迭代次数加1，清空路径记录表
            ljjl = zeros(m,n);
            t(ddcs) = toc;
        end
%         figure
%         plot(L_best)
%         hold on
%         plot(L_ave)
%         legend('迭代最优值','平均值')
%         title('蚁群算法寻优过程')
        [F,A]=my_mathmode(Lujin_best(end,:),W,H,R,w_max,h_max,changdu2);
        result(changdu-changdu2+1:changdu,1) = ii-1; % 原片序号：从0开始
        result(changdu-changdu2+1:changdu,2) = data(changdu-changdu2+1:changdu,1); % 产品id：从0开始
        result(changdu-changdu2+1:changdu,3) = A(:,1);
        result(changdu-changdu2+1:changdu,4) = A(:,3);
        result(changdu-changdu2+1:changdu,5) = A(:,2);
        result(changdu-changdu2+1:changdu,6) = A(:,4);
        P_f(ii) = F; % 产品利用率
        %%
        %下面进行制图
        figure
        rectangle('position',[0,0,w_max,h_max]);axis equal;
        hold on;
        x=zeros(1,changdu2);y=zeros(1,changdu2);
        for i=1:changdu2
            if A(i,3)>0&&A(i,4)>0
                rectangle('position',A(i,:),'facecolor','y');
            end
            x(1,i)=A(i,1)+A(i,3)/2;
            y(1,i)=A(i,2)+A(i,4)/2;
        end
        text(x,y,num2str(Lujin_best(end,:)'),'FontSize',8);
        hold off;
%         求板材利用率.
        aaa='%';
        fprintf('板材的利用率为：\n %f %s \n',F,aaa);
        title(['排样图  空间利用率为',num2str(F),'%',],'FontSize',18)
        saveas(figure(ii),strcat('dataA4/',num2str(ii),'.png')); % 保存每块板样的空间利用率图片
        %% 终止条件
        if changdu == M
            break;
        else
            ii = ii +1;
        end

    end

%     t1(cishu) = toc

% end
sz = [M,7];
VarTypes = {'table','int8','int8','int8','int8','single','single'};
Varnames = {'原片材质','原片序号','产品id','产品x坐标','产品y坐标',...
    '产品x方向长度','产品y方向长度'};
ypcz = {'YW10-0218S'};
ypcz = repmat(ypcz,[M,1]);
table1 = table(ypcz,result(:,1),result(:,2),result(:,3),result(:,4), ...
    result(:,5),result(:,6),'VariableNames',Varnames);
% writetable(table1,'resultA1.csv') % 保存数据

