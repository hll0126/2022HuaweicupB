clear;close all;clc
original_data = readtable(['E:\Documents\Research_Topic\比赛\2022华为杯\' ...
    'code\dataB1_1.csv']);
original_data2 = original_data(:,1:2);
original_data2(:,3:4) = original_data(:,4:5);
original_data2(:,5) = original_data(:,8); % 每一列分别代表产品id、产品材质、产品长度、产品宽度和批次序号；
data = table2cell(original_data2); % tabel转换为cell类型，方便比较
data = sortrows(data,5); % 按照订单分好的批次对矩阵从批次那列从小到大排序
[data_row,data_col] = size(data);
pici = data(:,5); % 批次
VarTypes = {'int8','table','int8','int8','int8','int8','single','single'};
Varnames = {'批次序号','原片材质','原片序号','产品id','产品x坐标','产品y坐标',...
    '产品x方向长度','产品y方向长度'};

%% 将材质转换为特定数值
for caizhi_num = 1:data_row
    caizhi = data(caizhi_num,2);
    caizhi_char = cell2mat(caizhi);
    caizhi_shuzhi(caizhi_num) = sum(caizhi_char);
end
caizhi_tab = table(data(:,2),caizhi_shuzhi', ...
    'VariableNames',{'caizhi','caizhi_shuzhi'});
data(:,6) = table2cell(caizhi_tab(:,2));
%%
data2 = data;
data2(:,2) = num2cell(caizhi_shuzhi'); % 将材质第二列转换成数值返回给data中
data2(:,end) = [];
data_shuzhi = cell2mat(data2);
table3 = {};
table3 = cell2table(table3);
P_f3 = [];
P_f4 = {};
P_f4 = cell2table(P_f4);
% data_shuzhi2 = data_shuzhi;
for pici_num = 0:max(cell2mat(pici))
    % 后续保存表格需要
    tic
    P_f2 = []; % 板材利用率
    ypcz = {};
    table2 = {};
    table2 = cell2table(table2);
    [pici_row,pici_col] = find(data_shuzhi(:,5)==pici_num); % 找到同一批次的产品
    data_shuzhi2 = data_shuzhi(pici_row(1):pici_row(end),:); % 全部数值的表格
    data_shuzhi22 = data(pici_row(1):pici_row(end),:); % 找到同一批次的所有table
    pici_length(pici_num+1) = size(data_shuzhi22,1); % 获取每一批次各有多少数量的item
    % 找到同一批次里面同一材质内的产品，包括数值和table
    data_shuzhi3 = sortrows(data_shuzhi2,2);
    data_shuzhi33 = sortrows(data_shuzhi22,6);
    data_pici_caizhi = unique(data_shuzhi3(:,2)); % 查看当前批次有多少种材质
    data_pici_caizhi1 = unique(data_shuzhi33(:,2));
    for pici_caizhi_sort = 1:size(data_pici_caizhi1,1)
        % 找到属于同一批次里面同种材质的产品的行和列,
        [pcs_r,pcs_c] = find(data_shuzhi3(:,2)==data_pici_caizhi(1));
        % 在table中寻找符合条件的值
        pat = data_pici_caizhi1(pici_caizhi_sort);
        TF = contains(data_shuzhi33(:,2),pat);
        [TF_r,TF_c] = find(TF);

        data_pici_caizhi_num = data_shuzhi3(TF_r(1):TF_r(end),:);
        [M,N] = size(data_pici_caizhi_num);


        %%%---- 套用第一问的算法进行求解 ----%%%
        W1 = data_pici_caizhi_num(:,4)';
        H1 = data_pici_caizhi_num(:,3)';
        w_max=1220;%宽
        h_max=2440;%长/高
        result = [];
        P_f = [];
        %% 找出符合面积约束的方形组件
        changdu = 0;
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
            ddcs_max = 1;          % 最大迭代次数
            Lujin_best = zeros(ddcs_max,n);      % 各代最佳路径
            L_best = zeros(ddcs_max,1);     % 各代最佳路径的长度
            L_ave = zeros(ddcs_max,1);      % 各代路径的平均长度
            for ddcs=1:ddcs_max%在ddcs小于ddcs_max前，一直循环
                %% 随机产生各个蚂蚁的起点
                %     tic
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
                %     t(ddcs) = toc;
            end
%             figure
%             plot(L_best)
%             hold on
%             plot(L_ave)
%             legend('迭代最优值','平均值')
%             title('蚁群算法寻优过程')
            [F,A]=my_mathmode(Lujin_best(end,:),W,H,R,w_max,h_max,changdu2);
            result(changdu-changdu2+1:changdu,1) = ii-1; % 原片序号：从0开始
            % 产品id：产品id相当于订单号，需要一一对应
            result(changdu-changdu2+1:changdu,2) = data_pici_caizhi_num(changdu-changdu2+1:changdu,1);
            result(changdu-changdu2+1:changdu,3) = A(:,1);
            result(changdu-changdu2+1:changdu,4) = A(:,3);
            result(changdu-changdu2+1:changdu,5) = A(:,2);
            result(changdu-changdu2+1:changdu,6) = A(:,4);
            P_f(ii) = F; % 产品利用率
            %%
            %下面进行制图
            figure_n = (pici_num+1)*100 + pici_caizhi_sort*10 + ii;
            figure(figure_n)
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

            %求板材利用率.
            aaa='%';
            fprintf('板材的利用率为：\n %f %s \n',F,aaa);
            title(['排样图  空间利用率为',num2str(F),'%'],'FontSize',18)
            saveas(figure(figure_n),strcat('dataB1/',num2str(pici_num), ...
                num2str(pici_caizhi_sort),num2str(ii),'.png')); % 保存每块板样的空间利用率图片


            %% 终止条件
            if changdu == M
                break;
            else
                ii = ii +1;
            end
        end
        % sz = [M,7];
        VarTypes = {'int8','table','int8','int8','int8','int8','single','single'};
        Varnames = {'批次序号','原片材质','原片序号','产品id','产品x坐标','产品y坐标',...
            '产品x方向长度','产品y方向长度'};
        ypcz = pat;
        ypcz = repmat(ypcz,[1,M]);
        ypcz = ypcz'; %原片材质
        pcxh = repmat(pici_num,[M,1]);
        table1 = table(pcxh,ypcz,result(:,1),result(:,2),result(:,3),result(:,4), ...
            result(:,5),result(:,6),'VariableNames',Varnames);
        table2 = [table2;table1]; % 同一批次内同种材质的所有产品
        P_f2 = [P_f2,P_f]; % 同一批次所有板材利用率
    end
    table3 = [table3;table2]; % 所有批次
    % 所有板材利用率
    P_f3 = table(P_f2');
    P_f4 = [P_f4;P_f3];
    P_f44(pici_num+1) = mean(P_f2);
    time_pici(pici_num+1) = toc;
end

P_f5 = table2array(P_f4); % 全部板材利用率
P_f5_max = max(P_f5); % 最大值
P_f5_avg = mean(P_f5); % 均值
%% 保存为表格
table33 = table2cell(table3);
table4 = table(table33(:,1),table33(:,2),table33(:,3),table33(:,4),table33(:,5), ...
    table33(:,6),table33(:,7),table33(:,8),'VariableNames',Varnames);
% writetable(table3,'resultB1.csv')


























% %% 按照板样材料类型来区分
% material_n = data(:,2);
% material_n2 = [];
% for i = 1:M
%     material_n1 = cell2mat(material_n(i));
%     material_n11 = sum(material_n1);
%     material_n2(i) = material_n11;
% end
% material_n3 = length(unique(material_n2)); % 材质类型数量
% data2 = data(:,1);
% data2(:,2:3) = data(:,4:5);
% data2 = cell2mat(data2);
% data2(:,4) = material_n2'; % 最后一列为材料类型，转换成了数字
% data22 = sortrows(data2,4); % 按照材料来分组的
%
% %% 求取不同订单数
% order_n = data(:,6);
% order_n2 = [];
% for i = 1:M
%     order_n11 = cell2mat(order_n(i));
%     order_n111 = str2num(order_n11(6:end));
%     order_n2(i) = order_n111;
% end
% order_n3 = length(unique(order_n2));
% data2(:,5) = order_n2';
% data3 = sortrows(data2,5); % 按照订单来分组的
%
%




