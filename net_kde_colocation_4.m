function [ Ps_whole ] = net_kde_colocation_4(E,ET,PairPointsData,distance, min_col)
% E={Instance-ID, Feature_type_ID, X,Y};
% ET={FeatureType1_ID...FeatureTypeN_ID};
% PairPointsData = {PairPointID, Instance-ID,Instance-ID2, distance }

ET = sort(ET);
Ps_whole ={};

tic;
CP_List_2_matrix = create_CP_List_2_matrix(ET, E,distance, PairPointsData);
T1  = gen_table_ins(ET,E);  
[candidate_col, CP_List_2_matrix] = kde_2d(CP_List_2_matrix,distance, min_col,T1);%CP_List_2_matrix, distance, min_col,T1,N_max
Ps_whole = get_colocationPs(candidate_col,1,CP_List_2_matrix,Ps_whole,T1, min_col);
toc;
end

%产生1阶表实例
% CurrentT={[1;[1,3,5]];[2;[2,7,9,15]]...}
function T1 = gen_table_ins(ET,E)   
   T1 = cell(length(ET),2);
   temp=size(length(ET),1);
   for FeatureTypeID=1:length(ET)
       T1{FeatureTypeID,1} = FeatureTypeID;% featureID
       T1{FeatureTypeID,2} = find(E(:,2)==ET(FeatureTypeID));%element
       temp(FeatureTypeID)=length(T1{FeatureTypeID,2});
   end 
end

% 由点对数据创建CP_List_2_matrix
function CP_List_2_matrix = create_CP_List_2_matrix(ET, E,r, PairPointsData)
  ET_L = size(ET,1);
   CP_List_2_matrix = cell(ET_L, ET_L);
   for i= 1:size(PairPointsData, 1)
       if PairPointsData(i,4)<r 
           before = E(PairPointsData(i,2),2);
           after = E(PairPointsData(i,3),2);
           instance_ETPair = [before,after]; %[E(instance_ETPair_memberID,2),E(instance_ETPair_memberID2,2)];
           if ~(issorted(instance_ETPair)) && ~isequal(before,after)
               %有方向性的连接不能用union，否则会过滤掉相同的项
                %CP_List_2_matrix{instance_ETPair(2),instance_ETPair(1)} = [CP_List_2_matrix{instance_ETPair(2),instance_ETPair(1)};PairPointsData(i,3),PairPointsData(i,2),PairPointsData(i,4)];
               %无方向性的连接用union
                 if isempty(CP_List_2_matrix{instance_ETPair(1),instance_ETPair(2)})
                     CP_List_2_matrix{instance_ETPair(2),instance_ETPair(1)} = [CP_List_2_matrix{instance_ETPair(2),instance_ETPair(1)};PairPointsData(i,3),PairPointsData(i,2),PairPointsData(i,4)];
                 else
                     CP_List_2_matrix{instance_ETPair(2),instance_ETPair(1)} = union(CP_List_2_matrix{instance_ETPair(2),instance_ETPair(1)},[PairPointsData(i,3),PairPointsData(i,2),PairPointsData(i,4)],'rows');
                 end             
           elseif issorted(instance_ETPair) && ~isequal(before,after)
              % CP_List_2_matrix{instance_ETPair(1),instance_ETPair(2)} = [CP_List_2_matrix{instance_ETPair(1),instance_ETPair(2)};PairPointsData(i,2:4)];
               if isempty(CP_List_2_matrix{instance_ETPair(1),instance_ETPair(2)})
                    CP_List_2_matrix{instance_ETPair(1),instance_ETPair(2)} = [CP_List_2_matrix{instance_ETPair(1),instance_ETPair(2)};PairPointsData(i,2:4)];
               else
                    CP_List_2_matrix{instance_ETPair(1),instance_ETPair(2)} = union(CP_List_2_matrix{instance_ETPair(1),instance_ETPair(2)},PairPointsData(i,2:4),'rows');
               end
           end   
       end
   end
end

% 利用高斯核函数求取二阶colocation patterns/同时去掉CP_List_2_matrix的点对方向，不同方向点对做平均化处理
function [candidate_col, CP_List_2_matrix] = kde_2d(CP_List_2_matrix, distance, min_col,T1)
   % cellfun(@(x) sum(exp(-(x(:,3)/distance).^2/sqrt(2*pi))),CP_List_2_matrix,'UniformOutput', false);
   size_CPmatrix = size(CP_List_2_matrix,1);
   cell_value = zeros(size_CPmatrix,size_CPmatrix);
   for i=1:size_CPmatrix
       for j=i+1:size_CPmatrix
           if ~isempty(CP_List_2_matrix{i,j})
%               CP_List_2_matrix{i,j} = projection([i,j],CP_List_2_matrix([i,j],[i,j]), [i,j]);         
%                closeIndices = find(CP_List_2_matrix{i,j}(:,3)<distance);
%                CP_List_2_matrix{i,j} = CP_List_2_matrix{i,j}(closeIndices,:);              
               CP_List_2_matrix{i,j}(:,3) = exp(-(CP_List_2_matrix{i,j}(:,3)).^2/(2*distance^2));
               proj_ei_w_A = Proj([i,j],i,CP_List_2_matrix([i,j],[i,j]), CP_List_2_matrix{i,j}(:,1:2));%projection([i,j], CP_List_2_matrix([i,j],[i,j]), i);
               proj_ei_w_B = Proj([i,j],j,CP_List_2_matrix([i,j],[i,j]), CP_List_2_matrix{i,j}(:,1:2)); %projection([i,j], CP_List_2_matrix([i,j],[i,j]), j);
               cell_value(i,j) = min(sum(proj_ei_w_A(:,2))/size(T1{i,2},1),sum(proj_ei_w_B(:,2))/size(T1{j,2},1));
           end
       end
   end
   cell_01 = cell_value>min_col;
   candidate_col_m = maximalCliques(cell_01+cell_01');
   candidate_col ={};   
   for i = 1:size(candidate_col_m,2)
    cell_item = find(candidate_col_m(:,i));
    if(length(cell_item)>1)
        candidate_col = [candidate_col,cell_item'];
    end
   end
end

%求取ei类型实例投影在CP上的w（权重），CP为候选模式，ei为投影类型，CP_matrix为对应在CP上的二阶实例表，insT_CP为CP的实例表
function proj_ei_w = Proj(CP,ei,CP_matrix,insT_CP)
Ps_2 = []; % Ps_2存储包含ei的类型对
if size(CP,2)<3
    Ps_2 = [Ps_2; CP];
else
    temp_Ps2 = nchoosek(CP,2);
    %CP剖分为包含ei的类型对，存储在Ps_2中
    for i = 1:size(temp_Ps2,1)
        P2 = temp_Ps2(i,:);
        if ismember(ei,P2)
            Ps_2 = [Ps_2; P2];
        end
    end
end

proj_ei = [];

for i= 1:size(Ps_2,1)
    %找到当前类型对在CP中对应的索引位置，存储在loc中
    [tf, loc] = ismember(Ps_2(i,:), CP);
    %在CP_matrix中提取当前类型对（记作P2）对应的实例对，存储在temp_list中
    temp_list = CP_matrix{loc(1),loc(2)};
    flag_ei = 1; % flag_ei表明ei类型对应temp_list上的index值，另外一个类型记作ej
    if isequal(Ps_2(i,2),ei)
        flag_ei = 2;
    end
    
    %---以ei类型的实例为基础，计算各实例对应的连接边的数量，存储在edgeNum_ei中
    temp_List_insP = temp_list(:,1:2);  %提取二阶实例表中P2对应的instance pair
    unique_insP = unique(temp_List_insP,'rows'); %过滤非重复的instance pair，存储在unique_insP中
    unique_ei = unique(temp_list(:,flag_ei));   %过滤非重复的类型为ei的实例
    unique_ei_cell = num2cell(unique_ei,2);   %将过滤后的ei实例扩展成cell格式，为了运行操作下面的arrayfun和cellfun函数
    %查找unique_ei_cell中每个实例在unique_insP对应的类型为ei列的实例位置，标记为1,否则为0
    independent_logical = arrayfun(@(x) ismember(unique_insP(:,flag_ei),x{:},'rows'), unique_ei_cell, 'UniformOutput', false);  
    %求得各个实例对应的连接边的数量，注意要*2，表明边是有方向的
    independent_edgeNum = cellfun(@(x) sum(x),independent_logical);
    %将各个类型为ei的实例和与其有连接关系的类型为ej的实例的边合并，存储在edgeNum_ei中
    edgeNum_ei = [unique_ei,independent_edgeNum];
    
    %---将temp_list2投影在temp_list的前两元上，记录temp_list2中每个实例对出现在temp_list上的index,并得到kernal值，连接在投影实例对的后面，作为第三列。
    %在insT_CP高阶实例表中提取类型对P2对应的实例对，存储在temp_list2中
    temp_list2 = unique(insT_CP(:,loc),'rows');
    temp_list2_cell = num2cell(temp_list2,2);  %将temp_list2扩展成cell格式，参数2表示以行为单位扩展成cell，比如A=[1,2;2,3],num2cell(A,2)= {[1,2];[2,3]}
    %查找temp_list2_cell中每个实例在temp_List对应的实例对的位置，标记为1,否则为0
    independent_logical = arrayfun(@(x) ismember(temp_List_insP,x{:},'rows'), temp_list2_cell, 'UniformOutput', false);
    %将求得各个非重复投影实例对的kernal值相加，因为存在方向性，所以必须要经过这一步
    independent_kernal = cellfun(@(x) sum(x.*temp_list(:,size(temp_list,2))),independent_logical);
    %将实例对和他们的kernal和进行连接合并，存储在temp_list2_kernal中
    temp_list2_kernal = [temp_list2,independent_kernal];
    
    %---以ei类型的instance为基础，计算各instance对应的kernal值之和，记作ei_K
    %从投影后的实例对中提取类型为ei的实例投影，记作temp_list2_ei
    temp_list2_ei = temp_list2_kernal(:,flag_ei);
    %在temp_list2_ei中寻找二阶实例表中ei类型对应的实例在其上的索引位置，标记为1，否则为0
    independent_logical = arrayfun(@(x) ismember(temp_list2_ei,x{:},'rows'), unique_ei_cell, 'UniformOutput', false);
    %每个ei类型实例对应的kernal做加和处理
    independent_ei_kernal = cellfun(@(x) sum(x.*temp_list2_kernal(:,size(temp_list2_kernal,2))),independent_logical);
    %连接edgeNum_ei和投影kernal之和
    ei_K = [edgeNum_ei,independent_ei_kernal];
    
    %proj_ei存储了三元组序列，分别为：ei类型要素，ei实例投影在二阶实例表的可连接边的个数，ei实例投影在高阶实例表的可连接实例边的kernal值之和
    proj_ei = [proj_ei; ei_K];   
end

%---以proj_ei第一列（类型为ei的实例）为基准，用第三项除以第二项，覆盖原proj_ei(点除)
proj_ei = [proj_ei(:,1),proj_ei(:,3)./proj_ei(:,2)];

%---提取出相同instance对应的权重最小值（w），得到proj_ei_w
unique_instace_ei = unique(proj_ei(:,1),'rows');
unique_instace_ei_cell = num2cell(unique_instace_ei,2); %提取可连接边的非重复项
independent_logical = arrayfun(@(x) find(ismember(proj_ei(:,1),x{:},'rows')), unique_instace_ei_cell, 'UniformOutput', false);
independent_ei_w = cellfun(@(x) min(proj_ei(x,size(proj_ei,2))),independent_logical);

proj_ei_w = [unique_instace_ei,independent_ei_w];

end

%利用高斯核函数求算CP_row是否为colocation pattern
function Ps_bool = kde_Nd(Ins_list,CP_row, CP_List_2_matrix,T1, min_col)
      % Dv_mean_coll = zeros(1,size(Ins_list,1));
       Ps_bool =0;
       Ps_values =[];
%       size_CP_row = size(CP_row, 2);
%       cell_matrix = cell(size_CP_row,size_CP_row);
%        for m= 1:size_CP_row
%            for n=m+1:size_CP_row
%                Pos_inslist = Ins_list(:,[m,n]);
%                %Dv_index = arrayfun(@(x) ismember(x{:},Pos_inslist,'rows'),PairPointsData(:,2),PairPointsData(:,3));
%                [tf, Dv_index] = ismember(Pos_inslist,CP_List_2_matrix{m,n}(:,1:2),'rows');
%                Dv_coll = CP_List_2_matrix{m,n}(Dv_index,3);
%                cell_matrix{m,n} = [Pos_inslist,Dv_coll];
%            end
%        end
       
       for i = 1:size(CP_row, 2)
           cleaned_List = Proj(CP_row,CP_row(i),CP_List_2_matrix,Ins_list);%projection(CP_row, cell_matrix, CP_row(i));
           ps = sum(cleaned_List(:,2))/size(T1{CP_row(i),2},1);
           Ps_values = [Ps_values, ps];
       end
       
       Ps_value = min(Ps_values);
       disp(CP_row);
       disp(Ps_value);
       if Ps_value>min_col
           Ps_bool=1;
       end
end

function Ps_whole = get_colocationPs(CP_List_clear,index,CP_List_2_matrix,Ps_whole,T1,min_col)
  if index<= size(CP_List_clear,2);
    CP_row = CP_List_clear{index};
    if size(CP_row,2)==2  && ~ifCellMember2(CP_row, Ps_whole)
       Ps_whole = [Ps_whole, CP_row];
    end
    if ~ifCellMember2(CP_row, Ps_whole)
        Ins_list = create_Ins_Tree( CP_row, CP_List_2_matrix(CP_row,CP_row)); 
        Ps_bool = kde_Nd(Ins_list,CP_row, CP_List_2_matrix(CP_row,CP_row),T1, min_col);%Ins_list,CP_row, PairPointsData, distance, T1, min_col,N_max
        if Ps_bool
            Ps_whole = [Ps_whole, CP_row];
        else
            CP_row_list = nchoosek(CP_row,size(CP_row,2)-1);
            CP_row_cells = cell(1,size(CP_row_list,1));
            for i = 1:size(CP_row_list,1)
               CP_row_cells{i} = CP_row_list(i,:); 
            end
            Ps_whole = get_colocationPs(CP_row_cells,1, CP_List_2_matrix,Ps_whole,T1, min_col);
        end
    end
    Ps_whole = get_colocationPs(CP_List_clear,index+1, CP_List_2_matrix,Ps_whole, T1,min_col);
  end
end

% 创建Instree树，得到Ins_list
% CP_row= 'ABC';
% CP_List_2_matrix是一个对应于ABC行列的3*3的cells，每个cell里面，如CP_List_2_matrix{1,2}，对应为类型A和类型B的距离小于某个阈值
% 的二阶实例对，并且保存有他们之间的距离
% 即，CP_List_2_matrix{1,2} = [1,4,48934;1,8,4833;7,2,83940;...]
function Ins_list = create_Ins_Tree( CP_row, CP_List_2_matrix)
insTree = tree('root');
i=1;
while i<size(CP_row,2)
    if i==1
       for j = 1:size(CP_List_2_matrix{1,2},1)
           instance_ChildernIds = insTree.getchildren(1);
           items = getTreeContent(insTree, instance_ChildernIds);
           targetnode = find(arrayfun(@(x) isequal(x,CP_List_2_matrix{1,2}(j,1)), items),1);
           if isempty(targetnode)
               sub_tree = tree(CP_List_2_matrix{1,2}(j,1));
               sub_tree = sub_tree.addnode(1,CP_List_2_matrix{1,2}(j,2));
               insTree = insTree.graft(1,sub_tree);
           else
               insTree = insTree.addnode(instance_ChildernIds(targetnode),CP_List_2_matrix{1,2}(j,2));               
           end           
       end
    else
        dt = insTree.depthtree;
        currLevelIDs = find(dt == i); % 获取当前i层级的所有元素
        currLevelItems = getTreeContent(insTree, currLevelIDs);
        CP_List_1 = CP_List_2_matrix{i,i+1}(:,1); 
        loc = cell(size(currLevelItems,2),1);
        parentItems_i =[];
        ConnectItems_i_1 = {};
        for m = 1:size(currLevelItems,2)
            loc{m} = find(currLevelItems(m)==CP_List_1);
            if ~isempty(loc{m})
                parentItems_i = [parentItems_i, currLevelIDs(m)];
                ConnectItems_i_1 = [ConnectItems_i_1, CP_List_2_matrix{i,i+1}(loc{m},2)]; %#ok<*AGROW>
            end
        end
%         [tf,loc] = ismember(currLevelItems,CP_List_1); 
%         parentItems_i = currLevelIDs(find(tf));
%         loc(loc==0)=[]; % 去掉0元素
%         ConnectItems_i_1 = CP_List_2_matrix{i,i+1}(loc,2); % 获取相同
        for j = 1:size(parentItems_i,2)
            for k = 1:size(ConnectItems_i_1{j},1)
                mark = i-1;
                currIns = insTree.getparent(parentItems_i(j));
                while mark>=1
                    targetBoolen = ismember([insTree.get(currIns), ConnectItems_i_1{j}(k)],CP_List_2_matrix{mark,i+1}(:,[1,2]),'rows');
                    if targetBoolen
                        currIns = insTree.getparent(currIns);
                    else
                        break;
                    end
                    mark = mark-1;
                end
                if mark ==0
                   insTree = insTree.addnode(parentItems_i(j),ConnectItems_i_1{j}(k)); 
                end
            end
        end
    end
    i = i+1;        
end
dt = insTree.depthtree;
Ins_list_index = find(dt == size(CP_row,2));
Ins_list = zeros(size(Ins_list_index,2),size(CP_row,2));
for k = 1:size(Ins_list_index,2)
     i=size(CP_row,2);
     currInt = Ins_list_index(k);
     while i>0
         Ins_list(k,i) = insTree.get(currInt);
         currInt = insTree.getparent(currInt);
         i=i-1;
     end
%      path = insTree.findpath(1,Ins_list_index(k));
%      Ins_list(k,:) = getTreeContent(insTree, path(1,2:size(path, 2)));
end
end

% 根据已经进行过清洗的临时实例要素，计算团关系实例
function [tree] = get_pattern_ins(culc_temp,befor_index,after_index, curr_R,tree)
  if(after_index<=size(curr_R,1))
      items = culc_temp{befor_index, after_index};
      if(befor_index ==1)
        tree=[tree,items(:,1:2)];% 若为初始的树，将culc_temp{1,2}赋值给tree
        tree = get_pattern_ins(culc_temp,after_index,after_index+1,curr_R, tree);
      else
          if(~isempty(tree))
            %items_2 = culc_temp{after_index,after_index+1};
            tree_temp=[];
            for i=1:size(tree,1)
               union_item = items(items(:,1)==tree(i,after_index-1),:);%新的clac_item与tree的最末位列的要素进行连接计算，
               if(~isempty(union_item))
                   for k=1:size(union_item,1)
                       flag =1;
    %                    maxlen =max(union_item(k,3),tree(i,after_index));
                       for m=1:befor_index-1
                           [Lia,Locb] = ismember([tree(i,m),union_item(k,2)],culc_temp{m,after_index}(:,1:2),'rows');                     
                           if(~Lia)
                               flag =0;                     
                               break;
    %                        else
    %                            maxlen= culc_temp{m,after_index}(Locb,3);
                           end
                       end
                       if(flag)
                          tree_temp =[tree_temp;[tree(i,1:befor_index),union_item(k,2)]];%maxlen]];% 树的最后一列存储团关系实例的最长边
                       end
                   end 
               end
            end           
          tree = get_pattern_ins(culc_temp,after_index,after_index+1,curr_R, tree_temp);
          end
      end
  end
end

% 判断item是否在targetcells存在包含关系的要素
function result = ifCellMember2(item, targetcells)
   result=0; 
   for i=1:size(targetcells,2)
       if all(ismember(item, targetcells{i}))
           result = result+1;
           break;
       end
   end
end

%从treeCase中批量获取indexs的node内容
function items = getTreeContent(treeCase, indexs)
     items = [];
     for i=indexs
        items = [items, treeCase.get(i)];
     end
end
