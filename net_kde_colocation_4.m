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

%����1�ױ�ʵ��
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

% �ɵ�����ݴ���CP_List_2_matrix
function CP_List_2_matrix = create_CP_List_2_matrix(ET, E,r, PairPointsData)
  ET_L = size(ET,1);
   CP_List_2_matrix = cell(ET_L, ET_L);
   for i= 1:size(PairPointsData, 1)
       if PairPointsData(i,4)<r 
           before = E(PairPointsData(i,2),2);
           after = E(PairPointsData(i,3),2);
           instance_ETPair = [before,after]; %[E(instance_ETPair_memberID,2),E(instance_ETPair_memberID2,2)];
           if ~(issorted(instance_ETPair)) && ~isequal(before,after)
               %�з����Ե����Ӳ�����union���������˵���ͬ����
                %CP_List_2_matrix{instance_ETPair(2),instance_ETPair(1)} = [CP_List_2_matrix{instance_ETPair(2),instance_ETPair(1)};PairPointsData(i,3),PairPointsData(i,2),PairPointsData(i,4)];
               %�޷����Ե�������union
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

% ���ø�˹�˺�����ȡ����colocation patterns/ͬʱȥ��CP_List_2_matrix�ĵ�Է��򣬲�ͬ��������ƽ��������
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

%��ȡei����ʵ��ͶӰ��CP�ϵ�w��Ȩ�أ���CPΪ��ѡģʽ��eiΪͶӰ���ͣ�CP_matrixΪ��Ӧ��CP�ϵĶ���ʵ����insT_CPΪCP��ʵ����
function proj_ei_w = Proj(CP,ei,CP_matrix,insT_CP)
Ps_2 = []; % Ps_2�洢����ei�����Ͷ�
if size(CP,2)<3
    Ps_2 = [Ps_2; CP];
else
    temp_Ps2 = nchoosek(CP,2);
    %CP�ʷ�Ϊ����ei�����Ͷԣ��洢��Ps_2��
    for i = 1:size(temp_Ps2,1)
        P2 = temp_Ps2(i,:);
        if ismember(ei,P2)
            Ps_2 = [Ps_2; P2];
        end
    end
end

proj_ei = [];

for i= 1:size(Ps_2,1)
    %�ҵ���ǰ���Ͷ���CP�ж�Ӧ������λ�ã��洢��loc��
    [tf, loc] = ismember(Ps_2(i,:), CP);
    %��CP_matrix����ȡ��ǰ���Ͷԣ�����P2����Ӧ��ʵ���ԣ��洢��temp_list��
    temp_list = CP_matrix{loc(1),loc(2)};
    flag_ei = 1; % flag_ei����ei���Ͷ�Ӧtemp_list�ϵ�indexֵ������һ�����ͼ���ej
    if isequal(Ps_2(i,2),ei)
        flag_ei = 2;
    end
    
    %---��ei���͵�ʵ��Ϊ�����������ʵ����Ӧ�����ӱߵ��������洢��edgeNum_ei��
    temp_List_insP = temp_list(:,1:2);  %��ȡ����ʵ������P2��Ӧ��instance pair
    unique_insP = unique(temp_List_insP,'rows'); %���˷��ظ���instance pair���洢��unique_insP��
    unique_ei = unique(temp_list(:,flag_ei));   %���˷��ظ�������Ϊei��ʵ��
    unique_ei_cell = num2cell(unique_ei,2);   %�����˺��eiʵ����չ��cell��ʽ��Ϊ�����в��������arrayfun��cellfun����
    %����unique_ei_cell��ÿ��ʵ����unique_insP��Ӧ������Ϊei�е�ʵ��λ�ã����Ϊ1,����Ϊ0
    independent_logical = arrayfun(@(x) ismember(unique_insP(:,flag_ei),x{:},'rows'), unique_ei_cell, 'UniformOutput', false);  
    %��ø���ʵ����Ӧ�����ӱߵ�������ע��Ҫ*2�����������з����
    independent_edgeNum = cellfun(@(x) sum(x),independent_logical);
    %����������Ϊei��ʵ�������������ӹ�ϵ������Ϊej��ʵ���ıߺϲ����洢��edgeNum_ei��
    edgeNum_ei = [unique_ei,independent_edgeNum];
    
    %---��temp_list2ͶӰ��temp_list��ǰ��Ԫ�ϣ���¼temp_list2��ÿ��ʵ���Գ�����temp_list�ϵ�index,���õ�kernalֵ��������ͶӰʵ���Եĺ��棬��Ϊ�����С�
    %��insT_CP�߽�ʵ��������ȡ���Ͷ�P2��Ӧ��ʵ���ԣ��洢��temp_list2��
    temp_list2 = unique(insT_CP(:,loc),'rows');
    temp_list2_cell = num2cell(temp_list2,2);  %��temp_list2��չ��cell��ʽ������2��ʾ����Ϊ��λ��չ��cell������A=[1,2;2,3],num2cell(A,2)= {[1,2];[2,3]}
    %����temp_list2_cell��ÿ��ʵ����temp_List��Ӧ��ʵ���Ե�λ�ã����Ϊ1,����Ϊ0
    independent_logical = arrayfun(@(x) ismember(temp_List_insP,x{:},'rows'), temp_list2_cell, 'UniformOutput', false);
    %����ø������ظ�ͶӰʵ���Ե�kernalֵ��ӣ���Ϊ���ڷ����ԣ����Ա���Ҫ������һ��
    independent_kernal = cellfun(@(x) sum(x.*temp_list(:,size(temp_list,2))),independent_logical);
    %��ʵ���Ժ����ǵ�kernal�ͽ������Ӻϲ����洢��temp_list2_kernal��
    temp_list2_kernal = [temp_list2,independent_kernal];
    
    %---��ei���͵�instanceΪ�����������instance��Ӧ��kernalֵ֮�ͣ�����ei_K
    %��ͶӰ���ʵ��������ȡ����Ϊei��ʵ��ͶӰ������temp_list2_ei
    temp_list2_ei = temp_list2_kernal(:,flag_ei);
    %��temp_list2_ei��Ѱ�Ҷ���ʵ������ei���Ͷ�Ӧ��ʵ�������ϵ�����λ�ã����Ϊ1������Ϊ0
    independent_logical = arrayfun(@(x) ismember(temp_list2_ei,x{:},'rows'), unique_ei_cell, 'UniformOutput', false);
    %ÿ��ei����ʵ����Ӧ��kernal���Ӻʹ���
    independent_ei_kernal = cellfun(@(x) sum(x.*temp_list2_kernal(:,size(temp_list2_kernal,2))),independent_logical);
    %����edgeNum_ei��ͶӰkernal֮��
    ei_K = [edgeNum_ei,independent_ei_kernal];
    
    %proj_ei�洢����Ԫ�����У��ֱ�Ϊ��ei����Ҫ�أ�eiʵ��ͶӰ�ڶ���ʵ����Ŀ����ӱߵĸ�����eiʵ��ͶӰ�ڸ߽�ʵ����Ŀ�����ʵ���ߵ�kernalֵ֮��
    proj_ei = [proj_ei; ei_K];   
end

%---��proj_ei��һ�У�����Ϊei��ʵ����Ϊ��׼���õ�������Եڶ������ԭproj_ei(���)
proj_ei = [proj_ei(:,1),proj_ei(:,3)./proj_ei(:,2)];

%---��ȡ����ͬinstance��Ӧ��Ȩ����Сֵ��w�����õ�proj_ei_w
unique_instace_ei = unique(proj_ei(:,1),'rows');
unique_instace_ei_cell = num2cell(unique_instace_ei,2); %��ȡ�����ӱߵķ��ظ���
independent_logical = arrayfun(@(x) find(ismember(proj_ei(:,1),x{:},'rows')), unique_instace_ei_cell, 'UniformOutput', false);
independent_ei_w = cellfun(@(x) min(proj_ei(x,size(proj_ei,2))),independent_logical);

proj_ei_w = [unique_instace_ei,independent_ei_w];

end

%���ø�˹�˺�������CP_row�Ƿ�Ϊcolocation pattern
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

% ����Instree�����õ�Ins_list
% CP_row= 'ABC';
% CP_List_2_matrix��һ����Ӧ��ABC���е�3*3��cells��ÿ��cell���棬��CP_List_2_matrix{1,2}����ӦΪ����A������B�ľ���С��ĳ����ֵ
% �Ķ���ʵ���ԣ����ұ���������֮��ľ���
% ����CP_List_2_matrix{1,2} = [1,4,48934;1,8,4833;7,2,83940;...]
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
        currLevelIDs = find(dt == i); % ��ȡ��ǰi�㼶������Ԫ��
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
%         loc(loc==0)=[]; % ȥ��0Ԫ��
%         ConnectItems_i_1 = CP_List_2_matrix{i,i+1}(loc,2); % ��ȡ��ͬ
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

% �����Ѿ����й���ϴ����ʱʵ��Ҫ�أ������Ź�ϵʵ��
function [tree] = get_pattern_ins(culc_temp,befor_index,after_index, curr_R,tree)
  if(after_index<=size(curr_R,1))
      items = culc_temp{befor_index, after_index};
      if(befor_index ==1)
        tree=[tree,items(:,1:2)];% ��Ϊ��ʼ��������culc_temp{1,2}��ֵ��tree
        tree = get_pattern_ins(culc_temp,after_index,after_index+1,curr_R, tree);
      else
          if(~isempty(tree))
            %items_2 = culc_temp{after_index,after_index+1};
            tree_temp=[];
            for i=1:size(tree,1)
               union_item = items(items(:,1)==tree(i,after_index-1),:);%�µ�clac_item��tree����ĩλ�е�Ҫ�ؽ������Ӽ��㣬
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
                          tree_temp =[tree_temp;[tree(i,1:befor_index),union_item(k,2)]];%maxlen]];% �������һ�д洢�Ź�ϵʵ�������
                       end
                   end 
               end
            end           
          tree = get_pattern_ins(culc_temp,after_index,after_index+1,curr_R, tree_temp);
          end
      end
  end
end

% �ж�item�Ƿ���targetcells���ڰ�����ϵ��Ҫ��
function result = ifCellMember2(item, targetcells)
   result=0; 
   for i=1:size(targetcells,2)
       if all(ismember(item, targetcells{i}))
           result = result+1;
           break;
       end
   end
end

%��treeCase��������ȡindexs��node����
function items = getTreeContent(treeCase, indexs)
     items = [];
     for i=indexs
        items = [items, treeCase.get(i)];
     end
end
