function [seeds, int, ints] = numseeds(PosList,dotlocations)

%load([PathName '\Pos' num2str(Posnum) '\pos' num2str(Posnum) 'Barcodes11092016.mat']);


seeds = PosList;
int = PosList;
ints = PosList;
% for each cell i
for i = 1:size(dotlocations,2)
    % for each gene j that appears in dotlocations (only genes that are
    % called)
    for j = 1:size(dotlocations(i).cell,1)
        % find the conversion from dotlocation id to global geneID (ind)
        ind = find(strcmpi(dotlocations(i).cell{j,1},PosList));
        % check if there's only one such gene called (i+1 because PosList
        % has gene names in the first column
        if size(PosList{ind,i+1},1) == 1
            % if this is the only time this gene is called, then all
            % seedings are associated with this one call and the number of
            % seeds for this RNA call = the number of seedings called for
            % this gene
            seeds{ind,i+1} = size(dotlocations(i).cell{j,4},1);
%             int{ind,i+1} = dotlocations(i).cell{j,7};
%             ints{ind,i+1} = dotlocations(i).cell{j,6};
        else
            % if a gene is called multiple times (has multiple centroids),
            % then for each FISH dot [dotlocations(i).cell{j,2}] it finds
            % the closest RNA centroid [PosList{ind,i+1}]
            idx = dsearchn(PosList{ind,i+1},dotlocations(i).cell{j,2});
%             holdint = zeros(size(PosList{ind,i+1},1),hybs);
%             holdints = zeros(size(PosList{ind,i+1},1),hybs);
%             for k = 1:size(PosList{ind,i+1},1)
%                 kee = idx ==k;
%                 hybnums = cell2mat(dotlocations(i).cell(j,4));
%                 hybnums = hybnums(kee);
% %                 inttemps = cell2mat(dotlocations(i).cell(j,6));
% %                 inttemp = cell2mat(dotlocations(i).cell(j,7));
% %                 inttemps = inttemps(kee);
% %                 inttemp = inttemp(kee);
% %                 holdint(k,hybnums) = inttemp;
% %                 holdints(k,hybnums) = inttemps;
%             end
            % uses hist to bin the centroid ids associated with each FISH
            % dot, so the result is the number of FISH dots associated with
            % each centroid
            seeds{ind,i+1} = hist(idx,size(PosList{ind,i+1},1));
%             int{ind,i+1} = holdint;
%             ints{ind,i+1} = holdints;
        end
    end
end
%save([PathName '\Pos' num2str(Posnum) '\pos' num2str(Posnum) 'Barcodes11092016.mat'],'seeds','-append');