function [numtotal final_PosList] = filterseeds_marked_up(numseeds,seedlist,orig_PosList)
b = cellfun(@(x) x >= numseeds,seedlist(:,2:end),'UniformOutput',0);
num = cellfun(@sum,b);
final_PosList = orig_PosList;
for i = 1:length(b)
    if ~isempty(b{i})
        dIdx = find(~b{i});
        final_PosList{i,2}(dIdx,:) = [];
    end
end
final_PosList = final_PosList(:,2);
numtotal = sum(num,2);
