function [dotlocations,copynumfinalrevised,PosList] = pointLocations(PathName,posnum,hyb, channum, points, foundbarcodes,barcodekey,copynumfinal,radius)

copynumfinalrevised = copynumfinal;
fullpath = [PathName filesep 'RoiSet_Pos' num2str(posnum) ];
vertex = selfseg(fullpath);
PosList = cell(size(barcodekey.names,1),length(vertex)+1);
PosList(:,1) = copynumfinal(:,1);

for i = 1:length(vertex)
    g = {};
    allcodes = [];
    allpoints = [];
    ptchannel = [];
    hybnum = [];
    ptint = [];
    ptintabs = [];
    for j = 1:hyb
        for k = 1:channum
            if ~isempty(foundbarcodes(j).found(k).consensus)
                include = inpolygon(points(j).dots(k).channels(:,1),points(j).dots(k).channels(:,2),vertex(i).x,vertex(i).y);
                foundint = foundbarcodes(j).found(k).consensus > 0;
                keepdgf = include & foundint;
                allcodes = [allcodes; foundbarcodes(j).found(k).consensus(keepdgf)];
                allpoints = [allpoints; points(j).dots(k).channels(keepdgf,:)];
                ptchannel = [ptchannel; k*ones(size(foundbarcodes(j).found(k).consensus(keepdgf),1),1)];
                ptint = [ptint;ones(size(foundbarcodes(j).found(k).consensus(keepdgf)))];
                ptintabs = [ptintabs;ones(size(foundbarcodes(j).found(k).consensus(keepdgf)))];
                hybnum = [hybnum; j*ones(size(foundbarcodes(j).found(k).consensus(keepdgf),1),1)];
            end
        end
    end
    fun = unique(allcodes);
    for l = 1:length(fun); 
        g{l,2} = allpoints(allcodes==fun(l),:); 
        g{l,1} = barcodekey.names{fun(l)};
        g{l,3} = ptchannel(allcodes==fun(l),:);
        g{l,4} = hybnum(allcodes==fun(l),:);
        g{l,6} = ptint(allcodes==fun(l),:);
        g{l,7} = ptintabs(allcodes==fun(l),:);
        clusnum = copynumfinal{fun(l),i+1};
        if size(g{l,2},1) < clusnum
            clusnum = size(g{l,2},1);
        end
            bobobo = 1;
            keepgoing = 1;
            while keepgoing == 1
                if clusnum <2
                    C = mean(g{l,2},1);
                    idx = ones(size(g{l,2},1),1);
                else
                    [idx,C] = kmeans(g{l,2},clusnum);
                end
                for m = 1:size(C,1); 
                    D{m,:} = pdist2(C(m,:),g{l,2}(idx == m,:)); 
                end
                s = sum(cell2mat(cellfun(@(x) sum(x > sqrt(radius))>0, D, 'UniformOutput',0)));
                clear D;
                if s > 0
                    clusnum = clusnum + 1;
                    copynumfinalrevised{fun(l),i+1} = clusnum;
                elseif pdist(C)<sqrt(radius)+.00001
                    clusnum = clusnum - 1;
                    copynumfinalrevised{fun(l),i+1} = clusnum;
                else
                    keepgoing = 0;
                end
                bobobo = bobobo + 1;
                if bobobo == 100;
                    keepgoing = 0;
                end
            end
        g{l,5} = C;
        PosList{fun(l),i+1} = C; 
    end
    dotlocations(i).cell = g;
end
