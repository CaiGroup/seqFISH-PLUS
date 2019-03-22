function vertex = selfseg(PathName)

files = dir([PathName filesep '*.roi']);
[sROI] = ReadImageJROI(strcat(PathName, filesep, {files.name}));


for i = 1:length(sROI)
    if strcmp(sROI{1,i}.strType, 'Oval')
        y0 = mean(sROI{1,i}.vnRectBounds(logical([1 0 1 0])));
        x0 = mean(sROI{1,i}.vnRectBounds(logical([0 1 0 1])));
        rx = abs(sROI{1,i}.vnRectBounds(logical([0 0 1 0])) - sROI{1,i}.vnRectBounds(logical([1 0 0 0])))/2;
        ry = abs(sROI{1,i}.vnRectBounds(logical([0 1 0 0])) - sROI{1,i}.vnRectBounds(logical([0 0 0 1])))/2;
        t = linspace(0, 2*pi, 360)';
        vertex(i).x = x0 + ry * cos(t);
        vertex(i).y = y0 + rx * sin(t);
    elseif strcmp(sROI{1,i}.strType, 'Rectangle')
        vertex(i).y = [sROI{1,i}.vnRectBounds(1) sROI{1,i}.vnRectBounds(3) sROI{1,i}.vnRectBounds(3) sROI{1,i}.vnRectBounds(1)];
        vertex(i).x = [sROI{1,i}.vnRectBounds(2) sROI{1,i}.vnRectBounds(2) sROI{1,i}.vnRectBounds(4) sROI{1,i}.vnRectBounds(4)];
    else
        vertex(i).x = sROI{1,i}.mnCoordinates(:,1);
        vertex(i).y = sROI{1,i}.mnCoordinates(:,2);
    end
end


