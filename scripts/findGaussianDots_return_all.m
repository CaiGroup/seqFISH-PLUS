function [gXf] = findGaussianDots_return_all(locDots,locIm,Pixels)

    options = optimset('Display','off');
    Xf = [];
    sigma = [];

    for dotID = 1:size(locDots,1)
        data = double(locIm([-Pixels:Pixels]+locDots(dotID,2),[-Pixels:Pixels]+locDots(dotID,1)));
        x0 = [min(min(data)) max(max(data))-min(min(data)) Pixels+1 Pixels+1 1 1 45];
        sup = size(data);
        f = @(x)D2gauss_wRot(x,sup,data);
        [xf, fva] = fminsearch(f,x0,options);
        Xf = [Xf; xf];

    end
    gXf = Xf;
    gXf(:,3:4) = locDots(:,1:2) + Xf(:,[4 3]) - [Pixels, Pixels]-1;
end
