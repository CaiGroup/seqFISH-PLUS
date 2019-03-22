function [FISH_only] = removeBeads(allDots,beadDots)

    refDots = beadDots;
    locDots = allDots;

    my_color = [0 0 1];
    myDots = [];
    badDots = [];

    % remove ALL beads, using bead only image as reference
    for dot = 1:length(refDots)
        currDot = refDots(dot,:);
        locDist = [];
        for i = 1:length(locDots)
            locDist(i) = norm(currDot-locDots(i,1:2));
        end
        [m, I] = min(locDist);
        I = find(locDist<1);
            myDots = [myDots;locDots(I,1:2)];
            badDots = [badDots, I];

    end

    locDotsX = locDots(:,1);
    locDotsY = locDots(:,2);
    IDs = locDots(:,3);
    locDotsX(badDots) = [];
    locDotsY(badDots) = [];
    IDs(badDots) = [];
    FISH_only = {[locDotsX locDotsY IDs]};

end