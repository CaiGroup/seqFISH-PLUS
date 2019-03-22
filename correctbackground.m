function hybnumcor = correctbackground(hybnum, corrections,channels)

    for j = 1:length(channels)
        B = repmat(corrections.corrections{j},1,1,size(hybnum.color{j},3));
        hybnumcor.color{j} = uint16(double(hybnum.color{j})./B);
    end