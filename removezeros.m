function cellin = removezeros(cellin)

log = cellin == 0;
cellin(log) = [];
if isempty(cellin) ==1
    cellin = 0;
end

