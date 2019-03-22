function [ cellin ] = removezeros2( cellin )
%REMOVEZEROS2 Summary of this function goes here
%   Detailed explanation goes here
log = cellin == 0;
cellin(log) = [];

end

