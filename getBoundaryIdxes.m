function[idxes] = getBoundaryIdxes(sizeA)
%% boundary indexes for 2D matrix
% sizeA: self-explained
%% Example
% getBoundaryIdxes([3, 3])
% 
% ans =
% 
%      1     4     7     3     6     9     2     8
    
    m = sizeA(1);
    n = sizeA(2);
    num = 2 * (m + n) - 4;
    
%     low up row
    r1 = ones(1, n);
    c1 = 1:n;
%     left right col
    r2 = 2:m-1;
    c2 = ones(1, m-2);
    
    idxes = sub2ind(sizeA, [r1, r1*m, r2, r2], [c1, c1, c2, c2*n]);
end