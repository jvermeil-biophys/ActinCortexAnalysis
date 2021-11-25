function bool = isodd(X)
% test if number is odd (1) or not (0)
% non-integers are not odd

bool = isequal((X+1)/2,floor((X+1)/2));


end