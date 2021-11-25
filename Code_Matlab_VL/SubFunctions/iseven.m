function bool = iseven(X)
% test if number is even (1) or not (0)
% non-integers are not even

bool = isequal(X/2,floor(X/2));


end