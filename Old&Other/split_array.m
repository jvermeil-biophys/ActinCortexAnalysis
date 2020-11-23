function matrix = split_array(array)
%This function splits an array into a number of colums. 
%This by creating a new row every time there is a large step in
%the values in the array. This function assumes the values to be sorted!
%Input:     -Array containing SORTED (intensity) values
%Output:    -Matrix created from the original array

array = double(array);
lengte = length(array);
previous_value = array(1);
teller2 = 1;
teller3 = 1;
for teller = 1:lengte
    if ((array(teller) == previous_value) | (array(teller) == (previous_value + 1)))
        matrix(teller3, teller2)=array(teller);
        teller = teller + 1;
        teller2 = teller2 + 1;
    else
        matrix(teller3 + 1, 1)=array(teller);
        teller = teller + 1;
        teller2 = 2;
        teller3 = teller3 + 1;
    end
    previous_value = array(teller - 1);
end
%matrix = unint8(matrix);
