function y = CosInt(A,wa,wb)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if A==0
    y=(1/pi)*(wb-wa);
else
    y=1/(pi*A)*(sin(wb*A)-sin(wa*A));
end

