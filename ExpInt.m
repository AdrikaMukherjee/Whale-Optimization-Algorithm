function y = ExpInt(p,q,n,wa,wb)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if(p ==0)
    y1=(wb-wa)/(2*pi);
else
    y1=1/(2*pi*p)*(sin(wb*p)-sin(wa*p));
end


if(q==n-1)
    y2=(wb-wa)/(2*pi);
else
    y2=1/(2*pi*(q-n+1))*(sin(wb*(q-n+1))-sin(wa*(q-n+1)));
end
y=y1+y2;
    
 end

