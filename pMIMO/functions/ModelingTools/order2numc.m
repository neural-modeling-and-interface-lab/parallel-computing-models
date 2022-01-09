function D = order2numc(L, ModelOrder)
% D = order2numc(L, ModelOrder)
%
% if ModelOrder==1,
%     D = L;
% elseif ModelOrder==2,
%     D = L+L*(L+1)/2;
% elseif ModelOrder==3,
%     D = L+L*(L+1)/2+L*(L+1)*(L+2)/6;
% else,
%     error('Invalid Model Order');
% end;

if ModelOrder == 0 % Modified by Xiwei Dec.13, 2021
    D = 0;
elseif ModelOrder==1
    D = L;
elseif ModelOrder==2
    D = L+L*(L+1)/2;
elseif ModelOrder==3
    D = L+L*(L+1)/2+L*(L+1)*(L+2)/6;
else
    error('Invalid Model Order');
end