function y = WC_downsampling(x,bin)

n = size(x,1);

if mod(size(x,2),bin) == 0
   T = size(x,2)/bin;
else
   T = fix(size(x,2)/bin)+1;
end

y = zeros(n,T);

for i = 1:T
    y(:,i) = x(:,bin*(i-1)+1);
end

end