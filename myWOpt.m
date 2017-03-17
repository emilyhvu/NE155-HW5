function [w_opt]=myWOpt(n, error,tol)

w_v=[1:0.01:1.99];
iter_v=[];
for t=1:numel(w_v)
    w=w_v(t);
    [SOR, SOR_iter]=mySOR(n,error,tol,w);
    iter_v=[iter_v SOR_iter];
end

w_opt_index=find(min(iter_v)==iter_v);
w_opt=w_v(w_opt_index);

plot(w_v,iter_v,'r+')
axis([1 1.5 0 30])
title('SOR Iterations vs. w')
xlabel('w')
ylabel('iterations')

end
