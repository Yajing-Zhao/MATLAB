
function outf = fdextension(k,j,xbar,u)

n = length(j);
if k>=n 
    error('***  length(j) must be larger than k');
end

c = fdcoeffF(k,0,j);      % coefficients for k'th derivative 

    function out = computeh(h)
        out = 0;
        for i=1:n
            out = out + c(i) * u(xbar + j(i) * h);
        end
        out = out  / (h.^k);
    end

outf = @computeh;
end 
