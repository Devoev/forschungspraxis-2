function [d_max,d_min] = calc_max_min_pos(y,L,d,lambda)
    j=1;
    for k = ceil(min(y)/(lambda*(L/d))):floor(max(y)/(lambda*(L/d)))
        d_max(j) = (k*lambda)*(L/d);
        d_min(j) = ((2*k+1)*lambda/2)*(L/d);
        j = j+1;
    end
end