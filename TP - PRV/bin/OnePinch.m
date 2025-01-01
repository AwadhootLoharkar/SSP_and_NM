function [y] = OnePinch(x, p, z)
    if (x(end) <= p)
        error("The pinch have to be performed on the string!")
    end
    for i=1:length(x) % eq 2.63
        if(x(i) < p)
            y(i) = x(i)*z/p;
        else
            y(i) = z*(x(i)-x(end))/(p-x(end));
        end
    end
end