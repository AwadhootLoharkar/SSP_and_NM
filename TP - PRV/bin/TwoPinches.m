function [y] = TwoPinches(x, p1, p2, z1, z2)
    if (z1 > z2 && x(end) <= p2)
        error("The pinches have to be performed on the string and p1 < p2");
    end
    for i=1:length(x) % eq 2.64
        if(x(i) < p1)
            y(i) = x(i)*z1/p1;
        else if(x(i) < p2)
            y(i) = (z1-z2)/(p1-p2)*(x(i)-p1)+z1;
        else
            y(i) = z2*(x(i)-x(end))/(p2-x(end));
        end
    end
end