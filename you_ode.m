function xdot= you_ode(t,x,p) 

%x_names= {'N','E','A'}; -> x(1) is N, x(2) is E, x(3) is A 
%p_names= {'k','Nm','d','ke','de','va','da'}; -> p(1) is k, p(2) is Nm, etc

xdot= zeros(size(x));

xdot(1) = p(1)*x(1) * (1-x(1)/p(2)) - (p(3)*x(2)*x(1));
xdot(2) = p(4)*x(3) - p(5)*x(2);
xdot(3) = p(6)*x(1) - p(7)*x(3);
end
