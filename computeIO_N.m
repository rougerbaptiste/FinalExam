function Neq= computeIO_N(m, tspan, x0, pe) 
%Neq is a matrix that stores the cell density at steady state for all circuits and all inducer concentrations.
%More precisely, Neq(i,j) stores the cell density obtained for circuit j and inducer concentration m(i)
%m is a vector of inducer concentrations

for i=1:length(m)
    pe(8) = m(i); % update pe using m(i)
    [t x]= ode15s(@you_ode,tspan, x0, [], pe); % num sim of you circuit
    Neq(i,1)= x(end,1); %stores N at steady state
    [t x]= ode15s(@you_odeR,tspan, x0, [], pe); % num sim of youR circuit
    Neq(i,2)= x(end,1); %stores N at steady state
    [t x]= ode15s(@you_odeI,tspan, x0, [], pe);% num sim of youI circuit
    Neq(i,3)= x(end,1); %stores N at steady state
    [t x]= ode15s(@you_odeRI,tspan, x0, [], pe); % num sim of youRI circuit
    Neq(i,4)= x(end,1); %stores N at steady state
end
