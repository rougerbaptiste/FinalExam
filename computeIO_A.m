function Aeq= computeIO_A(m, tspan, x0, pe) 
% see computeIO_N for comments

for i=1:length(m)
    pe(8) = m(i); % update pe using m(i)
    [t x]= ode15s(@you_ode,tspan, x0, [], pe); % num sim of you circuit
    Aeq(i,1)= x(end,3); %stores E at steady state
    [t x]= ode15s(@you_odeR,tspan, x0, [], pe); % num sim of youR circuit
    Aeq(i,2)= x(end,3); %stores E at steady state
    [t x]= ode15s(@you_odeI,tspan, x0, [], pe); % num sim of youI circuit
    Aeq(i,3)= x(end,3); %stores E at steady state
    [t x]= ode15s(@you_odeRI,tspan, x0, [], pe); % num sim of youRI circuit
    Aeq(i,4)= x(end,3); %stores E at steady state
end
