function cost= compute_cost(k,m, tspan, x0, pe_ref, data) 

pe = pe_ref; %here define pe using pe_ref and k 
pe([9 10]) = k';

for i=1:length(m)
    pe(8) = m(i); %here pe stores the ith inducer concentration
    [t x]= ode15s(@you_odeRI,tspan, x0, [], pe);
    sim(i,1)= x(end,3); %store the value at steady state of a chosen variable (either N, E, or A) 
end
cost= sqrt(sum(( sim - data ).^2)/length(m)); %here you need to define the cost
