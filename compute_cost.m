function cost= compute_cost(k,m, tspan, x0, pe_ref, data) 

.. %here define pe using pe_ref and k 
..

for i=1:length(m)
    .. %here pe stores the ith inducer concentration
    [t x]= ode15s(..);
    sim(i,1)= .. ; %store the value at steady state of a chosen variable (either N, E, or A) 
end
cost= .. %here you need to define the cost
