function msd= compute_msd(m,tspan,x0,pe, pe_nom,circuit,EorA)

%circuit: 1, 2, or 3 corresponds to R, I or RI
%EorA: 1 or 2 corresponds to E or A

for i=1:length(m)
    pe(8) = m(i); %update pe
    pe_nom(8) = m(i); %update pe_nom
    if circuit==1
        [t x_nom]= ode15s(@you_odeR,tspan, x0, [], pe_nom); %compute nominal behavior of circuit R
        [t x]= ode15s(@you_odeR,tspan, x0, [], pe); %compute behavior of circuit R for specific values of theta and eta, stored in pe
    elseif circuit==2
        [t x_nom]= ode15s(@you_odeI,tspan, x0, [], pe_nom); %the same for circuit I
        [t x]= ode15s(@you_odeI,tspan, x0, [], pe);
    elseif circuit==3
        [t x_nom]= ode15s(@you_odeRI,tspan, x0, [], pe_nom);%the same for circuit RI
        [t x]= ode15s(@you_odeRI,tspan, x0, [], pe);
    end
    
    if EorA==1
        output_nom(i)= x_nom(end, 2); %stores the steady state value of E in nominal case
        output(i)= x(end, 2); %stores the steady state value of E in specific case;
    else
        output_nom(i)= x_nom(end,3); % same but stores the value of A
        output(i)= x(end,3);
    end
end
msd= sqrt( (sum( (output - output_nom).^2))/length(m) ); %here is the square root of the sum of squared differences between the output and the nominal output, divided by the number of values
end