%had an extra initial condition for adhesion
function u0 = pdeic_two_species_1D(x)

%global max(x) 
global x_max;


%initialise
u0=[0;0];

if x >=0.4*x_max && x<=0.6*x_max
    u0(1)=1;
    u0(2)=0;
    
else
    u0(1)=0;
    u0(2)=0.5;
end

