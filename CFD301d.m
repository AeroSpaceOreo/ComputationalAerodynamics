%Student ID: 1001358558
%Name: Ching Huai Wang

clc;
clear;

%%Given condition
c = 0.1
x_min = 0;
x_max = 100;
t_max = 100; %Use 20/100/800s


%Discretization
n = 250; %number of cells
delx = (x_max-x_min)/n; %delta x
x = x_min-(2*delx):delx:x_max;

t = 0 %Starting from time = 0
delt = 2

lambda = 2 %CFL number

%%Boundary Conditions
%phi(0,t) = 1
phi0t = 1;
%phi(x,0) = 1/(1+exp(x-10))
phix0 = 1./(1+exp(x-10));

%Plot
plot(x,phix0)

phi = phix0;
phi_step = 1./(1+exp(x-10));

timestep = t_max/delt %number of time steps

for i = 1:timestep
    phi(3) = phi0t;
    phi(1) = phi(3);
    phi(2) = phi(1);
    
    for j = 3:n+1
        %The equation derived in homework given by
        phi_step(j) = phi(j)-...
                      0.5*lambda*(phi(j-2) -4*phi(j-1)+3*phi(j))+...
                      0.5*(lambda^2)*(phi(j-2)-2*phi(j-1)+phi(j));
    end
    
    t = t+delt; %adding delta(t) every loop
    phi = phi_step;
    phi_exact = 1./(1+exp((x-10)-(c*t)));
    
    %Plot the function
    plot(x,phi_exact,'b-')
    hold on
    grid on
    plot(x,phi,'r-')
    
    
end
