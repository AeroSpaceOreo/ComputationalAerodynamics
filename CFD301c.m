%Amplification Factor |G|
%Student ID: 1001358558
%Name: Ching Huai Wang

clc;
clear all;
%Solving for dphi/dt+c*(dphi/dx) = 0

%From the solution I solved for question (b) of  homework
%|G| = e^(at) = 1-(lambda/2)*...please refer to my derived equation
%lambda is the notation of the CFL number


theta = linspace(0,pi,100); %theta for angle matrix
i = sqrt(-1)

%Given Values
for lambda = 0.6
    G = 1-...
    (lambda/2)*(3-4*(cos(theta)-i*sin(theta))+2*(cos(2*theta)-i*sin(2*theta)))+...
    ((lambda^2)/2)*(1-2*(cos(theta)-i*sin(theta))+(cos(2*theta)-i*sin(2*theta)));
    
    modG = abs(G)
    
    plot(theta,modG) %Plot graphs
    grid on
    hold on %We are going to combine all graphs
     
end

hold off
xlabel('theta')
ylabel('|G|')

legend
