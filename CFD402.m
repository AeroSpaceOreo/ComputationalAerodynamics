%Computational Aerodynamics HW402-Exact Solution
%Student ID: 1001358558
%Name: Ching Huai Wang

clear;
clc;

grid = 5; %The number of grids
mn = 15; %Summation number from 1:mn to compute
xgrid = linspace(0,1,grid); %0<x<1
ygrid = linspace(0,1,grid); %0<y<1
    
%Create the mesh of grid, for MATLAB to compute
[x,y] = meshgrid(xgrid,ygrid);

for m = 1:mn
    for n = 1:mn
for  i = 1:grid
    for j = 1:grid
        %The Amn given
        A(i,j) = (4/(n*m*pi))*(cos(m*pi)*cos(n*pi)-cos(m*pi)-cos(n*pi)+1);
        %The exact solution given
        u(i,j) = -(A(i,j)/(pi^2*(n^2+m^2)))*sin(n*pi*x(i,j))*sin(m*pi*y(i,j));
    end
end
    end
end

%Plot 2-D Contour
contour3(x,y,u,100)
xlabel('X Direction')
ylabel('Y Direction')
title('Exact Solution')