%Computational Aerodynamics HW401-Numerical Solution + Iteration Plot
%Student ID: 1001358558
%Name: Ching Huai Wang

clear;
clc;

nn = [5 25 50];

for p = 1:3; %Number of grid
n = nn(1,p);
    
delx = 1/n;

xgrid = linspace(0,1,n); %0<x<1
ygrid = linspace(0,1,n); %0<y<1
    
%Create the mesh of grid, for MATLAB to compute
[x,y] = meshgrid(xgrid,ygrid);

tol = 1e-5; %The tolerance of the system
%The error, we are going to compute when error is greater than tolerance
E = 1; %Error

k = 0; %For 1st iteration, next iteration is k+1, so on and so for

phi = zeros(n,n);

%Condition to compute
while E>tol
    k = k+1; %Next iteration
    pphi = zeros(n,n); %Dimension matrix of the solution will be nxn
    
for j = 2:n-1
    
    for i = 1:n-2
        %Making of the coefficient matrix as follows
        %[ 1 -4  1                         ]
        %[    1 -4  1                      ]
        %[       1 -4  1                   ]
        %[           ...                   ]
        %[              ln-1 dn-1 un-1     ]
        %[                    ln   dn   un ]
        %[                        ln+1 dn+1]
        %Only using bi-diagonal matrix(neglect the l)
        u(i) = 1; %Upper diagonal
        d(i) = -4; %Center diagonal
        s(i) = delx^2-(phi(i+1,j-1)+phi(i+1,j+1)); %Solution matrix      
    end
    
    %Modify the solution matrix
    s(1) = 0.04 - (phi(i+1,j+1)+phi(i+1,j-1)) - phi(1,j);
    s(n-2) = 0.04 - (phi(i+1,j+1)+phi(i+1,j-1)) - phi(n,j);
    
    %Run the Bi-Diag Thomas Algorithm function, please see Talgo.m
    A = Talgo(u,d,s);
    %Changing the new value(excluding B.C.s in the solution)
    pphi(2:n-1,j) = A; %Replace from the 2 to (end-1)   
end 
 
%Adding the errors
summation = sum(sum((pphi-phi).^2));
E = sqrt(summation);

%Assign the new phi to old phi
phi = pphi;
end

%Plot 2-D Contour of Approximate Solution
figure
contour(x,y,phi,n*5)
xlabel('X direction')
ylabel('Y direction')
title('Approximate Solution')

%Store the iteration number k in a matrix every loop
itr = p+1;
kmat(itr) = k;
end

%Since by the storage of k matrix, we get one more column
select = size(nn);
select = select(1,2);

%Plot Iteration to Grid
figure
kmat = kmat(2:select+1); %Excluding kmat(1,1) = 0
plot(kmat,nn)
xlabel('Iterations')
ylabel('Grid')
title('Iteration to Grid Plot')