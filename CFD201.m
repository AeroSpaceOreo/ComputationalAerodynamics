%Computational Aerodynamics HW#2
%Name: Ching Huai Wang; Student ID: 1001358558
clc;
clear;

n = 250; %cells

%%-----Solve Equation-----
%First we solve the differential equations using dsolve function
syms phi(x) dphi(x) ddphi(x) x %declaring variables phi and x
dphi = diff(phi,x); %1st derivative of phi
ddphi = diff(dphi,x); %2nd derivative of phi
%using dsolve function solving the following diff eqn
%ddphi - 0.59*phi = x^2
%boundary conditions phi(0) = 0, phi(1) = 0
phi = dsolve(ddphi-0.5*phi == x^2, phi(0) == 0, phi(1) == 0)
l = linspace(0,1,n+1); %n+1 rows of solution
%the exact solution of y
yexact = single(subs(phi,x,l));

%%-----Create Matrix
%Tx=r
%u for upper, l for lower, d for center diagonal
%[d1 u1                            ][ phi1 ] [x1^2]
%[l2 d2 u2                         ][ phi2 ] [x2^2]
%[   l3 d3 u3                      ][  .   ] [ .  ]
%[           ...                   ][  .   ]=[ .  ]
%[              ln-1 dn-1 un-1     ][  .   ] [ .  ]
%[                    ln   dn   un ][ phin ] [ .  ]
%[                        ln+1 dn+1][phin+1] [ .  ]

%Since at the same row, l will be same as u

%for i=1 l_1=u_1, l_1 is out of column
%l(1,1) = 0; %doing bi-->tridiagonal, only need l or u
u(1,1) = 0;

%for i=n+1 l_n+1=u_n+1, u_n+1 is out of column
%l(n+1,1) = 0; %only need l or u
u(n+1,1) = 0; 

%Fill in the rest of the column
%The rest value lower and upper diagonal is the same 1/delta(x)
delx=1/n;
coeff = 1/(delx^2);
%l(2:n,1) = coeff; %only need l or u
u(2:n,1) = coeff;

%The diagonal matrix d
d(1,1) = coeff; d(n+1,1) = coeff;
beta = 2+(0.5*(delx^2));
coeff2 = beta/(delx^2);
d(2:n,1) = -1*coeff2;

%The solution matrix s, which is a column of x^2
s(1,1) = 0; s(n+1,1) = 0; %As the boundary condition states
for i = 0:1:(n-2)
    s(2+i,1) = (delx+(delx*i))^2;
end

%%Show Matrix
%Print out the matrix to see if it's correctly set
T = zeros(n+1,n+1); %Create zeros for T matrix
for k = 0:1:n
    T(1+k,1+k) = d(1+k,1); %Set center diagonal
    T(1+k,2+k) = u(1+k,1); %Set upper diagonal
end
disp('PRINT MATRIX T(BI-DIAGONAL), displaying 6x6 only:')
disp(T(1:6,1:6))

%%Thomas Algorithm
for j = 1:n
    %Divide by d to normalize
    u(j) = u(j)/d(j);
    s(j) = s(j)/d(j);
    d(j) = 1; %d(j)/d(j)=1, divide by itself
    %Setting the pivot element
    alpha = u(j+1);
    d(j+1) = d(j+1)-alpha*u(j);
    s(j+1) = s(j+1)-alpha*s(i);
end

%%Back solve of Thomas Algorithm
%y solution matrix
ynum = zeros(n+1,1);
ynum(n+1,1) = s(n+1)/d(n+1); %normalize the n+1 element
for h = n:-1:1
    ynum(h) = s(h)-u(h)*ynum(h+1);
end

plot(yexact)

%%Plot Solution
figure
plot(l,yexact,'k-')
hold on
plot(l,-ynum,'b-')
xlabel('Length')
ylabel('Black(exact) / Blue(approximate)')

%%Plot of L^2 Norm
sqr = 2:n+1;
dist = yexact-ynum;
LSNorm = norm(dist,2);
%Plot
figure
loglog(LSNorm)
xlabel('Least Square Norm')