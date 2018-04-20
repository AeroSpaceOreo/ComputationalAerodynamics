%Function for Thomas Algorithm, bi-diagonal 
function ans = Talgo(u,d,s) %u for upper-diag, d for cen-diag, s for sol

n = max(size(s)); %To get the length of the solution matrix

%[d1 u1                            ][ phi1 ] [x1^2]
%[l2 d2 u2                         ][ phi2 ] [x2^2]
%[   l3 d3 u3                      ][  .   ] [ .  ]
%[           ...                   ][  .   ]=[ .  ]
%[              ln-1 dn-1 un-1     ][  .   ] [ .  ]
%[                    ln   dn   un ][ phin ] [ .  ]
%[                        ln+1 dn+1][phin+1] [ .  ]

for i = 2:1:n
    d(i) = d(i)-(u(i-1)*u(i-1)/d(i-1));
    s(i) = s(i)-(s(i-1)*u(i-1)/d(i-1));
end

yn(n) = s(n)/d(n); %Divide by center diagonal to normalize
for j = n-1:-1:1
    yn(j)=(s(j) - u(j)*yn(j+1))/d(j);
end

ans = yn; %return the value


