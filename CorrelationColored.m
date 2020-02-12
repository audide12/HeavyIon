%cl=parcluster('local')
%pool = cl.parpool()

tic
datestr(now,'HH:MM:SS')

T = 6.352; N = 2^12; taui=0.5;L=N;
D = 0.162035; tauQ=0.48;
sf = 9.1518;
chi = (2/3)*Tf^2;
dt = (T-taui)/N; 
dx = (T-taui)/(N*taui);
Tf = 0.760148;tauf = T;taui = 0.5;
w2 = 8*D*((1/taui)-(1/tauf));
b=taui/tauQ;vQ2=D/tauQ;

Xi=18;M=200;dxi=Xi/M;

es = 5.5041e-04;
n=N;
k=(M/2)+1;
Cov4 = zeros(M,1);
J = 10000000;

parfor q = 1:J
    
F = zeros(M,N);

randn('state',q)
F(:,L) = sqrt(es*2*D/(2*tauQ*dxi))*randn(M,1);

randn('state',I*q+J)
dW = sqrt(es*2*D/(dt*dxi))*randn(M,N);

for e = L:-1:2
     
   F(:,e-1) = F(:,e) - dx*(F(:,e)-dW(:,e))*b;   
  
end

Z = zeros(M,L);

%SDE  Solver 
for i = 3:1:L
    
    x=(i*dx)+1;  
    i0 = i-1;i2 = i+1;
    if i0 == 0 
            i0=L;
    elseif i2 == L+1 
            i2=1;
    end
    
    for j=1:1:M  %this is the xi movement
    
        u0=j-1;u1=j+1;
        if u0 == 0     
             u0=M;
        elseif  u1==(M+1)
            u1=1;
        end
        
        Zdxidxi = (Z(u1,i)+Z(u0,i)-2*Z(j,i))/(dxi*dxi);
        fxi = (F(u1,i)-F(j,i))/(1.0*dxi);
        fxxi = (F(u1,i)+F(j,i0)-F(u1,i0)-F(j,i))/(dx*dxi);
         
        comp1 = (1/dx)+((x*tauQ)/(dx^2*(taui*x+2*tauQ)));
        comp2 = Z(j,i)*((1/dx)+((2*x*tauQ)/(dx^2*(taui*x+2*tauQ)))) - Z(j,i0)*((x*tauQ)/(dx^2*(taui*x+2*tauQ))) + Zdxidxi*(D/(x*(taui*x+2*tauQ))) - fxi*((taui*x+tauQ)/(x*(taui*x+2*tauQ))) - fxxi*(tauQ/(taui*x+2*tauQ));
         
        Z(j,i2) = comp2/comp1; 
        
    end
   
   
end

Cov4 = Cov4 + (Z(k,n).*Z(:,n));

end

Cov4=Cov4./(J);

z=taui:dt:T-dt;
xi = -Xi/2:dxi:(Xi/2)-dxi;
Cov = (sf^2).*Cov4;

toc
datestr(now,'HH:MM:SS')
