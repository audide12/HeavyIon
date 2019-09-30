%i m using this at the moment for white noise
tic
datestr(now,'HH:MM:SS')

Tf = 0.760148;tauf = 6.35185;d = 2;
chi = (2/3)*Tf^2;
sf = 9.1518;
%m = 2.5018;%Kaon
%m = 4.75853;%proton
m = 0.707292;%pion
w = 1.55;
a1=chi*Tf/(tauf*sqrt(3.14*2.38));

M = 2000;Xi = 18;dxi = Xi/M; xi = -(Xi/2):dxi:(Xi/2);
Gauss=a1*exp(-((xi+0.0)/1.55).^2);Gauss = Gauss';
Xinew = Xi/2;xinew = -(Xinew/2):dxi:(Xinew/2);
Kron = zeros(M+1,1);Kron(length(xinew),1) = 1*chi*Tf/(tauf*dxi);
Cov=Kron - Gauss;

M1 = 2000;Y = 5;Deltady = Y/M1;y1=0;

Prefactor = (chi*Tf)/sqrt(w^2*pi*tauf^2);
C = zeros(M1+1,1);

%Cov = COV14;tauQ=0.48;Cov(M+1,1)=0;

parfor dyi = 1:M1+1
    
    y2 = 0 + (dyi*Deltady);
    for xi1 = 1:length(xinew)
       
        Xi1 = -(Xinew/2) + (xi1*dxi);
        Fn1 = (1/(chi*cosh(y1-Xi1)^2))*gamma(3)*(1-gammainc((m/Tf)*cosh(y1-Xi1),3));
        
        for xi2 = 1:length(xinew)
           
            Xi2 = -(Xinew/2) + (xi2*dxi);
            Fn2 = (1/(chi*cosh(y2-Xi2)^2))*gamma(3)*(1-gammainc((m/Tf)*cosh(y2-Xi2),3));
            x = round(abs(Xi1-Xi2)/dxi)+length(xinew);
            C(dyi,1) = C(dyi,1) + dxi*dxi*Fn1*Fn2*Kron(x,1);
       end
        
    end
    
end
%Qfactor = 0.0989996;%proton
%Qfactor = 0.912565;%kaon
Qfactor = 3.33243;%pion

C = ((d*tauf*Tf)/(4*pi^2*Qfactor)).*C;
y = 0:Deltady:Y;

%C = (tauf-tauQ).*C; %unexplained required factor
%figure;scatter(y,C);
%figure;scatter(y,C,'Marker','d');hold on;scatter(A1,B1,'Marker','o')

toc