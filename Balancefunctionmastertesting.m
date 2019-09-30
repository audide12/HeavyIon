%Balance Function master -- I don't think it gets better than this.
tic
datestr(now,'HH:MM:SS')
Tf = 0.760148;tauf = 6.3244;d = 2;
chi = (2/3)*Tf^2;
sf = 9.1518;
%m = 2.5018;%Kaon
%m = 4.75853;%proton
m = 0.707292;%pion
w = 1.55;
a1=chi*Tf/(tauf*sqrt(3.14*2.38));

M = 300;Xi = 18;dxi = Xi/M; XI1 = -9.0:dxi:9.0;
XI2 = -20:dxi:20;
xi = -29:dxi:29;

M1 = 2000;Y = 5;Deltady = Y/M1;y1=0;

Prefactor = (chi*Tf)/sqrt(w^2*pi*tauf^2);
C = zeros(M1+1,1);
COVE = zeros(1,length(xi));
COVE(1,((length(xi)+1)/2)-(M/2):((length(xi)+1)/2)+(M/2)-1) = Cov(:,1);%this one
Cov = COVE';
%tauQ=4.0;
%
parfor dyi = 1:M1+1
    Ddy = 0 + (dyi*Deltady);
    Sum1 = 0;
    for xi1 = 1:1:length(XI1)   
        Xi1 = XI1(xi1);
        Sum2 = 0;
        for xi2 = 1:1:length(XI2)
            Xi2 = XI2(1,xi2);
            x = round(abs(Xi1-Xi2)/dxi)+(((length(xi)+1)/2));%this one
            Fn1 = (1/(chi*cosh(y1-Xi1)^2))*gamma(3)*(1-gammainc((m/Tf)*cosh(y1-Xi1),3));
            Fn2 = (1/(chi*cosh(Ddy-Xi2)^2))*gamma(3)*(1-gammainc((m/Tf)*cosh(Ddy-Xi2),3));
            %Fn2 = (1/(chi*cosh(Ddy-Xi1)^2))*gamma(3)*(1-gammainc((m/Tf)*cosh(Ddy-Xi1),3));
            %Fn1 = (1/(chi*cosh(-Xi2)^2))*gamma(3)*(1-gammainc((m/Tf)*cosh(-Xi2),3));
            
            Sum2 = Sum2 + (Cov(x,1)*Fn1*Fn2*dxi);
            
        end
        Sum1 = Sum1 + (Sum2*dxi);
    end
    C(dyi,1) = Sum1;
end
%Qfactor = 0.0989996;%proton
%Qfactor = 0.912565;%kaon
Qfactor = 3.33243;%pion

C = ((d*tauf*Tf)/(4*pi^2*Qfactor)).*C;
y = 0:Deltady:Y;

%figure;scatter(y,C);
%figure;scatter(y,C,'Marker','d');hold on;scatter(A1,B1,'Marker','o')
%}
%length(xi)/2 = 1389
%M/2 = 250
%COVE15(1,1389-250:1389+249) = COV15(:,1);
toc