function solver(X, Y, four)


%---------------------константы--------------------%
%Задание параметров Ламе через технические константы
E=2e+11;
w=0.33;
Ro = 7800;
alpha = 11.7e-6;
T = 0.001; %толщина
%Задание матрицы упругих констант
h=(E*w)/((1+w)*(1-2*w));
m=E/(2*(1+w));
De = [h+2*m h 0;h h+2*m 0;0 0 m];
%---------------------------------------------------%

g = - 0;%input("ускроние свободного падения: ");
dT = 0;%input("температурное изменение: ");
q = - 100;%input("распределенная нагрузка: ");

e0 = [alpha*dT; alpha*dT; 0];
Kg = zeros(2*max(max(four)));
Fg = zeros(2*max(max(four)),1);

c = [0.555555556 0.888888889 0.555555556];
z = [-0.774596669 0.000 0.774596669];

for p =1:size(four,1)
    Ke = zeros(8);
    Fet = zeros(8,1);
    Feg = zeros(8,1);
    Feq = zeros(8,1);

r = transpose([X(four(p,1)) X(four(p,2)) X(four(p,3)) X(four(p,4));Y(four(p,1)) Y(four(p,2)) Y(four(p,3)) Y(four(p,4))]);

Ae = (1/2)*(X(four(p,1))*Y(four(p,2))+X(four(p,2))*Y(four(p,3))+X(four(p,3))*Y(four(p,4))+X(four(p,4))*Y(four(p,1))-X(four(p,2))*Y(four(p,1))-X(four(p,3))*Y(four(p,2))-X(four(p,4))*Y(four(p,3))-X(four(p,1))*Y(four(p,4)));

for i=1:3
    for j=1:3
        xi = z(i);
        eta = z(j);

        N = transpose([0 1/4*(1-xi)*(1-eta) 0 1/4*(1+xi)*(1-eta) 0 1/4*(1+xi)*(1+eta) 0 1/4*(1-xi)*(1+eta)]);
        dN = [1/4*(-1+eta) 1/4*(1-eta) 1/4*(1+eta) 1/4*(-1-eta); 1/4*(-1+xi) 1/4*(-1-xi) 1/4*(1+xi) 1/4*(1-xi)];

        J = dN*r;
        dN1 = inv(J)*[dN(1,1);dN(2,1)];
        dN2 = inv(J)*[dN(1,2);dN(2,2)];
        dN3 = inv(J)*[dN(1,3);dN(2,3)];
        dN4 = inv(J)*[dN(1,4);dN(2,4)];

        Be1 = [dN1(1) 0; 0 dN1(2);dN1(2) dN1(1)];
        Be2 = [dN2(1) 0; 0 dN2(2);dN2(2) dN2(1)];
        Be3 = [dN3(1) 0; 0 dN3(2);dN3(2) dN3(1)];
        Be4 = [dN4(1) 0; 0 dN4(2);dN4(2) dN4(1)];
        Be = [Be1 Be2 Be3 Be4];
        BeT = transpose(Be);
        Ke = Ke + c(i)*c(j)*BeT*De*Be*det(J)*T;

        Fet = Fet + c(i)*c(j)*BeT*De*e0*det(J)*T;
        Feg = Feg + c(i)*c(j)*Ro*g*det(J)*T*N;
    end
end

delta = Ae/1e+6;

if and(abs(Y((four(p,1)))+sqrt(3)*X((four(p,1))) - tan(pi/3)*4) < delta, abs(Y((four(p,2)))+sqrt(3)*X((four(p,2))) - tan(pi/3)*4) < delta)
    eta = -1;

    for i=1:3
        xi = z(i);

        N = transpose([0 1/4*(1-xi)*(1-eta) 0 1/4*(1+xi)*(1-eta) 0 1/4*(1+xi)*(1+eta) 0 1/4*(1-xi)*(1+eta)]);
        dN = [1/4*(-1+eta) 1/4*(1-eta) 1/4*(1+eta) 1/4*(-1-eta); 1/4*(-1+xi) 1/4*(-1-xi) 1/4*(1+xi) 1/4*(1-xi)];

        J = dN*r;

        he = hh/n;%sqrt((SX1-SX2)^2+(SY1-SY2)^2);
        Feq = Feq + c(i)*det(J)*T*(2/he) * transpose([N(2)*q*sin(pi/3) N(2)*q*cos(pi/3) N(4)*q*sin(pi/3) N(4)*q*cos(pi/3) N(6)*q*sin(pi/3) N(6)*q*cos(pi/3) N(8)*q*sin(pi/3) N(8)*q*cos(pi/3)]);
    end
end

Fe = Fet + Feg + Feq;

ne = [2*four(p,1)-1;2*four(p,1); 2*four(p,2)-1;2*four(p,2);2*four(p,3)-1;2*four(p,3);2*four(p,4)-1;2*four(p,4)];

for i=1:8
    for j = 1:8
    Kg(ne(i,1),ne(j,1)) = Kg(ne(i,1),ne(j,1))+Ke(i,j);
    end
end

    for i=1:8
    Fg(ne(i,1),1) = Fg(ne(i,1),1)+Fe(i,1);
    end
end

% Проверка на проекции сил
TESTsum = sum(Fg(1:2:end));
TESTeq = q*sin(pi/3)*4*T;

%Блок проверки глобальной матрицы Kg
%Симметрия
Pr = Kg - transpose(Kg);
max(max(abs(Pr)))
%Сумма по строкам
max(abs(sum(Kg,2)))
%Сумма по столбцам
max(abs(sum(Kg,1)))

for i=1:max(max(four))
    if NLIST(i,2)==0
       Kg(2*i-1,2*i-1)=Kg(2*i-1,2*i-1)*1e+10;
       Kg(2*i,2*i)=Kg(2*i,2*i)*1e+10;
    end
end

U = Kg\Fg;
M = 0.3;
Un = U/max(abs(U))*M;

XN = X  + Un(1:2:end,1);
YN = Y  + Un(2:2:end,1);

figure
hold on
grid on
ylim([min(YN)-delta max(YN)+delta])
xlim([min(X)-delta max(X)+delta])
for i=1:length(four)
    
for j=1:3
    plot([XN(four(i,j)),XN(four(i,j+1))],[YN(four(i,j)),YN(four(i,j+1))],'b.-')
    if j == 3
       plot([XN(four(i,j+1)),XN(four(i,1))],[YN(four(i,j+1)),YN(four(i,1))],'b.-') 
    end
end
end

sigmax = zeros(size(NLIST,1),1);
sigmay = zeros(size(NLIST,1),1);
tayxy = zeros(size(NLIST,1),1);
ch = zeros(size(NLIST,1),1);

for p =1:size(four,1)
    ne = [2*four(p,1)-1;2*four(p,1); 2*four(p,2)-1;2*four(p,2);2*four(p,3)-1;2*four(p,3);2*four(p,4)-1;2*four(p,4)];
    for j = 1:8
    Ue(j,1)=U(ne(j),1);
    end
    xi = 0;
    eta = 0;
    
    r = transpose([X(four(p,1)) X(four(p,2)) X(four(p,3)) X(four(p,4));Y(four(p,1)) Y(four(p,2)) Y(four(p,3)) Y(four(p,4))]);
    N = transpose([0 1/4*(1-xi)*(1-eta) 0 1/4*(1+xi)*(1-eta) 0 1/4*(1+xi)*(1+eta) 0 1/4*(1-xi)*(1+eta)]);
    dN = [1/4*(-1+eta) 1/4*(1-eta) 1/4*(1+eta) 1/4*(-1-eta); 1/4*(-1+xi) 1/4*(-1-xi) 1/4*(1+xi) 1/4*(1-xi)];


    J = dN*r;
    dN1 = inv(J)*[dN(1,1);dN(2,1)];
    dN2 = inv(J)*[dN(1,2);dN(2,2)];
    dN3 = inv(J)*[dN(1,3);dN(2,3)];
    dN4 = inv(J)*[dN(1,4);dN(2,4)];

    Be1 = [dN1(1) 0; 0 dN1(2);dN1(2) dN1(1)];
    Be2 = [dN2(1) 0; 0 dN2(2);dN2(2) dN2(1)];
    Be3 = [dN3(1) 0; 0 dN3(2);dN3(2) dN3(1)];
    Be4 = [dN4(1) 0; 0 dN4(2);dN4(2) dN4(1)];
    Be = [Be1 Be2 Be3 Be4];
    
    sigmae = De*Be*Ue;
    
   sigmax(four(p,1))=sigmax(four(p,1))+sigmae(1,1);
   sigmax(four(p,2))=sigmax(four(p,2))+sigmae(1,1);
   sigmax(four(p,3))=sigmax(four(p,3))+sigmae(1,1);
   sigmax(four(p,4))=sigmax(four(p,4))+sigmae(1,1);
   
   sigmay(four(p,1))=sigmay(four(p,1))+sigmae(2,1);
   sigmay(four(p,2))=sigmay(four(p,2))+sigmae(2,1);
   sigmay(four(p,3))=sigmay(four(p,3))+sigmae(2,1);
   sigmay(four(p,4))=sigmay(four(p,4))+sigmae(2,1);
   
   tayxy(four(p,1))=tayxy(four(p,1))+sigmae(3,1);
   tayxy(four(p,2))=tayxy(four(p,2))+sigmae(3,1);
   tayxy(four(p,3))=tayxy(four(p,3))+sigmae(3,1);
   tayxy(four(p,4))=tayxy(four(p,4))+sigmae(3,1);
   
   ch(four(p,1))=ch(four(p,1))+1;
   ch(four(p,2))=ch(four(p,2))+1;
   ch(four(p,3))=ch(four(p,3))+1;
   ch(four(p,4))=ch(four(p,4))+1;
end

sigmax=sigmax./ch;
sigmay=sigmay./ch;
tayxy= tayxy./ch;

figure
trisurf(four,X,Y,sigmax)
title('\sigma_x')
view(0,90)
shading interp
colorbar;

figure
trisurf(four,X,Y,sigmay)
title('\sigma_y')
view(0,90)
shading interp
colorbar;

figure
trisurf(four,X,Y,tayxy)
title('\tau_x_y')
view(0,90)
shading interp
colorbar;

P2 = 0;
Nr0 = 0;

for i=1:max(max(four))
    if abs(NLIST(i,2)+sqrt(3)*NLIST(i,1) - tan(pi/3)*4) < delta && NLIST(i,2)>=0.2 && NLIST(i,2)<3.2
        P2 = P2 + (sigmax(i)*sin(pi/3)+tayxy(i)*cos(pi/3)-q*sin(pi/3))^2+(tayxy(i)*sin(pi/3)+sigmay(i)*cos(pi/3)-q*cos(pi/3))^2;
        Nr0 = Nr0+1;
    end
end

P2R = (1/(2*Nr0))*P2;
Po = sqrt(P2R)/abs(q)*100;

end
