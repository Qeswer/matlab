function solver(X, Y, four)


%---------------------константы--------------------%
%Задание параметров Ламе через технические константы
E=2e+11;
w=0.33;
Ro = 7800;
alpha = 11.7e-6;
T = 0.01; %толщина
%Задание матрицы упругих констант
h=(E*w)/((1+w)*(1-2*w));
m=E/(2*(1+w));
De = [h+2*m h 0;h h+2*m 0;0 0 m];
%---------------------------------------------------%

g = - input("ускроние свободного падения: ");
dT = input("температурное изменение: ");
q = - input("распределенная нагрузка: ");

e0 = [alpha*dT; alpha*dT; 0];
Kg = zeros(2*max(max(four)));
Fg = zeros(2*max(max(four)),1);

XY = [X Y];

for p =1:size(four,1)
    Ke = zeros(8);
    Fet = zeros(8,1);
    Feg = zeros(8,1);
    Feq = zeros(8,1);
    X1 = X(four(p,1));
    X2 = X(four(p,2));
    X3 = X(four(p,3));
    X4 = X(four(p,4));
    Y1 = Y(four(p,1));
    Y2 = Y(four(p,2));
    Y3 = Y(four(p,3));
    Y4 = Y(four(p,4));
    S = [X1-X3 X4-X2;Y1-Y3 Y4-Y2];
    Ae = abs((1/2)*det(S));

    c = [0.555555556 0.888888889 0.555555556];
    z = [-0.774596669 0.000 0.774596669];

    r = [X1 X2 X3 X4;Y1 Y2 Y3 Y4]';
    for i=1:3
        for j=1:3
            xi = z(i);
            eta = z(j);

            if and(X1<X2, X4<X3)
                N = [0 1/4*(1-xi)*(1-eta) 0 1/4*(1+xi)*(1-eta) 0 1/4*(1+xi)*(1+eta) 0 1/4*(1-xi)*(1+eta)]';
                dN = [1/4*(-1+eta) 1/4*(1-eta) 1/4*(1+eta) 1/4*(-1-eta); 1/4*(-1+xi) 1/4*(-1-xi) 1/4*(1+xi) 1/4*(1-xi)];
            end

            if and(Y1<Y2, Y4<Y3)
                N = [0 1/4*(1-xi)*(1+eta) 0 1/4*(1-xi)*(1-eta) 0 1/4*(1+xi)*(1-eta) 0 1/4*(1+xi)*(1+eta)]';
                dN = [1/4*(-1-eta) 1/4*(-1+eta) 1/4*(1-eta) 1/4*(1+eta); 1/4*(1-xi) 1/4*(-1+xi) 1/4*(-1-xi) 1/4*(1+xi)];
            end

            if and(X1>X2, X4>X3)
                N = [0 1/4*(1+xi)*(1+eta) 0 1/4*(1-xi)*(1+eta) 0 1/4*(1-xi)*(1-eta) 0 1/4*(1+xi)*(1-eta)]';
                dN = [1/4*(1+eta) 1/4*(-1-eta) 1/4*(-1+eta) 1/4*(1-eta); 1/4*(1+xi) 1/4*(1-xi) 1/4*(-1+xi) 1/4*(-1-xi)];
            end

            if and(Y1>Y2, Y4>Y3)
                N = [0 1/4*(1+xi)*(1-eta) 0 1/4*(1+xi)*(1+eta) 0 1/4*(1-xi)*(1+eta) 0 1/4*(1-xi)*(1-eta)]';
                dN = [1/4*(1-eta) 1/4*(1+eta) 1/4*(-1-eta) 1/4*(-1+eta); 1/4*(-1-xi) 1/4*(1+xi) 1/4*(1-xi) 1/4*(-1+xi)];
            end


            J = dN*r;
            dN1 = J\[dN(1,1);dN(2,1)];
            dN2 = J\[dN(1,2);dN(2,2)];
            dN3 = J\[dN(1,3);dN(2,3)];
            dN4 = J\[dN(1,4);dN(2,4)];

            Be1 = [dN1(1) 0; 0 dN1(2);dN1(2) dN1(1)];
            Be2 = [dN2(1) 0; 0 dN2(2);dN2(2) dN2(1)];
            Be3 = [dN3(1) 0; 0 dN3(2);dN3(2) dN3(1)];
            Be4 = [dN4(1) 0; 0 dN4(2);dN4(2) dN4(1)];
            Be = [Be1 Be2 Be3 Be4];
            BeT = Be';
            Ke = Ke + c(i)*c(j)*BeT*De*Be*det(J)*T;

            Fet = Fet + c(i)*c(j)*BeT*De*e0*det(J)*T;
            Feg = Feg + c(i)*c(j)*Ro*g*det(J)*T*N;
        end
    end

    delta = Ae/1e+6;

     if and((Y1+sqrt(3)*X1 - tan(pi/3)*4) < delta, abs(Y2+sqrt(3)*X2 - tan(pi/3)*4) < delta)
        eta = -1;
        dN = [1/4*(1+eta) 1/4*(-1-eta) 1/4*(-1+eta) 1/4*(1-eta); 1/4*(1+xi) 1/4*(1-xi) 1/4*(-1+xi) 1/4*(-1-xi)];
        for i=1:3
            xi = z(i);
            J = dN*r;
            h = sqrt((X1-X2)^2+(Y1-Y2)^2);
            Feq = Feq + c(i)*det(J)*T*h/2*[N(2)*q*sin(pi/3) N(2)*q*cos(pi/3) N(4)*q*sin(pi/3) N(4)*q*cos(pi/3) N(6)*q*sin(pi/3) N(6)*q*cos(pi/3) N(8)*q*sin(pi/3) N(8)*q*cos(pi/3)]';
        end
     end
 
    if and(abs(Y2+sqrt(3)*X2 - tan(pi/3)*4) < delta, abs(Y3+sqrt(3)*X3 - tan(pi/3)*4) < delta)
        xi = 1;
        dN = [1/4*(-1-eta) 1/4*(-1+eta) 1/4*(1-eta) 1/4*(1+eta); 1/4*(1-xi) 1/4*(-1+xi) 1/4*(-1-xi) 1/4*(1+xi)];
        for i=1:3
            eta = z(i);
            J = dN*r;
            h = sqrt((X2-X3)^2+(Y2-Y3)^2);
            Feq = Feq + c(i)*det(J)*T*h/2 * [N(2)*q*sin(pi/3) N(2)*q*cos(pi/3) N(4)*q*sin(pi/3) N(4)*q*cos(pi/3) N(6)*q*sin(pi/3) N(6)*q*cos(pi/3) N(8)*q*sin(pi/3) N(8)*q*cos(pi/3)]';
        end
    end
 
    if and(abs(Y3+sqrt(3)*X3 - tan(pi/3)*4) < delta, abs(Y4+sqrt(3)*X4 - tan(pi/3)*4) < delta)
        eta = 1;
        dN = [1/4*(-1+eta) 1/4*(1-eta) 1/4*(1+eta) 1/4*(-1-eta); 1/4*(-1+xi) 1/4*(-1-xi) 1/4*(1+xi) 1/4*(1-xi)];
        for i=1:3
            xi = z(i);
            J = dN*r;
            h = sqrt((X3-X4)^2+(Y3-Y4)^2);
            Feq = Feq + c(i)*det(J)*T*h/2 * [N(2)*q*sin(pi/3) N(2)*q*cos(pi/3) N(4)*q*sin(pi/3) N(4)*q*cos(pi/3) N(6)*q*sin(pi/3) N(6)*q*cos(pi/3) N(8)*q*sin(pi/3) N(8)*q*cos(pi/3)];
        end
    end
 
    if and(abs(Y4+sqrt(3)*X4 - tan(pi/3)*4) < delta, abs(Y1+sqrt(3)*X1 - tan(pi/3)*4) < delta)
        xi = -1;
        dN = [1/4*(1-eta) 1/4*(1+eta) 1/4*(-1-eta) 1/4*(-1+eta); 1/4*(-1-xi) 1/4*(1+xi) 1/4*(1-xi) 1/4*(-1+xi)];
        for i=1:3
            eta = z(i);
            J = dN*r;
            h = sqrt((X1-X4)^2+(Y4-Y1)^2);
            Feq = Feq + c(i)*det(J)*T*h/2 * [N(2)*q*sin(pi/3) N(2)*q*cos(pi/3) N(4)*q*sin(pi/3) N(4)*q*cos(pi/3) N(6)*q*sin(pi/3) N(6)*q*cos(pi/3) N(8)*q*sin(pi/3) N(8)*q*cos(pi/3)]';
        end
    end

    Fe = Fet + Feg + Feq;

    ne = [2*four(p,1)-1;2*four(p,1); 2*four(p,2)-1;2*four(p,2);2*four(p,3)-1;2*four(p,3);2*four(p,4)-1;2*four(p,4)];

    for i = 1:8
        for j = 1:8
            Kg(ne(i,1),ne(j,1)) = Kg(ne(i,1),ne(j,1))+Ke(i,j);
        end
    end

    for i=1:8
        Fg(ne(i,1),1) = Fg(ne(i,1),1)+Fe(i,1);
    end
end

%Блок проверки глобальной матрицы Kg
%Симметрия
Pr = Kg - Kg';
max(max(abs(Pr)))
%Сумма по строкам
max(abs(sum(Kg,2)))
%Сумма по столбцам
max(abs(sum(Kg,1)))

%ГУ
for i=1:max(max(four))
    if Y(i) == 0
       Kg(2*i-1,2*i-1)=Kg(2*i-1,2*i-1)*1e+10;
       Kg(2*i,2*i)=Kg(2*i,2*i)*1e+10;
    end
end

U = Kg\Fg;
M = 0.3;
Un = U/max(abs(U))*M;

XYPX = zeros(length(XY),1);
XYPY = zeros(length(XY),1);

for i=1:length(XY)
    XYPX(i) = XY(i,1);
    XYPY(i) = XY(i,2);
end

XYPX = XYPX + Un(1:2:end,1);
XYPY = XYPY + Un(2:2:end,1);
XYNEW = [XYPX XYPY];

figure
hold on
grid on
for i=1:length(four)
    for j=1:3
        plot([XYNEW(four(i,j),1),XYNEW(four(i,j+1),1)],[XYNEW(four(i,j),2),XYNEW(four(i,j+1),2)],'b.-')
        if j == 3
           plot([XYNEW(four(i,j+1),1),XYNEW(four(i,1),1)],[XYNEW(four(i,j+1),2),XYNEW(four(i,1),2)],'b.-') 
        end
    end
end

sigmax = zeros(size(XY,1),1);
sigmay = zeros(size(XY,1),1);
tayxy = zeros(size(XY,1),1);
ch = zeros(size(XY,1),1);

for p =1:size(four,1)
    ne = [2*four(p,1)-1;2*four(p,1); 2*four(p,2)-1;2*four(p,2);2*four(p,3)-1;2*four(p,3);2*four(p,4)-1;2*four(p,4)];
    for j = 1:8
        Ue(j,1)=U(ne(j),1);
    end
    xi = 0;
    eta = 0;
    
    r = [X1 X2 X3 X4;Y1 Y2 Y3 Y4]';
    if and(X1<X2, X4<X3)
        N = [0 1/4*(1-xi)*(1-eta) 0 1/4*(1+xi)*(1-eta) 0 1/4*(1+xi)*(1+eta) 0 1/4*(1-xi)*(1+eta)]';
        dN = [1/4*(-1+eta) 1/4*(1-eta) 1/4*(1+eta) 1/4*(-1-eta); 1/4*(-1+xi) 1/4*(-1-xi) 1/4*(1+xi) 1/4*(1-xi)];
    end

    if and(Y1<Y2, Y4<Y3)
        N = [0 1/4*(1-xi)*(1+eta) 0 1/4*(1-xi)*(1-eta) 0 1/4*(1+xi)*(1-eta) 0 1/4*(1+xi)*(1+eta)]';
        dN = [1/4*(-1-eta) 1/4*(-1+eta) 1/4*(1-eta) 1/4*(1+eta); 1/4*(1-xi) 1/4*(-1+xi) 1/4*(-1-xi) 1/4*(1+xi)];
    end

    if and(X1>X2, X4>X3)
        N = [0 1/4*(1+xi)*(1+eta) 0 1/4*(1-xi)*(1+eta) 0 1/4*(1-xi)*(1-eta) 0 1/4*(1+xi)*(1-eta)]';
        dN = [1/4*(1+eta) 1/4*(-1-eta) 1/4*(-1+eta) 1/4*(1-eta); 1/4*(1+xi) 1/4*(1-xi) 1/4*(-1+xi) 1/4*(-1-xi)];
    end

    if and(Y1>Y2, Y4>Y3)
        N = [0 1/4*(1+xi)*(1-eta) 0 1/4*(1+xi)*(1+eta) 0 1/4*(1-xi)*(1+eta) 0 1/4*(1-xi)*(1-eta)]';
        dN = [1/4*(1-eta) 1/4*(1+eta) 1/4*(-1-eta) 1/4*(-1+eta); 1/4*(-1-xi) 1/4*(1+xi) 1/4*(1-xi) 1/4*(-1+xi)];
    end


    J = dN*r;
    dN1 = J\[dN(1,1);dN(2,1)];
    dN2 = J\[dN(1,2);dN(2,2)];
    dN3 = J\[dN(1,3);dN(2,3)];
    dN4 = J\[dN(1,4);dN(2,4)];

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
trisurf(four,XYPX,XYPY,sigmax)
view(0,90)
ylim([min(Y)-delta max(XYPY)+delta])
xlim([min(X)-delta max(X)+delta])
shading interp
colorbar;

figure
trisurf(four,XYPX,XYPY,sigmay)
view(0,90)
ylim([min(Y)-delta max(XYPY)+delta])
xlim([min(X)-delta max(X)+delta])
shading interp
colorbar;

figure
trisurf(four,XYPX,XYPY,tayxy)
view(0,90)
ylim([min(Y)-delta max(XYPY)+delta])
xlim([min(X)-delta max(X)+delta])
shading interp
colorbar;


P2 = 0;
Nr0 = 0;

for i=1:max(max(four))
    if abs(XY(i,2)+sqrt(3)*XY(i,1) - tan(pi/3)*4) < delta && XY(i,2)>=0.2 && XY(i,2)<3.2
        P2 = P2 + (sigmax(i)*sin(pi/3)+tayxy(i)*cos(pi/3)+q*sin(pi/3))^2+(tayxy(i)*sin(pi/3)+sigmay(i)*cos(pi/3)+q*cos(pi/3))^2;
        Nr0 = Nr0+1;
    end
end
    

P2R = (1/(2*Nr0))*P2;
Po = sqrt(P2R)/abs(q)*100

end