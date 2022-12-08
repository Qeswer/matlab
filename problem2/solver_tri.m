function [tri, Po, sum_stroca, sum_stolb] = solver_tri(tri, X, Y,Fg, Kg)

%---------------------константы--------------------%
%Задание параметров Ламе через технические константы
E=2e+11;
w=0.33;
Ro = 7800;
alpha = 11.7e-6;
t = 0.001; %толщина
A = 4;
%Задание матрицы упругих констант
h=(E*w)/((1+w)*(1-2*w));
m=E/(2*(1+w));
De = [h+2*m h 0;h h+2*m 0;0 0 m];
%---------------------------------------------------%

g = - 0;%input("ускроние свободного падения: ");
dT = 0;%input("температурное изменение: ");
q = 100;%input("распределенная нагрузка: ");

e0 = [alpha*dT; alpha*dT; 0];
% Kg = zeros(2*max(max(tri)));
% Fg = zeros(2*max(max(tri)),1);
p = 1;
% for p =1:size(tri,1)
% figure
% while p<size(tri,1)
%     
% %    S = [1 X(tri(p,1)) Y(tri(p,1)); 1 X(tri(p,2)) ...
% %        Y(tri(p,2));1 X(tri(p,3)) Y(tri(p,3))]; 
% %    Ae =1/2*det(S);
%     x1 = X(tri(p,1));
%     x2 = X(tri(p,2));
%     x3 = X(tri(p,3));
%     y1 = Y(tri(p,1));
%     y2 = Y(tri(p,2));
%     y3 = Y(tri(p,3));
%     Ae = 1/2*abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
%    
%    if abs(Ae)<1e-5
%     tri(p,:)=[];
%     continue
%    end   
%    
%    if Ae<0
%        S = [1 X(tri(p,2)) Y(tri(p,2));1 X(tri(p,1)) Y(tri(p,1));1 X(tri(p,3)) Y(tri(p,3))];
%        Ae = det(S);
%    end
%    
%    AE(p) = Ae;
% %    ai = X(tri(p,2))*Y(tri(p,3)) - X(tri(p,3))*Y(tri(p,2));
% %    aj = X(tri(p,3))*Y(tri(p,1)) - X(tri(p,1))*Y(tri(p,3));
% %    am = X(tri(p,1))*Y(tri(p,2)) - X(tri(p,2))*Y(tri(p,1));
% %    
% %    bi = Y(tri(p,2)) - Y(tri(p,3));
% %    bj = Y(tri(p,3)) - Y(tri(p,1));
% %    bm = Y(tri(p,1)) - Y(tri(p,2));
% %    
% %    ci = X(tri(p,3)) - Y(tri(p,2));
% %    cj = X(tri(p,1)) - Y(tri(p,3));
% %    cm = X(tri(p,2)) - Y(tri(p,1));
% %    
% %    Be1 = [bi 0; 0 ci; ci bi];
% %    Be2 = [bj 0; 0 cj; cj bj];
% %    Be3 = [bm 0; 0 cm; cm bm];
%    Be1 = [Y(tri(p,2))-Y(tri(p,3)) 0; 0 X(tri(p,3))-X(tri(p,2)); X(tri(p,3))-X(tri(p,2)) Y(tri(p,2))-Y(tri(p,3))];
%    Be2 = [Y(tri(p,3))-Y(tri(p,1)) 0; 0 X(tri(p,1))-X(tri(p,3)); X(tri(p,1))-X(tri(p,3)) Y(tri(p,3))-Y(tri(p,1))];
%    Be3 = [Y(tri(p,1))-Y(tri(p,2)) 0; 0 X(tri(p,2))-X(tri(p,1)); X(tri(p,2))-X(tri(p,1)) Y(tri(p,1))-Y(tri(p,2))];
%    Be = [Be1 Be2 Be3]/(2*Ae);
%    BeT = transpose(Be);
% 
%    Ke = BeT*De*Be*t*Ae;
%    
% %    e0 = [alpha*dT; alpha*dT; 0];
% %    fet = BeT*De*e0*t*Ae;
% %    feg = [0; Ro*g*t*Ae/3; 0; Ro*g*t*Ae/3; 0; Ro*g*t*Ae/3];
%     delta = Ae/1e+6;
% 
% %     Spov = sqrt((X(tri(p,2))-X(tri(p,1)))^2+(Y(tri(p,2))-Y(tri(p,1)))^2);
% % %      Spov = A/n;
% %      feq = zeros(6,1);
% % 
% %      if and(abs(Y(tri(p,1))+sqrt(3)*X(tri(p,1)) - A*tan(pi/3)) < delta, abs(Y(tri(p,2))+sqrt(3)*X(tri(p,2)) - A*tan(pi/3)) < delta)
% % %          plot(X(tri(p,1)), Y(tri(p,1)), "*",X(tri(p,2)), Y(tri(p,2)), "o")
% %      feq = -q*[(sin(pi/3))*Spov/2; (cos(pi/3))*Spov/2;(sin(pi/3))*Spov/2;(cos(pi/3))*Spov/2;0;0];
% %      end 
% %      
% %      if and(abs(Y(tri(p,1))+sqrt(3)*X(tri(p,1)) - A*tan(pi/3)) < delta, abs(Y(tri(p,3))+sqrt(3)*X(tri(p,3)) - A*tan(pi/3)) < delta)
% % %         plot(X(tri(p,1)), Y(tri(p,1)), "*",X(tri(p,1)), Y(tri(p,1)), "o")
% %      feq = -q*[(sin(pi/3))*Spov/2; (cos(pi/3))*Spov/2;0;0;(sin(pi/3))*Spov/2;(cos(pi/3))*Spov/2];
% %      end
% %      
% %      if and(abs(Y(tri(p,2))+sqrt(3)*X(tri(p,2)) - A*tan(pi/3)) < delta, abs(Y(tri(p,3))+sqrt(3)*X(tri(p,3)) - A*tan(pi/3)) < delta)
% % %          plot(X(tri(p,2)), Y(tri(p,2)), "+")
% %      feq = -q*[0;0;(sin(pi/3))*Spov/2; (cos(pi/3))*Spov/2;(sin(pi/3))*Spov/2;(cos(pi/3))*Spov/2];
% %      end
%      
% %      fe = fet + feg + feq;
%      
%      
%      ne = [2*tri(p,1)-1;2*tri(p,1); 2*tri(p,2)-1;2*tri(p,2);2*tri(p,3)-1;2*tri(p,3)];
%      
%      for i=1:6
%          for j = 1:6
%             Kg(ne(i,1),ne(j,1)) = Kg(ne(i,1),ne(j,1))+Ke(i,j);
%          end
%      end
%      
% %      for i=1:6
% %         Fg(ne(i,1),1) = Fg(ne(i,1),1)+fe(i,1);
% %      end
%      p = p+1;
% end

%Блок проверки глобальной матрицы Kg
%Симметрия
Pr = Kg - transpose(Kg);
max(max(abs(Pr)));
%Сумма по строкам
sum_stroca = max(abs(sum(Kg,2)));
%Сумма по стролбцам
sum_stolb = max(abs(sum(Kg,1)));

for i=1:max(max(tri))
    if Y(i)==0
       Kg(2*i-1,2*i-1)=Kg(2*i-1,2*i-1)*1e+10;
       Kg(2*i,2*i)=Kg(2*i,2*i)*1e+10;
    end
end

U = Kg\Fg;
M = 0.3;
Un = U/max(abs(U))*M;

Xn = X + Un(1:2:end,1);
Yn = Y + Un(2:2:end,1);
figure
triplot(tri, Xn, Yn)
grid on

sigmax = zeros(length(X),1);
sigmay = zeros(length(X),1);
tayxy = zeros(length(X),1);
ch = zeros(length(X),1);

for i=1:size(tri,1)
    ne = [2*tri(i,1)-1;2*tri(i,1); 2*tri(i,2)-1;2*tri(i,2);2*tri(i,3)-1;2*tri(i,3)];
    for j = 1:6
        Ue(j,1)=U(ne(j),1);
    end
    
%    S = [X(tri(i,1)) Y(tri(i,1)) 1;X(tri(i,2)) Y(tri(i,2)) 1; X(tri(i,3)) Y(tri(i,3)) 1]; 
%    Ae = det(S);
    x1 = X(tri(i,1));
    x2 = X(tri(i,2));
    x3 = X(tri(i,3));
    y1 = Y(tri(i,1));
    y2 = Y(tri(i,2));
    y3 = Y(tri(i,3));
    Ae = 1/2*abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
   if Ae<0
       S = [X(tri(i,2)) Y(tri(i,2)) 1;X(tri(i,1)) Y(tri(i,1)) 1; X(tri(i,3)) Y(tri(i,3)) 1];
       Ae = det(S);
   end
%    ai = X(tri(i,2))*Y(tri(i,3)) - X(tri(i,3))*Y(tri(i,2));
%    aj = X(tri(i,1))*Y(tri(i,3)) - X(tri(i,3))*Y(tri(i,1));
%    am = X(tri(i,2))*Y(tri(i,1)) - X(tri(i,1))*Y(tri(i,2));
%    
%    bi = Y(tri(i,2)) - Y(tri(i,3));
%    bj = Y(tri(i,1)) - Y(tri(i,3));
%    bm = Y(tri(i,2)) - Y(tri(i,1));
%    
%    ci = X(tri(i,3)) - Y(tri(i,2));
%    cj = X(tri(i,3)) - Y(tri(i,1));
%    cm = X(tri(i,1)) - Y(tri(i,2));
%    
%    Be1 = [bi 0; 0 ci; ci bi];
%    Be2 = [bj 0; 0 cj; cj bj];
%    Be3 = [bm 0; 0 cm; cm bm];
   Be1 = [Y(tri(i,2))-Y(tri(i,3)) 0; 0 X(tri(i,3))-X(tri(i,2)); X(tri(i,3))-X(tri(i,2)) Y(tri(i,2))-Y(tri(i,3))];
   Be2 = [Y(tri(i,3))-Y(tri(i,1)) 0; 0 X(tri(i,1))-X(tri(i,3)); X(tri(i,1))-X(tri(i,3)) Y(tri(i,3))-Y(tri(i,1))];
   Be3 = [Y(tri(i,1))-Y(tri(i,2)) 0; 0 X(tri(i,2))-X(tri(i,1)); X(tri(i,2))-X(tri(i,1)) Y(tri(i,1))-Y(tri(i,2))];
   Be = [Be1 Be2 Be3]/(2*Ae);
   
   sigmae = De*Be*Ue;
   
   sigmax(tri(i,1))=sigmax(tri(i,1))+sigmae(1,1);
   sigmax(tri(i,2))=sigmax(tri(i,2))+sigmae(1,1);
   sigmax(tri(i,3))=sigmax(tri(i,3))+sigmae(1,1);
   
   sigmay(tri(i,1))=sigmay(tri(i,1))+sigmae(2,1);
   sigmay(tri(i,2))=sigmay(tri(i,2))+sigmae(2,1);
   sigmay(tri(i,3))=sigmay(tri(i,3))+sigmae(2,1);
   
   tayxy(tri(i,1))=tayxy(tri(i,1))+sigmae(3,1);
   tayxy(tri(i,2))=tayxy(tri(i,2))+sigmae(3,1);
   tayxy(tri(i,3))=tayxy(tri(i,3))+sigmae(3,1);
   
   ch(tri(i,1))=ch(tri(i,1))+1;
   ch(tri(i,2))=ch(tri(i,2))+1;
   ch(tri(i,3))=ch(tri(i,3))+1;
end

sigmax=sigmax./ch;
sigmay=sigmay./ch;
tayxy= tayxy./ch;

figure
trisurf(tri,X,Y,sigmax)
title('\sigma_x tri')
view(0,90)
shading interp
colorbar;

figure
trisurf(tri,X,Y,sigmay)
title('\sigma_y tri')
view(0,90)
shading interp
colorbar;

figure
trisurf(tri,X,Y,tayxy)
title('tay_xy tri')
view(0,90)
shading interp
colorbar;

P2 = 0;
Nr0 = 0;
q = -100;
for i=1:max(max(tri))
    if abs(Y(i)+sqrt(3)*X(i) - tan(pi/3)*4) < 0.001 && Y(i)>=0.2 && Y(i)<3.2
        P2 = P2 + (sigmax(i)*sin(pi/3)+tayxy(i)*cos(pi/3)-q*sin(pi/3))^2+(tayxy(i)*sin(pi/3)+sigmay(i)*cos(pi/3)-q*cos(pi/3))^2;
        Nr0 = Nr0+1;
    end
end

P2R = (1/(2*Nr0))*P2;
Po = sqrt(P2R)/abs(q)*100-100;
sU3 = size(U);
   
