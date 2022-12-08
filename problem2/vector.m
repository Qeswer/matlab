function [tri] = vector(X,Y)

tri = delaunay(X,Y);

    xA = 1.5;
    yA = 0.866025;
    xB = 2;
    yB = 1.73205;
    xC = 2.5;
    yC = yA;
    i=1;
    n = 1;
    while i<size(tri,1)
        if (xA-X(tri(i,1)))*(yB-yA)-(xB-xA)*(yA-Y(tri(i,1)))<=0.01 &&...
                (xB-X(tri(i,1)))*(yC-yB)-(xC-xB)*(yB-Y(tri(i,1)))<=0.01 &&...
                (xC-X(tri(i,1)))*(yA-yC)-(xA-xC)*(yC-Y(tri(i,1)))<=0.01
            if (xA-X(tri(i,2)))*(yB-yA)-(xB-xA)*(yA-Y(tri(i,2)))<=0.01 &&...
                (xB-X(tri(i,2)))*(yC-yB)-(xC-xB)*(yB-Y(tri(i,2)))<=0.01 &&...
                (xC-X(tri(i,2)))*(yA-yC)-(xA-xC)*(yC-Y(tri(i,2)))<=0.01
                    if (xA-X(tri(i,3)))*(yB-yA)-(xB-xA)*(yA-Y(tri(i,3)))<=0.01 &&...
                        (xB-X(tri(i,3)))*(yC-yB)-(xC-xB)*(yB-Y(tri(i,3)))<=0.01 &&...
                        (xC-X(tri(i,3)))*(yA-yC)-(xA-xC)*(yC-Y(tri(i,3)))<=0.01
                            n = n + 1;
                            tri(i,:)=[];
                            i = i - 1;
                    end
            end
        end
        i = i + 1;
    end


% for p = 1:size(tri,1)
%     if abs(Y(tri(p,1))+sqrt(3)*X(tri(p,1)) - 4*tan(pi/3)) < 0.001 &&...
%             abs(Y(tri(p,2))+sqrt(3)*X(tri(p,2)) - 4*tan(pi/3)) < 0.001
% %         plot(X(tri(p,1)),Y(tri(p,1)),"o")
% %         plot(X(tri(p,2)),Y(tri(p,2)),"*")
%     end
%     if abs(Y(tri(p,1))+sqrt(3)*X(tri(p,1)) - 4*tan(pi/3)) < 0.001 &&...
%             abs(Y(tri(p,3))+sqrt(3)*X(tri(p,3)) - 4*tan(pi/3)) < 0.001
% %         plot(X(tri(p,1)),Y(tri(p,1)),"+")
% %         plot(X(tri(p,3)),Y(tri(p,3)),"-*")
%     end
%     if abs(Y(tri(p,2))+sqrt(3)*X(tri(p,2)) - 4*tan(pi/3)) < 0.001 &&...
%             abs(Y(tri(p,3))+sqrt(3)*X(tri(p,3)) - 4*tan(pi/3)) < 0.001
%         plot(X(tri(p,2)),Y(tri(p,2)),"+")
%         plot(X(tri(p,3)),Y(tri(p,3)),"o")
%     end
p=1;
while p<size(tri,1)
    
%    S = [1 X(tri(p,1)) Y(tri(p,1)); 1 X(tri(p,2)) ...
%        Y(tri(p,2));1 X(tri(p,3)) Y(tri(p,3))]; 
%    Ae =1/2*det(S);
    x1 = X(tri(p,1));
    x2 = X(tri(p,2));
    x3 = X(tri(p,3));
    y1 = Y(tri(p,1));
    y2 = Y(tri(p,2));
    y3 = Y(tri(p,3));
    Ae = 1/2*abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
   
   if abs(Ae)<1e-5
    tri(p,:)=[];
    continue
   end   
   p = p+1;
end
figure
hold on
triplot(tri,X,Y)
