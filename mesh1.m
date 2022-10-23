function [X, Y, four1] = mesh1(X,Y,A,n)
    alfa1=120;
    X1 = (X-A)*cosd(alfa1)+Y*sind(alfa1);
    Y1 = -(X-A)*sind(alfa1)+Y*cosd(alfa1);
    
    X2 = (X1-A)*cosd(alfa1)+Y1*sind(alfa1);
    Y2 = -(X1-A)*sind(alfa1)+Y1*cosd(alfa1);
    %удаление пересекаемых узлов
    for i = 1:length(X1)-n-1
        for j = 1:length(X)
            if abs(X(j)-X1(i))<0.01 && abs(Y(j)-Y1(i))<0.01
                X1(i) = [];
                Y1(i) = []; 
            end
        end
    end
    X1(end) = [];
    Y1(end) = [];
    
    for i = 1:length(X2)-n-1
        for j = 1:length(X)
            if abs(X(j)-X2(i))<0.01 && abs(Y(j)-Y2(i))<0.01
                X2(i) = [];
                Y2(i) = []; 
            end
        end
    end
    X2(end) = [];
    Y2(end) = [];
    
    for i = 1:length(X2)-n
       if X2(i) - A/2 <0.0000001
           Y2(i) = [];
           X2(i) = [];
       end
    end
    
    % figure
    % plot(X,Y,'.')
    % hold on
    % plot(X1,Y1,'.')
    % hold on 
    % xlim([0 A])
    % ylim([0 A*cosd(30)])
    % plot(X,Y,'.')
    % grid on 
    k1 = length(X);
    k2 = length(X1);
    k3 = length(X2); %количество узлов
    X = [X X1 X2]';
    Y = [Y Y1 Y2]';
    figure
    plot(X,Y,'*')
    hold on
    grid on
    
    
    four1 = [];
    i = 1;
    k = 1;
    %получение матрицы нумерации узлов для нижней трапеции
    while k < k1-n
        if  mod(k,n+1) ~= 0
            four1(i,1)=k;
            four1(i,2)=k+1;
            four1(i,3)=k+n+2;
            four1(i,4)=k+n+1;
            i = i + 1;
        else
            i = i;
        end
            k = k + 1;   
    end
    
    %левая часть
    if1 = 1;
    iter = 1; %слой матрицы four
    j = length(four1)+1;
    k=k1+1;
    for i=k1:k1+k2
        if if1 ~= n
            four1(j,1)=k;
            four1(j,2)=k+1;
            four1(j,3)=k+n+1;
            four1(j,4)=k+n;
            if1 = if1+1;
        else
            four1(j,1)=k;
            four1(j,2)=(n+1)*(iter-1)+1;
            four1(j,3)=(n+1)*iter+1;
            four1(j,4)=k+n;
            iter = iter + 1;
            if1 = 1;
        end
        k = k+1;
        j=j+1;
    end
    
    %удаление лишних узлов
    for i=1:n+1
       four1(end,:)=[];
    end
    
    %правая часть
    if2 = 1;
    if3 = 1;
    iter2 = 1;
    iter1 = 1;
    j = length(four1)+1;
    k = k1+k2+1;
    for i=k1+k2:k1+k2+k3
        if if2 == 1
            four1(j,1)=(n+1)*iter1;
            four1(j,2)=k;
            four1(j,3)=k+n-1;
            four1(j,4)=(n+1)*(iter1+1);
            k = k+1;
            iter1 = iter1 + 1;
            j=j+1;
            if2 = 0;
            if3 = if3 + 1;
        end
        
        if if3 == n
            four1(j,1)=k-1;
            four1(j,2)=k1+n*(iter2-1)+1;
            four1(j,3)=k1+n*iter2+1;
            four1(j,4)=k+n-2;
            iter2 = iter2 + 1;
            if2 = 1;
            if3 = 1;
            j = j + 1;
            k = k;
        else
            four1(j,1)=k-1;         
            four1(j,2)=k;
            four1(j,3)=k+n-1;
            four1(j,4)=k+n-2;
            if3 = if3 + 1;
            j = j + 1;
            k = k+1;
        end
    end
    
    % %удаление лишних узлов
    for i=1:n+2
       four1(end,:)=[];
    end
    
    %Отображение сетки
    figure
    hold on
    for i=1:length(four1) 
        for j=1:3
            plot([X(four1(i,j)),X(four1(i,j+1))],[Y(four1(i,j)),Y(four1(i,j+1))],'b.-')
            if j == 3
               plot([X(four1(i,j+1)),X(four1(i,1))],[Y(four1(i,j+1)),Y(four1(i,1))],'b.-') 
            end
        end
    end
    grid on

end