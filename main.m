clear all
close all
clc

%----------------Геометрические характеристики фигуры-------------%

a=1;        %сторона внутреннего треугольника
A=1.5*a+a+1.5*a; %сторона внешнего треугольника
n = input("кол-во делений: "); %деление
h = (A-a)/2*tand(30); %высота до выреза
X = [];
Y = [];
hn = h/n; %высота элемента
len = A;
s = (A-a)/2/n; %сдвиг
for i=0:n
    dlen = len/n; %длина стороны элемента
    for j = 1:n+1
        X(length(X)+1)=dlen*(j-1)+s*i;
        Y(length(Y)+1)=hn*i;
    end
    len = len-2*s;
end    

[X, Y, four1] = mesh1(X,Y,A,n);

% solve1(X, Y, four1, z, c);

solver(X, Y, four1)


