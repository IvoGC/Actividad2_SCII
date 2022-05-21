clc; clear ; close all
%{
---------------------------------------------------------------------------
Implementar un sistema en variables de estado que controle el ángulo del 
motor, para consignas de pi/2 y –pi/2 cambiando cada 300 mili segundos y 
que el TL de 1,15 10-3 aparece sólo para pi/2, para –pi/2 es nulo. 
Hallar el valor de integración Euler adecuado. El objetivo es mejorar 
la dinámica del controlador que muestra la Fig. 1. 
---------------------------------------------------------------------------
%}
%DEFINO PARAMETROS
LAA = 366*10^-6;
J = 5*10^-9;
RA = 55.6;
Bm = 0;
Ki = 6.49*10^-3;
Km = 6.53*10^-3;

%DEFINO MATRICES
%X=[ia ; tita ; w];
A=[-RA/LAA 0 -Km/LAA  ; 0 0 1 ; Ki/J 0 -Bm/J];
B=[1/LAA; 0; 0];
C=[0 1 0];
D=[0];

%comprobamos CONTROLABILIDAD del sistema
M_controlable=[B A*B A*A*B];
rank(M_controlable);%me devuelve un estimado del numero de columnas o filas
%linealmente independientes que tiene esta matriz como me da 3 esto quiere 
%decir que el sistema es controlable

%implementacion de funciones a usar
tf=2; dt=1*10^-6; t=0:dt:(tf-dt); periodo=0.6;%[seg]
torq=1.15*10^-3;

Ref=pi/2*square(2*pi*t/periodo);%funcion de referencia que varia entre pi/2 y -pi/2
TL=torq/2*square(2*pi*t/periodo)+torq/2;%Funcion torque que varia entre 0 y 1.15*10^-3
% plot(RT);
% hold on
% plot(TL);

%CALCULO DEL CONTROLADOR K
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos
Q=diag([1 1 1]); R=1;

H=[A -B*inv(R)*B' ; -Q -A'];
[vects,autovals]=eig(H);  %columnas de vects: autovectores
%Debo extraer solo los autovectores cuyos autovalores son negativos:
autovects_neg=[];
for i=1:1:length(autovals)
    if (real(autovals(i,i)))<0
        autovects_neg=[autovects_neg vects(:,i)];
    end
end    

%divido la matriz de autovectores en 2 matrices:
[filas,colums]=size(autovects_neg);
M=autovects_neg(1:(filas/2),:);
PM=autovects_neg((filas/2+1):filas,:);
P=real(PM*inv(M));

K_h=inv(R)*B'*P;
K=lqr(A,B,Q,R);
%comparo resultados de utilizar el metodo del calculo del hamiltoniano o
%utilizando la funcion lqr() y verifico que son iguales

%Ganancia de prealimentacion

G=-INV(C*inv(A-B*K)*B);











