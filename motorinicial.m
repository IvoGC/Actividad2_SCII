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
tf=1.5; dt=1*10^-5; t=0:dt:(tf-dt); periodo=0.6;%[seg]
torq=1.15*10^-3;

Ref=pi/2*square(2*pi*t/periodo);%funcion de referencia que varia entre pi/2 y -pi/2
TL=torq/2*square(2*pi*t/periodo)+torq/2;%Funcion torque que varia entre 0 y 1.15*10^-3
% plot(Ref);
% hold on
% plot(TL);

%CALCULO DEL CONTROLADOR K
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos
%Q=diag([1/10000 1/7000000000 1/50]); R=1;
Q=diag([1 10/1.5e-1 1/1.5e2]); R=0.01;
H=[A -B*inv(R)*B' ; -Q -A'];
[V,D]=eig(H);  %columnas de vects: autovectores
%Debo extraer solo los autovectores cuyos autovalores son negativos:
Mx1x2=[];
for i=1:1:length(D)
    if (real(D(i,i)))<0
        Mx1x2=[Mx1x2 V(:,i)];
    end
end    

%divido la matriz de autovectores en 2 matrices:
[filas,colums]=size(Mx1x2);
M=Mx1x2(1:(filas/2),:);
PM=Mx1x2((filas/2+1):filas,:);
P=real(PM*inv(M));

K_h=inv(R)*B'*P;
K=lqr(A,B,Q,R);
%comparo resultados de utilizar el metodo del calculo del hamiltoniano o
%utilizando la funcion lqr() y verifico que son iguales

%Ganancia de prealimentacion

G=-inv(C*inv(A-B*K)*B);

%iteracion
n=round(tf/dt);
X=zeros(3,n);
X(1,1)=0; %ia inicial
X(2,1)=0; %tita inicial
X(3,1)=0; %w iicial
%DEFINO CONDICIONES INICIALES

for i=1:1:n-1
    X_a=[X(1,i); X(2,i) ; X(3,i)];%[ia ; tita ; w]
    U=-K*X_a+Ref(i)*G;
    
    Xp_1=-RA/LAA*X_a(1)-Km/LAA*X_a(3)+1/LAA*U;  %ia_p
    Xp_2= X_a(3);                               %tita_p
    Xp_3=Ki/J*X_a(1)-Bm/J*X_a(3)-1/J*TL(i);%    %W_p
    
    Xp_a=[Xp_1 ; Xp_2 ; Xp_3];
    
    Xf= X_a+ dt*Xp_a; % Realizamos la integracion de euler y actualizamos matriz X
    
    X(1,i+1)=Xf(1);
    X(2,i+1)=Xf(2);
    X(3,i+1)=Xf(3);
    
    %X(:,i+1)=Xf;
end
%ploteo de entrada con ganancia de prealimentacion U  y perturbacion TL
figure
subplot(2,1,1)
hold on; grid on;
plot(t,Ref);title('Referencia de Entrada');xlabel('tiempo[s]');ylabel('angulo');
subplot(2,1,2)
hold on;grid on;
plot(t,TL);title('Torque de perturbación');xlabel('Tiempo');ylabel('Torque');

figure

subplot(2,1,1)
plot(t,Ref);
grid on
hold on
plot(t,X(2,:),'r');title('angulo tita sin observador');xlabel('tiempo[s]');ylabel('angulo[rad]');legend('ref','var de estado')
subplot(2,1,2)
plot(t,X(1,:),'r');title('corriente sin observador');xlabel('tiempo[s]');ylabel('angulo[rad]');legend('var de estado')


