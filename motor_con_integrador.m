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
%EN ESTE CASO VAMOS A AGREGAR UN INTEGRADOR Y UNA GANANCIA K_I
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
%pero ahora como vamos a agregar integrador las matrices A y B se ven
%modificadas
An=[A zeros(3,1); -C 0];
Bn=[B ; 0];
Cn=[1 0];

%implementacion de funciones a usar
tf=2; dt=1*10^-5; t=0:dt:(tf-dt); periodo=0.6;%[seg]
torq=1.15*10^-3;

Ref=pi/2*square(2*pi*t/periodo);%funcion de referencia que varia entre pi/2 y -pi/2
TL=torq/2*square(2*pi*t/periodo)+torq/2;%Funcion torque que varia entre 0 y 1.15*10^-3

%CALCULO DEL CONTROLADOR K
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos
Q=diag([1 1/800000 1/9000000000000000 100000]); R=0.01;

Ka=lqr(An,Bn,Q,R);
K_i= -Ka(4);
K=Ka(1:3);

Kn=[K -K_i];

%iteracion
n=round(tf/dt);
X=zeros(3,n);
X(1,1)=0; %ia inicial
X(2,1)=0; %tita inicial
X(3,1)=0; %w inicial
psi(1)=0; %psi inicial
%DEFINO CONDICIONES INICIALES

for i=1:1:n-1
    X_a=[X(1,i); X(2,i) ; X(3,i)];%[ia ; tita ; w]
    psi_p=Ref(i)-C*X(:,i);
    psi(i+1)=psi(i)+psi_p*dt;
    U=-K*X_a+K_i*psi(i+1);
    
    Xp_1=-RA/LAA*X_a(1)-Km/LAA*X_a(3)+1/LAA*U;  %ia_p
    Xp_2= X_a(3);                               %tita_p
    Xp_3=Ki/J*X_a(1)-Bm/J*X_a(3)-1/J*TL(i);%    %W_p
    
    Xp_a=[Xp_1 ; Xp_2 ; Xp_3];
    
    Xf= X_a+ dt*Xp_a; % Realizamos la integracion de euler y actualizamos matriz X
    
    X(1,i+1)=Xf(1);
    X(2,i+1)=Xf(2);
    X(3,i+1)=Xf(3);
end
%ploteo de entrada con ganancia de prealimentacion U  y perturbacion TL
figure
subplot(2,1,1)
hold on; grid on;
plot(Ref);title('Referencia de Entrada');xlabel('tiempo[s]');ylabel('angulo');
subplot(2,1,2)
hold on;grid on;
plot(TL);title('Torque de perturbación');xlabel('Tiempo');ylabel('Torque');

figure


plot(t,Ref);
grid on
hold on

plot(t,X(2,:));title('angulo tita sin observador con integrador');xlabel('tiempo[s]');ylabel('angulo[rad]');legend('REF','var de estado')







