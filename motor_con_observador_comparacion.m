clc; clear ; close all
%{
---------------------------------------------------------------------------
Considerar que no puede medirse la corriente y sólo pueda medirse el ángulo, por lo 
que debe implementarse un observador. Obtener la simulación en las mismas condiciones que en 
el punto anterior, y superponer las gráficas para comparar.
---------------------------------------------------------------------------
%}
%este es el sistema con las 3 variables de estado calculadas con el sistema
%observado
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
%Dejamos explicitas las matrices del sistema dual qeu es el "observado"
An=[A zeros(3,1); -C 0];
Bn=[B ; 0];
Cn=[1 0];
Ao=A';
Bo=C';
Co=B';

%CALCULO DEL CONTROLADOR K
%para el calculo del mismo se utiliza el metodo LQR para lo cual definimos
Q=diag([1 1/800000 1/9000000 10000]); R=0.01;
K4=lqr(An,Bn,Q,R);
K=K4(1:3);

K_i=-K4(4);



Qo=diag([1 1 1/1000]); Ro=1;
Ko=lqr(Ao,Bo,Qo,Ro);
%Ko=1e4*[-0.0001 0.7639 3.1623];


%implementacion de funciones a usar
tf=2; dt=1*10^-5; t=0:dt:(tf-dt); periodo=0.6;%[seg]
torq=1.15*10^-3;

Ref=pi/2*square(2*pi*t/periodo);%funcion de referencia que varia entre pi/2 y -pi/2
TL=torq/2*square(2*pi*t/periodo)+torq/2;%Funcion torque que varia entre 0 y 1.15*10^-3

%iteracion
n=round(tf/dt);
X=zeros(3,n);
X(1,1)=0; %ia inicial
X(2,1)=0; %tita inicial
X(3,1)=0; %w inicial
psi(1)=0; %psi inicial
Xcom=zeros(3,n);
Xcom(1,1)=0; %ia_com inicial
Xcom(2,1)=0; %tita_com inicial
Xcom(3,1)=0; %w_com inicial
psi(1)=0; %psi inicial
psicom(1)=0;
Xhat=zeros(3,n);
Xhat(1,1)=0; %ia_hat inicial
Xhat(2,1)=0; %tita_hat inicial
Xhat(3,1)=0; %wr_hat inicial
U(1)=0;
Ucom(1)=0;

%DEFINO CONDICIONES INICIALES

for i=1:1:n-1
    X_a=X(:,i);%[ia ; tita ; w]
    Xhat_a=Xhat(:,i);%[ia_hat ; tita_hat ; w_hat]
    Y=C*X_a;
    Yhat= Co*Xhat_a;
    err=(Y-Yhat);
    psi_p=Ref(i)-Y;
    psi(i+1)=psi(i)+psi_p*dt;
    U_a=-K(2:3)*X_a(2:3)+K(1)*Xhat_a(1)+K_i*psi(i+1);% U estimando parametros
    %U_a=-K*Xhat_a+K_i*psi(i+1); % U con las 3 variables de estado estimadas
    U=[U U_a];
    
    
    Xp_1=-RA/LAA*X_a(1)-Km/LAA*X_a(3)+1/LAA*U_a;  %ia_p
    Xp_2= X_a(3);                               %tita_p
    Xp_3=Ki/J*X_a(1)-Bm/J*X_a(3)-1/J*TL(i);%    %W_p
    Xp_a=[Xp_1 ; Xp_2 ; Xp_3];
    % Realizamos la integracion de euler y actualizamos matriz X
    Xf= X_a+ dt*Xp_a; 
    
    X(1,i+1)=Xf(1);
    X(2,i+1)=Xf(2);
    X(3,i+1)=Xf(3);
    %Observador
    
    Xhat_p=U_a*B+Ko'*err+A*Xhat_a;
    Xhatf=Xhat_a + dt*Xhat_p;
    Xhat(1,i+1)=Xhatf(1);
    Xhat(2,i+1)=Xhatf(2);
    Xhat(3,i+1)=Xhatf(3);
    
    %creamos la integracion sin observador a comparar
    Xcom_a=Xcom(:,i);%[ia_com ; tita_com ; w_com]
    Ycom=C*Xcom_a;
    psicom_p=Ref(i)-Ycom;
    psicom(i+1)=psicom(i)+psicom_p*dt;
    Ucom_a=-K*Xcom_a+K_i*psicom(i+1);
    Ucom=[Ucom Ucom_a];
    
    Xcom_p_1=-RA/LAA*Xcom_a(1)-Km/LAA*Xcom_a(3)+1/LAA*Ucom_a;  %ia_p
    Xcom_p_2= Xcom_a(3);                               %tita_p
    Xcom_p_3=Ki/J*Xcom_a(1)-Bm/J*Xcom_a(3)-1/J*TL(i);%    %W_p
    Xcom_p_a=[Xcom_p_1 ; Xcom_p_2 ; Xcom_p_3];
    % Realizamos la integracion de euler y actualizamos matriz X
    Xcom_f= Xcom_a+ dt*Xcom_p_a; 
    
    Xcom(1,i+1)=Xcom_f(1);
    Xcom(2,i+1)=Xcom_f(2);
    Xcom(3,i+1)=Xcom_f(3);    
   
end

figure
subplot(2,1,1);
plot(t,Ref,'k');
hold on
grid on
plot(t,Xcom(2,:),'g');title('angulo tita sin observador y con integrador');xlabel('tiempo[s]');ylabel('angulo[rad]');
plot(t,X(2,:),'r');title('angulo tita observado vs orginial y con integrador');xlabel('tiempo[s]');ylabel('angulo[rad]');legend('tita original','tita observado')
subplot(2,1,2);
grid on
plot(t,Xcom(1,:),'g');title('corriente ia sin observador y con integrador');xlabel('tiempo[s]');ylabel('angulo[rad]');
hold on
plot(t,X(1,:),'r');title('corriente ia observada vs original con integrador');xlabel('tiempo[s]');ylabel('angulo[rad]');legend('Ia original','Ia observada')


