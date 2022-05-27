clc; clear ; close all
%{
---------------------------------------------------------------------------
Para el caso del avión, emplear un tiempo de integración por Euler adecuado y un tiempo 
de simulación de 70seg. Los parámetros son a=0.05; w=3; b=5; c=100, hallar un controlador para 
que los polos de lazo cerrado se ubican en ui=-15?15j; -0.5?0.5j, para referencias de 100 y -100 
metros en altura, ambas con alturas iniciales de -500 y 500.
---------------------------------------------------------------------------
%}
%DEFINO PARAMETROS
w=3; a=0.05; b=5; c=100;

%DEFINO MATRICES
%X=[alfa phi phi_p h]
A = [ -a,a,0,0 ; 0,0,1,0 ; w^2,-w^2,0,0 ; c,0,0,0 ];
B = [ 0; 0 ; w^2*b ; 0 ];
C = [ 0,0,0,1 ];
D = 0;

%comprobamos CONTROLABILIDAD del sistema
M_controlable=[B A*B A*A*B A*A*A*B];
rank(M_controlable)% nos da 4 por lo tanto es controlable

%CALCULO DEL CONTROLADOR K
%en este caso tenemos el valor de los polos deseados y en base a esto
%obtenemos el valor del controlador por lo tanto utilizando la formula de
%akerman lo impementamos

p1=-15+15i;
p2=-15-15i;
p3=-0.5+0.5i;
p4=-0.5-0.5i;

%syms p1 p2 p3 p4 s
%expand((s-p1)*(s-p2)*(s-p3)*(s-p4))
a0=1;
a1=-p1-p2-p3-p4;
a2=p1*p2+p1*p3+p1*p4+p2*p3+p2*p4+p3*p4;
a3=-p1*p2*p3-p1*p2*p4-p1*p3*p4-p2*p3*p4;
a4=p1*p2*p3*p4;
phi_A=a0*A^4+a1*A^3+a2*A^2+a3*A^1+a4*A^0;

AUX=[B A*B A*A*B A*A*A*B];
K=[0 0 0 1]*inv(AUX)*phi_A;

%podemos verificar con la funcion place que coloca los polos deseados y nos
%devuelve el controlador buscado
polos=[p1 p2 p3 p4];
K_verif=place(A, B, polos);

%Ganancia de prealimentacion
G=-inv(C*inv(A-B*K)*B);

%simulacion Para altura inicial 
tf=70; dt=1*10^-3; t=0:dt:(tf-dt); 
h_ref=4000;     %altura de referencia
h_o=-500;       %Altura inicial
referencia=100; %Valor de establecimiento
n=round(tf/dt);
X=zeros(4,n);
X_inicial=[0 0 0 h_o]; %[alfa phi phi_p h] ho=4000
X(1,1)=X_inicial(1);   %alfa   inicial
X(2,1)=X_inicial(2);   %phi    inicial
X(3,1)=X_inicial(3);   %phi_p  inicial
X(4,1)=X_inicial(4);   %h      inicial
Ua1(1)=0;
%DEFINO CONDICIONES INICIALES
Ref=referencia*ones(1,n);% altura de referencia a la cual se pretende llegar (altura referenciada al ho)
for i=1:1:n-1
    X_a=X(:,i);%[alfa ; phi ; phi_p ; h]
    U=-K*X_a+Ref(i)*G;
    Ua1=[Ua1 U];
    Xp_1=a*(X_a(2)-X_a(1));         %alfa_p
    Xp_2= X_a(3);                   %phi_p
    Xp_3=-w^2*(X_a(2)-X_a(1)-b*U);  %phi_pp
    Xp_4=c*X_a(1);                  %h_p
    
    Xp_a=[Xp_1 ; Xp_2 ; Xp_3 ; Xp_4];
    
    Xf= X_a+ dt*Xp_a; % Realizamos la integracion de euler y actualizamos matriz X
    
    X(:,i+1)=Xf;
end
figure
subplot(2,2,1);
plot(t,Ref+h_ref,'b');
grid on
hold on
plot(t,X(4,:)+h_ref,'g');title('altura h con h_o=-500 y ref=100');xlabel('tiempo[s]');ylabel('posicion[m]');legend('Referencia','altura')
subplot(2,2,3);
grid on
hold on
plot(t,Ua1,'r');title('accion de control caso 1')


%Simulamos el otro caso que tiene referencia de -100 y arranca en 500
h_ref=4000;         %atura de referencia
h_o=500;            %Altura inicial
referencia=-100;    %Valor de establecimiento
n=round(tf/dt);
X=zeros(4,n);
X_inicial=[0 0 0 h_o]; %[alfa phi phi_p h] ho=4000
X(1,1)=X_inicial(1);   %alfa   inicial
X(2,1)=X_inicial(2);   %phi    inicial
X(3,1)=X_inicial(3);   %phi_p  inicial
X(4,1)=X_inicial(4);   %h      inicial
Ua2(1)=0;
%DEFINO CONDICIONES INICIALES
Ref=referencia*ones(1,n);% altura de referencia a la cual se pretende llegar (altura referenciada al ho)
for i=1:1:n-1
    X_a=X(:,i);%[alfa ; phi ; phi_p ; h]
    U=-K*X_a+Ref(i)*G;
    Ua2=[Ua2 U];
    Xp_1=a*(X_a(2)-X_a(1));         %alfa_p
    Xp_2= X_a(3);                   %phi_p
    Xp_3=-w^2*(X_a(2)-X_a(1)-b*U);  %phi_pp
    Xp_4=c*X_a(1);                  %h_p
    
    Xp_a=[Xp_1 ; Xp_2 ; Xp_3 ; Xp_4];
    
    Xf= X_a+ dt*Xp_a; % Realizamos la integracion de euler y actualizamos matriz X
    
    X(:,i+1)=Xf;
end

subplot(2,2,2);
plot(t,Ref+h_ref,'b');
grid on
hold on
plot(t,X(4,:)+h_ref,'g');title('altura h con h_o=500 y ref=-100');xlabel('tiempo[s]');ylabel('posicion[m]');legend('Referencia','altura')
subplot(2,2,4);
grid on
hold on
plot(t,Ua2,'r');title('accion de control caso 2')





