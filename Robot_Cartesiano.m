%Limpieza de pantalla
clear all
close all
clc

tic
%Declaración de variables simbólicas
syms l1(t) l2(t) l3(t) t
syms l1p(t) l2p(t) l3p(t)    %Velocidades de cada articulación
syms l1pp(t) l2pp(t) l3pp(t)  %Aceleraciones de cada articulación
syms m1 m2 m3 Ixx1 Iyy1 Izz1 Ixx2 Iyy2 Izz2 Ixx3 Iyy3 Izz3  %Masas y matrices de Inercia
syms lc1 lc2 lc3 %lc=distancia al centro de masa de cada eslabón
syms pi g a cero

%Creamos el vector de coordenadas articulares
Q= [l1; l2; l3];
%disp('Coordenadas generalizadas');
%pretty (Q);

%Creamos el vector de velocidades articulares
Qp= [l1p; l2p; l3p];
%disp('Velocidades generalizadas');
%pretty (Qp);
%Creamos el vector de aceleraciones articulares
Qpp= [l1pp; l2pp; l3pp];
%disp('Aceleraciones generalizadas');
%pretty (Qpp);

%Configuración del robot, 0 para junta rotacional, 1 para junta prismática
RP=[1 1 1];

%Número de grado de libertad del robot
GDL= size(RP,2);
GDL_str= num2str(GDL);

%Articulación 1 
%Posición de la articulación 1 respecto a 0
P(:,:,1)= [0; 0; l1];
%Matriz de rotación de la junta 1 respecto a 0.... 
R(:,:,1)= [1 0 0; %Rotx90
           0 0 -1;
           0 1 0];

%Articulación 2
%Posición de la articulación 2 respecto a 1
P(:,:,2)= [0; 0; l2];
%Matriz de rotación de la junta 1 respecto a 0
R(:,:,2)= [0 0 1; %Roty90
           0 1 0;
          -1 0 0];

%Articulación 3
%Posición de la articulación 3 respecto a 2
P(:,:,3)= [0; 0; l3];
%Matriz de rotación de la junta 1 respecto a 0
R(:,:,3)= [1 0 0; %Identidad
           0 1 0;
           0 0 1];


%Creamos un vector de ceros
Vector_Zeros= zeros(1, 3);

%Inicializamos las matrices de transformación Homogénea locales
A(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las matrices de transformación Homogénea globales
T(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las posiciones vistas desde el marco de referencia inercial
PO(:,:,GDL)= P(:,:,GDL); 
%Inicializamos las matrices de rotación vistas desde el marco de referencia inercial
RO(:,:,GDL)= R(:,:,GDL); 


for i = 1:GDL
    i_str= num2str(i);
   %disp(strcat('Matriz de Transformación local A', i_str));
    A(:,:,i)=simplify([R(:,:,i) P(:,:,i); Vector_Zeros 1]);
   %pretty (A(:,:,i));

   %Globales
    try
       T(:,:,i)= T(:,:,i-1)*A(:,:,i);
    catch
       T(:,:,i)= A(:,:,i);
    end
%     disp(strcat('Matriz de Transformación global T', i_str));
    T(:,:,i)= simplify(T(:,:,i));
%     pretty(T(:,:,i))

    RO(:,:,i)= T(1:3,1:3,i);
    PO(:,:,i)= T(1:3,4,i);
    %pretty(RO(:,:,i));
    %pretty(PO(:,:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CALCULAMOS LAS VELOCIDADES PARA CADA ESLABÓN%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%  VELOCIDADES PARA ESLABÓN 3  %%%%%%%%%%

%Jacobiano eslabón 1
Jv_a1(:,1)=PO(:,:,1);
Jw_a1(:,1)=PO(:,:,1);

for k=1:1
    if RP(k)==0
        try
            Jv_a1(:,k)=cross(RO(:,3,k-1),PO(:,:,1)-PO(:,:,k-1));
            Jw_a1(:,k)=RO(:,3,k-1);
        catch
            Jv_a1(:,k)=cross([0;0;1],PO(:,:,1));
            Jw_a1(:,k)=[0;0;1];
        end
    else
        try
            Jv_a1(:,k)=RO(:,3,k-1);
        catch
            Jv_a1(:,k)=[0;0;1];
        end
        Jw_a1(:,k)=[0;0;0];
    end
end

Jv_a1=simplify(Jv_a1);
Jw_a1=simplify(Jw_a1);

Jac1=[Jv_a1;Jw_a1];

disp('Velocidad lineal eslabón 1')
Qp=Qp(t);
V1=simplify(Jv_a1*Qp(1:1));
pretty(V1)

disp('Velocidad angular eslabón 1')
W1=simplify(Jw_a1*Qp(1:1));
pretty(W1)

%Jacobiano eslabón 2
Jv_a2(:,2)=PO(:,:,2);
Jw_a2(:,2)=PO(:,:,2);

for k=1:2
    if RP(k)==0
        try
            Jv_a2(:,k)=cross(RO(:,3,k-1),PO(:,:,2)-PO(:,:,k-1));
            Jw_a2(:,k)=RO(:,3,k-1);
        catch
            Jv_a2(:,k)=cross([0;0;1],PO(:,:,2));
            Jw_a2(:,k)=[0;0;1];
        end
    else
        try
            Jv_a2(:,k)=RO(:,3,k-1);
        catch
            Jv_a2(:,k)=[0;0;1];
        end
        Jw_a2(:,k)=[0;0;0];
    end
end

Jv_a2=simplify(Jv_a2);
Jw_a2=simplify(Jw_a2);

Jac2=[Jv_a2;Jw_a2];

disp('Velocidad lineal eslabón 2')
V2=simplify(Jv_a2*Qp(1:2));
pretty(V2)

disp('Velocidad angular eslabón 2')
W2=simplify(Jw_a2*Qp(1:2));
pretty(W2)

%Jacobiano eslabón 3
Jv_a3(:,3)=PO(:,:,3);
Jw_a3(:,3)=PO(:,:,3);

for k=1:3
    if RP(k)==0
        try
            Jv_a3(:,k)=cross(RO(:,3,k-1),PO(:,:,3)-PO(:,:,k-1));
            Jw_a3(:,k)=RO(:,3,k-1);
        catch
            Jv_a3(:,k)=cross([0;0;1],PO(:,:,3));
            Jw_a3(:,k)=[0;0;1];
        end
    else
        try
            Jv_a3(:,k)=RO(:,3,k-1);
        catch
            Jv_a3(:,k)=[0;0;1];
        end
        Jw_a3(:,k)=[0;0;0];
    end
end

Jv_a3=simplify(Jv_a3);
Jw_a3=simplify(Jw_a3);

Jac3=[Jv_a3;Jw_a3];

disp('Velocidad lineal eslabón 3')
V3=simplify(Jv_a3*Qp(1:3));
pretty(V3)

disp('Velocidad angular eslabón 3')
W3=simplify(Jw_a3*Qp(1:3));
pretty(W3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Energía Cinética
%%%%%%%%%%%%%%%%%%%%%%%%%%Omitimos la división de cada lc%%%%%%%%%%%%%%%
%Distancia del origen del eslabón a su centro de masa
%Vectores de posición respecto al centro de masa
P01=subs(P(:,:,1)/2, l1, lc1);%La función subs sustituye l1 por lc1 en 
P12=subs(P(:,:,2)/2, l2, lc2); %la expresión P(:,:,1)/2
P23=subs(P(:,:,3)/2, l3, lc3);

%Creamos matrices de inercia para cada eslabón

I1=[Ixx1 0 0; 
    0 Iyy1 0; 
    0 0 Izz1];

I2=[Ixx2 0 0; 
    0 Iyy2 0; 
    0 0 Izz2];

I3=[Ixx3 0 0; 
    0 Iyy3 0; 
    0 0 Izz3];

%FUNCIÓN DE ENERGÍA CINÉTICA
%Calculamos la energía cinética para cada uno de los eslabones

%Eslabón 1
V1_Total= V1+cross(W1,P01);
K1= (1/2*m1*(V1_Total))'*((V1_Total)) + (1/2*W1)'*(I1*W1);
%disp('Energía Cinética en el Eslabón 1');
K1= simplify (K1);
%pretty (K1);

%Eslabón 2
V2_Total= V2+cross(W2,P12);
K2= (1/2*m2*(V2_Total))'*((V2_Total)) + (1/2*W2)'*(I2*W2);
%isp('Energía Cinética en el Eslabón 2');
K2= simplify (K2);
%pretty (K2);

%Eslabón 3
V3_Total= V3+cross(W3,P23);
K3= (1/2*m3*(V3_Total))'*((V3_Total)) + (1/2*W3)'*(I3*W3);
%isp('Energía Cinética en el Eslabón 3');
K3= simplify (K3);
%pretty (K3);

K_Total= simplify (K1+K2+K3);
disp('Energía Cinética Total');
pretty (K_Total);

