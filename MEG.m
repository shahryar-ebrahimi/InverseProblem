
clear all
close all
clc

%% Shahryar Ebrahimi
%% S.N = 810196093
%% Constants

r3 = 9; rq = 7;  q0  = [0,0,1];
m  = 33; n = 65; miu = 4*pi*(10^-7); 

% =========================================================================
% rq Vectors

Rq      = zeros(n,3);
Rq(1,:) = [0,0,rq];

for i = 1:4
    
    teta = (i*pi)/8 ;
    
    for j = 0:15

       phi                = j*pi/8 ;
       [x,y,z]            = sph2cart(phi,(pi/2)-teta,rq);
       Rq(16*(i-1)+j+2,:) = [x,y,z];
        
    end

end

% =========================================================================
% R Vectors

R      = zeros(m,3);
R(1,:) = [0,0,r3];


for i = 1:4
    
    teta = (i*pi)/8;
    
    for j = 0:7

       phi              = j*pi/4;
       [x,y,z]          = sph2cart(phi,(pi/2)-teta,r3);
       R(8*(i-1)+j+2,:) = [x,y,z];
        
    end

end

%% Part 1

G = zeros(m,3*n);                      % Gain matrix

for i = 1:m
    for j = 1:n
        
        G(i,(j-1)*3+1:j*3) = 1*cross(R(i,:),Rq(j,:))/(norm(R(i,:))*(norm(Rq(j,:)-R(i,:)))^3);
        
    end
end


figure, plot(1:33,G(:,33),'r'); grid;
xlabel('Sensor Number'); ylabel('Gain of z elemet'); title('Column No.33 of G matrix');


%% Part 2 

i  = 45/22.5 ;
j  = 180/22.5 ;
Q  = zeros(n*3,1) ;
Q( (16*(i-1)+j+1)*3 +1 : (16*(i-1)+j+1)*3 +3 ) = q0 ;

b  = G*Q;

X  = R(:,1);
Y  = R(:,2);

figure, colormap('jet'); scatter(X,Y,36,b,'MarkerFaceColor','flat'); grid; colorbar;
xlabel('X'); ylabel('Y'); title('radial value of B');


X  = Rq(:,1);
Y  = Rq(:,2);
qm = zeros(n,1);

for i = 1:n
    
   qm(i) = norm(Q(3*(i-1)+1:3*(i-1)+3)); 
   
end

figure, colormap('jet'); scatter(X,Y,36,qm,'MarkerFaceColor','flat'); grid; colorbar;
xlabel('X'); ylabel('Y'); title('Amplitude of dipoles');


%% Part 3 

q = G'*( (G*G')^-1 )*b ;
J = zeros(n,1) ;

for i = 1:n
    
    J(i) = norm(q((i-1)*3+1:(i )*3));
    
end

X = Rq(:,1);
Y = Rq(:,2);

figure, colormap('jet'); scatter(X,Y,36,J,'MarkerFaceColor','flat'); grid; colorbar;
xlabel('X'); ylabel('Y'); title('Amplitude of dipoles by MN ');


%% Part 4

s     = mean(abs(b));
Sig   = s/(10^2);
noise = random('Normal',0,Sig,m,1);

B     = b + noise ; 

q     = G'*( (G*G')^-1 )*B ;
J     = zeros(n,1) ;

for i = 1:n
    
    J(i) = norm(q((i-1)*3+1:(i )*3));
    
end

X=Rq(:,1);
Y=Rq(:,2);


figure, colormap('jet'); scatter(X,Y,36,J,'MarkerFaceColor','flat'); grid; colorbar; 
xlabel('x'); ylabel('Y'); title('Amplitude of dipoles with noise by MN');


%% Part 5


lambda = [.001 .01 .1] ;     % arbitrary values for lambda

q1 = G'*((G*G'+lambda(1)*eye(m))^-1)*B;
q2 = G'*((G*G'+lambda(2)*eye(m))^-1)*B;
q3 = G'*((G*G'+lambda(3)*eye(m))^-1)*B;

J  = zeros(n,3);

for i = 1:n
    
    J(i,1) = norm(q1((i-1)*3+1:(i )*3)); 
    J(i,2) = norm(q2((i-1)*3+1:(i )*3));
    J(i,3) = norm(q3((i-1)*3+1:(i )*3));
    
end

X = Rq(:,1);
Y = Rq(:,2);

% =========================================================================
% Comparing the results ...
                    
% lambda = 0.001

figure, colormap('jet'); scatter(X,Y,36,J(:,1),'MarkerFaceColor','flat'); grid; colorbar;
xlabel('X'); ylabel('Y'); title('Amplitude of dipoles with noise by RMN with lambda=0.001');

% lambda = 0.01

figure, colormap('jet'); scatter(X,Y,36,J(:,2),'MarkerFaceColor','flat'); grid; colorbar;
xlabel('x'); ylabel('y'); title('Amplitude of dipoles with noise by RMN with lambda=0.01');

% for lambda = 0.1

figure, colormap('jet'); scatter(X,Y,36,J(:,3),'MarkerFaceColor','flat'); grid; colorbar;
xlabel('x'); ylabel('y'); title('Amplitude of dipoles with noise by RMN with lambda=0.1');
