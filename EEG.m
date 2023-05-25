
clear all
close all
clc


%% Shahryar Ebrahimi
%% S.N = 810196093
%% Constants

r1   = 8;   r2     = 8.5;  r3  = 9; 
m    = 33;  rq     = 7;    n   = 65; 
nLim = 50 ; sig_sk = 1 ;   sig = 1; 

q0   = [0 0 1] ;

% =====================================================================
% rq Vectors

Rq      = zeros(n,3);
Rq(1,:) = [0,pi/2,rq];


for i=1:4
    
       theta = (i*pi)/8 ;
       
       for j=0:15

            phi                = j*pi/8 ;
            Rq(16*(i-1)+j+2,:) = [phi,(pi/2-theta),rq] ;
        
       end

end

% =====================================================================
% R vectors

R      = zeros(m,3);
R(1,:) = [0,pi/2,r3];

for i=1:4
    
       theta = (i*pi)/8 ;
       
       for j=0:7

            phi              = j*pi/4 ;
            R(8*(i-1)+j+2,:) = [phi,(pi/2 - theta),r3];
        
       end

end

%% Part 7

L = zeros(m,3*n);

for i = 1:m
    for j = 1:n
        
        [Lx, Ly, Lz]       = pot_eeg(Rq(j,:),R(i,:),sig,sig_sk,r1,r2,r3,nLim);
        L(i,(j-1)*3+1:j*3) = [Lx,Ly,Lz] ;
        
    end
end


figure, plot(1:1/3:size(L,2)/3+2/3,L(23,:),'m'), grid;
xlabel('dipole number (x:+0 , y:+1/3, z:+2/3 )'); ylabel('23rd sensor`s gain');   title('23rd row of the L matrix');


%% Part 8

Q = zeros(n*3,1);
i = 45/22.5 ;
j = 180/22.5 ;
Q((16*(i-1)+j+1)*3+1:(16*(i-1)+j+1)*3+3) = q0 ;

v = L*Q ;                                         % the main equation of EEG problem

X = r3.*cos(R(:,1)).*cos(R(:,2));
Y = r3.*sin(R(:,1)).*cos(R(:,2));

figure, colormap('jet'); scatter(X,Y,36,v,'MarkerFaceColor','flat'); grid ; colorbar ;
xlabel('X');  ylabel('Y');  title('measured V');


%% Part 9

q = pinv(L)*v ;
J = zeros(n,1);

for i = 1:n   
    
    J(i) = norm(q((i-1)*3+1:(i )*3));
    
end

X = rq.*cos(Rq(:,1)).*cos(Rq(:,2));
Y = rq.*sin(Rq(:,1)).*cos(Rq(:,2));


figure, colormap('jet'); scatter(X,Y,36,J,'MarkerFaceColor','flat') ;  grid; colorbar;
xlabel('X'); ylabel('Y'); title('Amplitude of dipoles by MN');

%% Part 10

fs = 500 ;
E1 = triang(fs)' ;
E2 = [ones(1,fs/2),-1*ones(1,fs/2)] ;

Q  = zeros(3*n,fs);

i  = 45/22.5 ;
j  = 180/22.5 ;
Q( (16*(i-1)+j+1)*3 +1 : (16*(i-1)+j+1)*3 +3 , :) = q0'*E1 ;

i  = 45/22.5 ;
j  = 90/22.5 ;
Q( (16*(i-1)+j+1)*3 +1 : (16*(i-1)+j+1)*3 +3 , :) = q0'*E2 ;

v = L*Q ;

% =========================================================================
% noise adding

s     = mean(abs(v(:,1))) ;
sig   = s/(10^2.5) ;
noise = random('Normal',0,sig,m,fs) ; 

V     = v + noise ;

% =========================================================================
% MUSIC

[U,S,v] = svd(V) ; 
lambda  = diag(S)/S(1,1) ;
p       = sum(lambda>.01);
Pperp   = (U(:,p+1:end)) * (U(:,p+1:end))' ;

J       = zeros(1,3*n);

for i = 1:3*n
    
   J(i) = norm(Pperp*L(:,i))/norm(L) ; 
    
end

[~,idx] = sort(J,'ascend');
dipole  = (floor((idx-1)/3)+1);
p_idx   = dipole(1:p);

A       = zeros(m,3*p);

for i = 1:p
    
    A(:, 3*(i-1)+1:3*(i-1)+3 ) = L(:, 3*(p_idx(i)-1)+1:3*(p_idx(i)-1)+3 ) ;
    
end

q = pinv(A)*V ;

J = zeros(n,fs);

for i = 1:p
    for t = 1:fs
        
        J(p_idx(i),t) = norm(q(3*(i-1)+1:3*(i-1)+3,t));  
        
    end
end


figure, colormap('jet'); scatter(X,Y,36,J(:,fs*.25),'MarkerFaceColor','flat'); grid; colorbar;
xlabel('X'); ylabel('Y'); title('Amplitude of dipoles by MUSIC t=0.25');


figure, plot(1/fs:1/fs:1 , J(p_idx(1) , :) , 'm' , 1/fs:1/fs:1 , J(p_idx(2) , :) , 'k'); grid;
xlabel('time (sec)');ylabel('Amplitude of active dipoles');

