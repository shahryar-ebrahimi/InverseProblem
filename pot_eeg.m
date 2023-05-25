

%% Shahryar Ebrahimi
%% S.N = 810196093
%% pot_eeg function --> returning Lx,Ly,Lz

function [Lx, Ly, Lz] = pot_eeg(Rm,Rs,Sig,Sig_SK,r1,r2,r3,nLim,type)

    [x,y,z] = sph2cart(Rm(3),Rm(2),Rm(1));
    Rmc     = [x,y,z];

    [x,y,z] = sph2cart(Rs(3),Rs(2),Rs(1));
    Rsc     = [x,y,z];

    deltaR  = Rsc - Rmc + [0 0 Rmc(3)] ;

    [phi,theta,~] = cart2sph(deltaR(1),deltaR(2),deltaR(3));
    theta         = pi/2 -theta;

%==========================================================================

    LP  = zeros(nLim,1);
    ALP = zeros(nLim,1);

    for i = 1:nLim
    
        p      = legendre(i,cos(theta));
        LP(i)  = p(1);
        ALP(i) = p(2);
   
    end


%==========================================================================

    zita = Sig_SK/Sig ;
    f1   = r1/r3 ;
    f2   = r2/r3 ;

%==========================================================================

    dn   = zeros(nLim,1);
    
    for i = 1:nLim
        
        dn(i) = ((i+1)*zita+i)+(1-zita)*((i+1)*zita + i)*(f1^(2*i+1)-f2^(2*i+1))-i*(1-zita)^2*(f1/f2)^(2*i+1) ;
    
    end
    
%==========================================================================

b  = Rm(1)/r3;
N  = cumsum(ones(nLim,1));
Lz = (zita/(4*pi*Sig))*sum(LP.*N.*b.^(N-1).*(2*N+1)./(dn.*(N+1).*N));
Lx = (zita*cos(phi)/(4*pi*Sig))*sum(ALP.*b.^(N-1).*(2*N+1).^3./(dn.*(N+1).*N));
Ly = (zita*sin(phi)/(4*pi*Sig))*sum(ALP.*b.^(N-1).*(2*N+1)./(dn.*(N+1).*N));

end



