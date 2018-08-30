function Res = Radial_Load_Flow(Sistema);
% Flujo de carga radial para sistemas de distribucion
% balanceados.  El sistema debe estar ordenado
% Sistema:
%     N1 N2 R X B/2 P Q alfa
NL = length(Sistema(:,1));  % numero de lineas
NN = NL + 1;
S = zeros(NN,1);
P0 = zeros(NN,1);
Q0 = zeros(NN,1);
alp = zeros(NN,1);
P0(2:NN,1) = Sistema(:,6);
Q0(2:NN,1)= Sistema(:,7);
alp(2:NN,1) = Sistema(:,8);
V = ones(NN,1);
err = 1000;
% barrido iterativo
iter = 0;
while (err>1E-10)
    Va = V;
    Pn = P0.*(abs(Va)).^alp;
    Qn = Q0.*(abs(Va)).^alp;
    S = Pn+j*Qn;
% corrientes en las cargas
    Ix = conj(S./V);
% efecto capacitivo de las lineas
    for k = NL:-1:1
	    N1 = Sistema(k,1);
        N2 = Sistema(k,2);
        Ix(N1) = Ix(N1) + j*Sistema(k,5)*(V(N1));
        Ix(N2) = Ix(N2) + j*Sistema(k,5)*(V(N2));
    end
	
% barrido de corrientes
    for k = NL:-1:1
        N1 = Sistema(k,1);
        N2 = Sistema(k,2);
        Ix(N1) = Ix(N1) + Ix(N2) ;
    end
% barrido de voltajes
    Perdidas = 0;
    for k = 1:NL
        N1 = Sistema(k,1);
        N2 = Sistema(k,2);
        ZX = Sistema(k,3) + j*Sistema(k,4);
        V(N2) = V(N1) - ZX.*Ix(N2);
        Perdidas = Perdidas + Sistema(k,3)*abs(Ix(N2))^2;
    end	
    err = max(abs(Va-V));
    iter = iter + 1;
    if(iter>50)
      disp('Error'); 
      err = 0;
    end
end
%% calculo de las perdidas
Res.Vs = V;
Res.Is = Ix;
Res.Pper = Perdidas;
Res.iter = iter;
