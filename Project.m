%% Project 12

clc, clear

% Parameter values

f = 100 ; % from 100 to 400 [Hz]
N = 10 ;

% Data
H = 0.05 ; % [m]
R0 = 0.2 ; % [m]
R1 = 0.5 ; % [m]
Rs = 0.3 ; % [m]

c0 = 343 ; % [m/s]
rho = 1.225 ; % [kg/m^3]
omega = 0.04 ; % [m]
q = 6e-7; %[m^3/s]

k = 2*pi*f/c0 ; % [-]

% Functions

g = @(r, theta) rho*c0*k*q*delta_omega(r, theta) ;
delta_omega = @(r, theta)  (1/((omega^2)*H*pi))*exp((-1)*((Rs/omega)^2)*(((1-r)/Rs)^2) + (pi - theta)^2)  ;



A = Create_A_Matrix(N, k) ;

function A = Create_A_Matrix(N, k)

    R0 = 0.2 ; % [m]
    R1 = 0.5 ; % [m]
    
    M = round(N*2*pi/(R1-R0)) ;
    
    % Step sizes
    h = (R1 - R0)/(N) ;
    r = R0:h:R1 ;

    % Sub Matrices
    a = (1/h^2 + (1/(2*h))*(1./r)) ;
    b = -2/(h^2) -2/(h^2)*(1./(r.^2)) + k^2 ;
    c = 1/h^2 - 1/(2*h)*1./r ;
    d = (1/h^2)*(1./(r.^2)) ; 
    e = d ; 

    
    % Diagonal Block
    A_Diagonal = diag(b) + diag(a(1:end-1),1) + diag(c(2:end),-1) ;
    A_Diagonal(end,1) = a(end) ;
    A_Diagonal(1,end) = c(1) ;

    % Creating big Matrix A
    A = kron(eye(M+1,M+1),A_Diagonal)  ;

    % Adding additional elements and BC:s
    for j = 1:(N+1):(M)*N
        if j == 1 
            A(j:(j + N ), N + j + 1:(j + 1+ 2*(N ))) = 2*diag(d) ;
        else
            A(j:(j + N ), N + j + 1:(j + 1+ 2*(N ))) = diag(d) ;
        end
        if j == M*N
            A(N + j + 1:(j + 1+ 2*(N )), j:(j + N )) = 2*diag(e) ;
        else
            A(N + j + 1:(j + 1+ 2*(N )), j:(j + N )) = diag(e) ;
        end
    end
end
  