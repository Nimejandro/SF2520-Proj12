%% Project 12

clc, clear

% Parameter values

f = 1000 ; % from 100 to 4000 [Hz]
N = 20 ;

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

% g = @(r, theta) rho*c0*k*q*delta_omega(r, theta) ;
delta_omega = @(r, theta)  (1/((omega^2)*H*pi))*exp((-1)*((Rs/omega)^2)*(((1-r)/Rs)^2) + (pi - theta)^2)  ;


M = floor(N*2*pi/(R1-R0)) ;


A = sparse(Create_A_Matrix(N, k)) ;

g_vec = sparse(Create_g_vector(N, f)) ;



U_matrix = sparse( A\g_vec) ;
U_matrix  = sparse(reshape(U_matrix, N+1,M+2));

h = (R1 - R0)/(N) ;

r = sparse(R0:h:R1) ;
r = r' ;

x_vec = [] ;
y_vec = [];


for theta = 0:h:2*pi + h
    x_vec = [x_vec, r*cos(theta)] ;
    y_vec = [y_vec, r*sin(theta)] ;
end



figure(1)
mesh(x_vec,y_vec,U_matrix)
xlabel("x")
ylabel("y")
zlabel("u(x,y)")




function g_vec = Create_g_vector(N, f)

    g_vec = sparse([]) ;
    
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
    
    M = floor(N*2*pi/(R1-R0)) ;
    
    % Step sizes
    h = (R1 - R0)/(N) ;
    
    
    g = @(r, theta) rho*c0*k*q*((1/((omega^2)*H*pi))*exp((-1)*((Rs/omega)^2)*(((1-r)/Rs)^2) + (pi - theta)^2)) ;
    
    
    for theta = 0:h:2*pi + h
        for r = R0:h:R1 
            g_vec = [g_vec;g(r,theta)];
        end
    end   
end

function A = Create_A_Matrix(N, k)

    R0 = 0.2 ; % [m]
    R1 = 0.5 ; % [m]
    
    M = floor(N*2*pi/(R1-R0)) ;
    
    % Step sizes
    h = (R1 - R0)/(N) ;
    r = sparse(R0:h:R1) ;

    % Sub Matrices
    a = sparse((1/h^2 + (1/(2*h))*(1./r))) ;
    b = sparse(-2/(h^2) -2/(h^2)*(1./(r.^2)) + k^2) ;
    c = sparse(1/h^2 - 1/(2*h)*1./r) ;
    d = sparse(((1/h^2)*(1./(r.^2)))) ; 
    e = sparse(d) ; 

    
    % Diagonal Block
    A_Diagonal = sparse(diag(b) + diag(a(1:end-1),1) + diag(c(2:end),-1)) ;
    A_Diagonal(end,1) = a(end) ;
    A_Diagonal(1,end) = c(1) ;

    % Creating big Matrix A
    A = sparse(kron(sparse(eye(M+1,M+1)),A_Diagonal))  ;

    % Adding additional elements and BC:s
    for j = 1:(N+1):(M+1)*(N+1)
        if j == 1 
            A(j:(j + N ), N + j + 1:(j + 1+ 2*(N ))) = sparse(2*diag(d)) ;
        else
            A(j:(j + N ), N + j + 1:(j + 1+ 2*(N ))) = sparse(diag(d)) ;
        end
        if j == (M+1)*(N+1)
            A(N + j + 1:(j + 1+ 2*(N )), j:(j + N )) = sparse(2*diag(e)) ;
        else
            A(N + j + 1:(j + 1+ 2*(N )), j:(j + N )) = sparse(diag(e)) ;
        end
    end
end
  