%% Project 12

clc, clear, clf

% Data
R0 = 0.2 ; % [m]
R1 = 0.5 ; % [m]
Rs = 0.3 ; % [m]


% Parameter values:

f = 100 ; % from 100 to 4000 [Hz]

% Steps in radial direction
N = 30 ;
% Steps in angular direction
M = floor(N*2*pi/(R1-R0));

% Resulting Uniform Step Size
h = (R1 - R0)/(N) ;


% Creating first U matrix
tic
U_matrix = Create_U(N, M, f, 0) ;
toc

% Creating Coordinate System
r = R0:h:R1 ;

x_vec = [] ;
y_vec = [];
for theta = 0:h:2*pi
    x_vec = [x_vec, r'*cos(theta)] ;
    y_vec = [y_vec, r'*sin(theta)] ;
end


% Plotting u(x,y) in (x,y) space
figure(1)
mesh(x_vec,y_vec,U_matrix)
grid on
xlabel("x")
ylabel("y")
zlabel("u(x,y)")
title("u(x,y) with frequency f = " + f)


% Frequency Sweep
Ep_f = frequency_sweep(N, M, 0);

% Frequency vector
f_vector = Ep_f(:,1) ;
% Potential Energy vector
Ep_vector = Ep_f(:,2) ;


figure(4)
semilogy(f_vector, Ep_vector)
grid on
xlabel("Frequency (f)")
ylabel("Potential Energy (Ep)")
title("Potential Energy as a function of frequency")


% Extracting which frequency gives the maximum Potential Energy
[max_Ep, Index] = max(Ep_vector)   ;
max_f = f_vector(Index)  ;

% Solving system for the resonant frequency
U_matrix_resonant = Create_U(N, M, max_f, 0) ;

figure(3)
mesh(x_vec,y_vec, U_matrix_resonant)
grid on
xlabel("x")
ylabel("y")
zlabel("u(x,y)")
title("u(x,y) at resonant frequency f = " + max_f)



% Angular size of the hole, from theta = 0 to theta = theta_end
theta_end = pi/4 ;

U = Create_U(N, M, f, theta_end) ; 


% Frequency Sweep
Ep_f = frequency_sweep(N, M, theta_end);

% Frequency vector
f_vector = Ep_f(:,1) ;
% Potential Energy vector
Ep_vector = Ep_f(:,2) ;


figure(5)
semilogy(f_vector, Ep_vector)
grid on
xlabel("Frequency (f)")
ylabel("Potential Energy (Ep)")
title("Potential Energy as a function of frequency, hole in cavity")



% Sound Intensity level as a function of frequency
L_f_vec = [] ;
df = 100;

%
for f = 100:df:4000
    U_matrix = Create_U(N, M, f, theta_end) ;
    L = Sound_Intensity(U_matrix, M, theta_end) ;
    L_f_vec = [L_f_vec; f, L] ;
end

figure(6)
plot(L_f_vec(:,1), L_f_vec(:,2))
yline(0)
grid on
xlabel("Frequency [Hz]")
ylabel("Sound Intensity Level [dB]")
title("Sound Intensity Level at 10 [m] from the cavity")



% Functions

function L = Sound_Intensity(U, M, theta_end)
    % Data
    H = 0.05 ; % [m]
    R0 = 0.2 ; % [m]
    R1 = 0.5 ; % [m]
    c0 = 343 ; % [m/s]
    rho = 1.225 ; % [kg/m^3]
    
    h = 2*pi/M ;
    P = 0 ;
    U_R1 = U(end,:) ;
    
    theta_end_index = floor((M/(2*pi/theta_end))) ;
    
    U_R1 = U_R1(1:theta_end_index) ;
    for l = 1:theta_end_index
        
        P = P + ((R1*H)/(2*rho*c0))*(U_R1(l))^2 ;

        % Element from the matrix need to be U(R1, theta), might have to
        % look closer at line above
    end
    P = P*h ;
    d = 10; % distance from the cavity [m]

    I0 = 10^-12 ;
    Ir = P/(4*pi*d^2) ;
    
    % Sound Intensity level
    L = 10*log10(Ir/I0) ; % [dB]
end



function EP_f = frequency_sweep(N, M, theta_end)
    
    % Step Size
    df = 10;
   
    EP_f = [] ;
    for f = 1000:df:1300

        U_matrix = Create_U(N, M, f, theta_end) ;
        
        Ep = Potential_Energy(U_matrix, N, M) ;
        
        EP_f = [EP_f; f, Ep] ;
    end
end


function Ep = Potential_Energy(U, N, M)

    % Data
    H = 0.05 ; % [m]
    R0 = 0.2 ; % [m]
    R1 = 0.5 ; % [m]
    c0 = 343 ; % [m/s]
    rho = 1.225 ; % [kg/m^3]
    h = (R1 - R0)/(N) ;
    r = R0:h:R1 ;
    
    Ep = 0 ;
    
    for j = 1:N+1
        for l = 1:M
            Ep = Ep + (abs(U(j,l))^2)*r(j)*h^2 ;
        end
    end
    Ep = Ep*H/(2*rho*c0^2) ;
end

function U = Create_U(N, M, f, theta_end)

    if theta_end == 0 % No hole
        A = sparse(Create_A_Matrix(N, M, f)) ;
    else
        A = Create_A_Matrix_hole(N, M, f, theta_end) ;
    end
    
    g_vec = sparse(Create_g_vector(N, f)) ;
    U = sparse(A\g_vec) ;
    U  = sparse(reshape(U, N+1,M));
    U = [U,U(:,1)];
end



function g_vec = Create_g_vector(N, f)

    % Data
    H = 0.05 ; % [m]
    R0 = 0.2 ; % [m]
    R1 = 0.5 ; % [m]
    Rs = 0.3 ; % [m]

    c0 = 343 ; % [m/s]
    rho = 1.225 ; % [kg/m^3]
    omega = 0.04 ; % [m]
    q = 6*10^-7; %[m^3/s]
    k = 2*pi*f/c0 ; % [-]
    
    % Step size
    h = (R1 - R0)/(N) ;
    
    g = @(r, theta) rho*c0*k*q*((1/((omega^2)*H*pi))*exp((-1)*((Rs/omega)^2)*((1-r/Rs)^2 + (pi - theta)^2))) ;
    
    g_vec = sparse([]) ;
    
    % Creating discretized source function function
    for theta = 0:h:2*pi-h 
        for r = R0:h:R1 
            g_vec = [g_vec;g(r,theta)];
        end
    end  

end
function A = Create_A_Matrix(N, M, f)
    % Data
    c0 = 343 ; % [m/s]
    k = 2*pi*f/c0 ; % [-]
    R0 = 0.2 ; % [m]
    R1 = 0.5 ; % [m]
    
    % Step sizes
    h = (R1 - R0)/(N) ;
    r = sparse(R0:h:R1) ;

    % Sub Matrices
    a = sparse(1/h^2 + (1/(2*h))*(1./r)) ;
    b = sparse(-2/(h^2) -(2/(h^2))*(1./(r.^2)) + k^2) ;
    c = sparse(1/h^2 - (1/(2*h))*1./r) ;
    d = sparse(((1/h^2)*(1./(r.^2)))) ; 
    e = sparse(((1/h^2)*(1./(r.^2)))) ;

    % Diagonal Block
    A_Diagonal = sparse(diag(b) + diag(a(1:end-1),1) + diag(c(2:end),-1)) ;
    A_Diagonal(1,2) = a(1)+c(1) ;
    A_Diagonal(end,end-1) = c(end)+a(end)  ;
    % Creating big Matrix A
    A = sparse(kron(sparse(eye(M,M)),A_Diagonal))  ;

    % Adding additional elements and BC:s
    for j = 1:(N+1):(M-1)*(N+1)
        if j == (M-1)*(N+1)
            A(1:N, (N+1)*(M-1):((N+1)*(M-1) + N)) = sparse(2*diag(d)) ;
        else
            A(j:(j + N ), N + j + 1:(j + 1+ 2*(N ))) = sparse(diag(d)) ;
        end
        if j == 1
            A(((N+1)*(M-1)+1):((N+1)*(M-1)+ N)+1, 1:N+1) = sparse(2*diag(e)) ;
        else
            A(N + 1 + j :(j + 1+ 2*(N )), j:(j + N )) = sparse(diag(e)) ;
        end
    end

end
function A = Create_A_Matrix_hole(N, M, f, theta_end)

    % Data
    c0 = 343 ; % [m/s]
    k = 2*pi*f/c0 ; % [-]
    R0 = 0.2 ; % [m]
    R1 = 0.5 ; % [m]
    
    
    % Step sizes
    h = (R1 - R0)/(N) ;
    r = sparse(R0:h:R1) ;


    % Sub Matrices
    a = sparse(1/h^2 + (1/(2*h))*(1./r)) ;
    b = sparse(-2/(h^2) -(2/(h^2))*(1./(r.^2)) + k^2) ;
    c = sparse(1/h^2 - (1/(2*h))*1./r) ;
    d = sparse(((1/h^2)*(1./(r.^2)))) ; 
    e = sparse(((1/h^2)*(1./(r.^2)))) ;

    
    theta_end_index = floor((M/(2*pi/theta_end))) ;
    
    
    % Diagonal Block
    A_Diagonal = sparse(diag(b) + diag(a(1:end-1),1) + diag(c(2:end),-1)) ;
    A_Diagonal(1,2) = a(1) + c(1) ;
    A_Diagonal(end,end-1) = a(end) + c(end);
    
    
    A_diagonal_nohole = A_Diagonal;
    A_diagonal_hole = A_Diagonal;
    A_diagonal_hole(end-1,end) = c(end) + a(end)  ;
    A_diagonal_hole(end,end) = b(end) - 2*a(end)*h*1i*k ;


    A1 = zeros(M,M) ;
    A1(1:theta_end_index,1:theta_end_index) = diag(diag(ones(theta_end_index))) ;
    A2 = zeros(M,M) ;
    A2(theta_end_index+1:end,theta_end_index+1:end) = diag(diag(ones(M- theta_end_index))) ;
    
    A = kron(A1, A_diagonal_hole) + kron(A2, A_diagonal_nohole);
    

    % Adding additional elements and BC:s
    for j = 1:(N+1):(M-1)*(N+1)
        if j == (M-1)*(N+1)
            A(1:N, (N+1)*(M-1):((N+1)*(M-1) + N)) = 2*sparse(diag(d)) ;
        else
            A(j:(j + N ), N + j + 1:(j + 1+ 2*(N ))) = sparse(diag(d)) ;
        end
        if j == 1
            A(((N+1)*(M-1)+1):((N+1)*(M-1)+ N)+1, 1:N+1) = 2*sparse(diag(e)) ;
        else
            A(N + 1 + j :(j + 1+ 2*(N )), j:(j + N )) = sparse(diag(e)) ;
        end
    end
 
end
