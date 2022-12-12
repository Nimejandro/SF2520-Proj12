%% Project 12

clc, clear, clf

% Data
R0 = 0.2 ; % [m]
R1 = 0.5 ; % [m]

% Parameter values:
f = 1000 ; % from 100 to 4000 [Hz]
N = 20;
theta_end = pi/4;

M = floor(N*2*pi/(R1-R0)) ;

% Resulting Uniform Step Size
h = (R1 - R0)/(N) ;

% Solving system for given frequency
U_matrix = Create_U(N, M, f) ;

% Creating Coordinate System
r = R0:h:R1 ;

x_vec = [] ;
y_vec = [] ;
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

%% Frequency Sweep
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
title("Potential Energy as a function of frequency, no hole")


%% Resonant

% Extracting which frequency gives the maximum Potential Energy
[max_Ep, Index] = max(Ep_vector)   ;
max_f = f_vector(Index)  ;

% Solving system for the resonant frequency
U_matrix_resonant = Create_U(N, M, max_f + 2) ;
 

figure(3)
mesh(x_vec,y_vec,U_matrix_resonant)
grid on
xlabel("x")
ylabel("y")
zlabel("u(x,y)")
title("u(x,y) close to resonant frequency f = " + (max_f + 2))



%% Hole

U_hole = real(Create_U_hole(N, M, f, theta_end)) ;


figure(5)
mesh(x_vec,y_vec,U_hole)
grid on
xlabel("x")
ylabel("y")
zlabel("u(x,y)")
title("u(x,y) with frequency f = " + f + " hole of size: " + theta_end)


%% Frequency Sweep
Ep_f = frequency_sweep(N, M, theta_end);

% Frequency vector
f_vector = Ep_f(:,1) ;
% Potential Energy vector
Ep_vector = Ep_f(:,2) ;

figure(6)
semilogy(f_vector, Ep_vector)
grid on
xlabel("Frequency (f)")
ylabel("Potential Energy (Ep)")
title("Potential Energy as a function of frequency, hole size: " + theta_end)


%% Resonant 

% Extracting which frequency gives the maximum Potential Energy
[max_Ep, Index] = max(Ep_vector)   ;
max_f = f_vector(Index)  ;


U_hole_res = real(Create_U_hole(N, M, max_f, theta_end)) ;


% Hole, resonant f
figure(7)
mesh(x_vec,y_vec,U_hole_res)
grid on
xlabel("x")
ylabel("y")
zlabel("u(x,y)")
title("u(x,y), hole,  at resonant frequency f = " + max_f)


% Close to Resonant
U_hole_res_close = real(Create_U_hole(N, M, max_f + 2, theta_end)) ;

figure(8)
mesh(x_vec,y_vec,U_hole_res_close)
grid on
xlabel("x")
ylabel("y")
zlabel("u(x,y)")
title("u(x,y), hole,  close to resonant frequency f = " + (max_f + 2))


[min_Ep, Index] = min(Ep_vector)   ;
min_f = f_vector(Index)  ;

% Far away from Resonant
U_hole_res_far = real(Create_U_hole(N, M, min_f + 2, theta_end)) ;

figure(9)
mesh(x_vec,y_vec,U_hole_res_far)
grid on
xlabel("x")
ylabel("y")
zlabel("u(x,y)")
title("u(x,y), hole, far away from resonant frequency f = " + (min_f + 2))

%% Sound Intensity level as a function of frequency
L_f_vec = [] ;
df = 30;


for f = 100:df:4000
    U_matrix = Create_U_hole(N, M, f, theta_end) ;
    L = Sound_Intensity(U_matrix, M, theta_end) ;
    L_f_vec = [L_f_vec; f, L] ;
end


figure(10)
plot(L_f_vec(:,1), L_f_vec(:,2))
yline(0)
grid on
xlabel("Frequency [Hz]")
ylabel("Sound Intensity Level [dB]")
title("Sound Intensity Level at 10 [m] from the cavity")

%% Functions

% Sound intensity 
function L = Sound_Intensity(U, M, theta_end)
    % Data
    H = 0.05 ; % [m]
    R1 = 0.5 ; % [m]
    c0 = 343 ; % [m/s]
    rho = 1.225 ; % [kg/m^3]
    
    h = 2*pi/M ;
    P = 0 ;
    U_R1 = U(end,:) ;
    
    theta_end_index = floor((M/(2*pi/theta_end))) ;
    
    U_R1 = U_R1(1:theta_end_index) ;
    for l = 1:theta_end_index
        P = P + h*((R1*H)/(2*rho*c0))*(abs(U_R1(l)))^2 ;
    end
    
    d = 10; % distance from the cavity [m]
    I0 = 10^-12 ;
    Ir = P/(4*pi*d^2) ;
    
    % Sound Intensity level
    L = 10*log10(Ir/I0) ; % [dB]
end


% Generating potential energy as function of frequency
function EP_f = frequency_sweep(N, M, theta_end)
    
    % Step Size
    df = 10;
    
    EP_f = [] ;
    for f = 1000:df:1300
        if theta_end == 0 
            U_matrix = Create_U(N, M, f) ;
        else
            U_matrix = Create_U_hole(N, M, f, theta_end);
        end
        Ep = Potential_Energy(U_matrix, N, M) ;
        
        EP_f = [EP_f; f, Ep] ;
    end
end

% Calculating potential energy for given frequency
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

% Solving system without hole
function U = Create_U(N, M, f)
    A = sparse(Create_A_Matrix(N, M, f)) ;
    g_vec = sparse(Create_g_vector(N, f)) ;
    
    U = sparse(A\g_vec) ;
    U  = sparse(reshape(U, N+1,M));
    U = [U,U(:,1)];
end

% Solving system with a hole
function U = Create_U_hole(N, M, f, theta_end)
    A = sparse(Create_A_Matrix_hole(N, M, f, theta_end)) ;
    g_vec = sparse(Create_g_vector(N, f)) ;
    
    U = sparse(A\g_vec) ;
    U  = sparse(reshape(U, N+1,M));
    U = [U,U(:,1)];
end

% Creating g vector
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
    
    g = @(r, theta) rho*c0*k*q*((1/((omega^2)*H*pi))*exp((-1)*((Rs/omega)^2)*(((1-r/Rs))^2 + (pi - theta)^2))) ;
    
    g_vec = sparse([]) ;
    
    % Creating discretized source function function
    for theta = 0:h:2*pi-h 
        for r = R0:h:R1 
            g_vec = [g_vec;g(r,theta)];
        end
    end  
end

% Creating A matrix with hole
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
    for index = 1:(M)
        if index == (M)
            A(end-N:end,1:N+1) = sparse(diag(d));
            A((index-1)*(N+1) + 1 :(index)*(N+1), (index-2)*(N+1) + 1 :(index-1)*(N+1)) = sparse(diag(e));
        elseif index == 1
            A(1:N+1,end-N:end) = sparse(diag(e));
            A((index-1)*(N+1) + 1 :(index)*(N+1), (index)*(N+1) + 1 :(index+1)*(N+1)) = sparse(diag(d));
        else
            A((index-1)*(N+1) + 1 :(index)*(N+1), (index-2)*(N+1) + 1 :(index-1)*(N+1)) = sparse(diag(e));
            A((index-1)*(N+1) + 1:(index)*(N+1), (index)*(N+1) + 1 :(index+1)*(N+1)) = sparse(diag(d));
        end
    end
end

% Creating A matrix with no hole
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
    A_Diagonal(1,2) = a(1) + c(1) ;
    A_Diagonal(end,end-1) = c(end) + a(end) ;
    %A_Diagonal(end, 1) = a(end) + c(end) ;
    %A_Diagonal(1, end) = c(1) + a(1);

    % Creating big Matrix A
    A = sparse(kron(sparse(eye(M,M)),A_Diagonal)) ;

    % Adding additional elements and BC:s
    for index = 1:(M)
        if index == (M)
            A(end-N:end,1:N+1) = sparse(diag(d));
            A((index-1)*(N+1) + 1 :(index)*(N+1), (index-2)*(N+1) + 1 :(index-1)*(N+1)) = sparse(diag(e));
        elseif index == 1
            A(1:N+1,end-N:end) = sparse(diag(e));
            A((index-1)*(N+1) + 1 :(index)*(N+1), (index)*(N+1) + 1 :(index+1)*(N+1)) = sparse(diag(d));
        else
            A((index-1)*(N+1) + 1 :(index)*(N+1), (index-2)*(N+1) + 1 :(index-1)*(N+1)) = sparse(diag(e));
            A((index-1)*(N+1) + 1:(index)*(N+1), (index)*(N+1) + 1 :(index+1)*(N+1)) = sparse(diag(d));
        end
    end
end