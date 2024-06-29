% Load the robot model
urfive = loadrobot("universalUR5", DataFormat="row", Gravity=[0 0 -9.81]);

% Generate a random valid configuration
q = randomConfiguration(urfive);

% Calculate the gravity vector
G = gravityTorque(urfive, q);

% Display the gravity vector
% disp('Gravity Vector G:');
% disp(G);

% Generate random joint velocities
qd = rand(1, 6); % Random joint velocities

% Calculate the mass matrix at the given configuration
M = massMatrix(urfive, q);

% Calculate the gravity matrix at the given configuration
G = gravityTorque(urfive,q);

% Initialize Coriolis matrix
C = zeros(size(M));

% Compute the Coriolis matrix
n = length(q);
for k = 1:n
    for j = 1:n
        for i = 1:n
            % Compute the Christoffel symbols
            C(k, j) = C(k, j) + 0.5 * (partialDerivative(M, k, i, q, urfive) + ...
                partialDerivative(M, k, j, q, urfive) - partialDerivative(M, i, j, q, urfive)) * qd(i);
        end
    end
end

% Function to compute the partial derivatives of the mass matrix
function pd = partialDerivative(M, row, col, q, urfive)
    delta = 1e-5;
    q_temp = q;
    q_temp(col) = q_temp(col) + delta;
    M1 = massMatrix(urfive, q_temp);
    q_temp(col) = q_temp(col) - 2 * delta;
    M2 = massMatrix(urfive, q_temp);
    pd = (M1(row, col) - M2(row, col)) / (2 * delta);
end

% Display the Coriolis matrix
%disp('Coriolis Matrix C:');
%disp(C);

% Define joint accelerations (for example purposes)
qdd = rand(1, 6); % Random joint accelerations

% Calculate the torques
tau = M * qdd' + C * qd' + G';

% Display the torques
%disp('Joint Torques tau:');
%disp(tau);

% State-space representation
A = [zeros(n) eye(n); zeros(n) -M\C];
B = [zeros(n); M\eye(n)];
C_ss = eye(2*n);
D = zeros(2*n, n);

% Create the state-space model
sys = ss(A, B, C_ss, D);

% Display the state-space model
%disp('State-Space Model:');
%disp(sys);

% Define LQR weight matrices
Q = eye(2*n); % State weighting matrix
R = eye(n);   % Control effort weighting matrix

% Design the LQR controller
K = lqr(sys, Q, R);

% Closed-loop system
Ac = A - B * K;
Bc = B;
Cc = C_ss;
Dc = D;

sys_cl = ss(Ac, Bc, Cc, Dc);

% Initial state
x0 = [q, qd]'; % Initial condition for states

% Time vector for simulation
t = linspace(0, 10, 1000);

% Simulate the closed-loop response
[y, t, x] = initial(sys_cl, x0, t);

% Plot the results
figure;
subplot(2, 1, 1);
plot(t, x(:, 1:n)); % Plot joint positions
title('Joint Positions');
xlabel('Time (s)');
ylabel('Position (rad)');
legend(arrayfun(@(i) sprintf('q%d', i), 1:n, 'UniformOutput', false));

subplot(2, 1, 2);
plot(t, x(:, n+1:end)); % Plot joint velocities
title('Joint Velocities');
xlabel('Time (s)');
ylabel('Velocity (rad/s)');
legend(arrayfun(@(i) sprintf('qd%d', i), 1:n, 'UniformOutput', false));