%{
function uv = QDEopticdecross_sir(pop, dv, cr, m, n)
% Quantum-inspired crossover using superposition and entanglement

u = m;
v = n;
xmin = 150;
xmax = v - 150;
ymin = 150;
ymax = u - 150;
ir = randi(4);

% Define Hadamard gate transformation
H = (1/sqrt(2)) * [1 1; 1 -1];

% Define Pauli-X gate transformation
X = [0 1; 1 0];

for i = 1:4
    %if rand(1) <= cr || i == ir
    if rand(1) <= cr + 0.1 * rand
    
        % Apply quantum superposition via Hadamard gate
if size(dv, 2) == 1  % If dv is a column vector
    qbit = [pop(1, i); dv(i)];  
else
    qbit = [pop(1, i); dv(1, i)];
end
        qbit = H * qbit;  % Apply Hadamard transform

        % Apply Pauli-X transformation
        qbit = X * qbit;
        
        % Measure quantum state (collapse)
        if rand < 0.5
            uv(1, i) = qbit(1);
        else
            uv(1, i) = qbit(2);
        end
        
        % Ensure the values remain in valid range
        if i == 1
            uv(1, i) = min(max(uv(1, i), xmin), xmax);
        elseif i == 2
            uv(1, i) = min(max(uv(1, i), ymin), ymax);
        elseif i == 3 || i == 4
            uv(1, i) = min(max(uv(1, i), -10), 10);
        end
    else
        uv(1, i) = pop(1, i);
    end
end

end
%}


function uv = QDEopticdecross_sir(pop, dv, cr, m, n)
% Quantum-inspired crossover using probability-based selection


ir = randi(4); % Randomly force crossover for at least one element
%ir = randperm(4, randi([2, 3])); % Force crossover for 2 or 3 elements

% Define Hadamard and Pauli-X gates
H = (1/sqrt(2)) * [1 1; 1 -1]; % Hadamard gate
X = [0 1; 1 0]; % Pauli-X gate
theta = (0.1 + 0.4 * rand()) * pi; % Controlled randomness for rotation
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];


for i = 1:4
    if rand(1) <= cr || i == ir  % Perform crossover based on cr
    %if rand(1) <= cr || ismember(i, ir)
        % Create a quantum state using Hadamard (superposition)
        %{
        if size(dv, 2) == 1  % If dv is a column vector
            qbit = [pop(1, i); dv(i)];
        else
            qbit = [pop(1, i); dv(1, i)];
        end
        %}

        %qbit = [rand(); sqrt(1 - rand()^2)];

        qbit = [1; 0]; % Assume an initial state |0⟩

        qbit = H * qbit;  % Apply Hadamard transform
        qbit = X * qbit;  % Apply Pauli-X transformation
        %qbit=R * qbit;

        % Measure quantum state (collapse)
        alpha = abs(qbit(1))^2; % Probability of first state
        beta = abs(qbit(2))^2; % Probability of second state
        
        % Use probability to select from parent or mutant
        if rand(1) <= alpha
            uv(1, i) = dv(1, i); % Take from mutant (dv)
        else
            uv(1, i) = pop(1, i); % Take from parent (pop)
        end
    else
        uv(1, i) = pop(1, i); % No crossover, inherit from parent
    end
    
   
end

end

