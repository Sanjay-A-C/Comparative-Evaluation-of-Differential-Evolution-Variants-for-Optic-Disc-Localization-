%{
function [B, q_population] = QDEsolset_sir(npop, u, v)
    rng(sum(100 * clock), 'twister');  % Random seed initialization
    xmin = 150;
    xmax = v - 150;
    ymin = 150;
    ymax = u - 150;

    % Quantum-Inspired Population Initialization
    B = zeros(npop, 4);  
    q_population = rand(npop, 4) * 2 * pi;  % Quantum angles (theta)

    % Convert Qubit angles to classical values
    for i = 1:npop
        B(i,1) = xmin + (xmax - xmin) * abs(sin(q_population(i,1))); % Ensure within [xmin, xmax]
        B(i,2) = ymin + (ymax - ymin) * abs(sin(q_population(i,2))); % Ensure within [ymin, ymax]
        B(i,3) = -10 + 20 * abs(sin(q_population(i,3))); % Ensure within [-10,10]
        B(i,4) = -10 + 20 * abs(sin(q_population(i,4))); % Ensure within [-10,10]
    end

    % Ensure (x, y) stays strictly within bounds
    B(:,1) = max(min(B(:,1), xmax), xmin);
    B(:,2) = max(min(B(:,2), ymax), ymin);

    % Ensure r1 and r2 stay within [-10, 10] (should already be, but extra check)
    B(:,3) = max(min(B(:,3), 10), -10);
    B(:,4) = max(min(B(:,4), 10), -10);
end
%}

function [B, q_population] = QDEsolset_sir(npop, u, v)
    rng(sum(100 * clock), 'twister');  % Random seed initialization
    xmin = 150;
    xmax = v - 150;
    ymin = 150;
    ymax = u - 150;

    % Quantum-Inspired Population Initialization
    B = zeros(npop, 4);  
    q_population = rand(npop, 4) * 2 * pi;  % Quantum angles (theta)

    % Convert Qubit angles to classical values (and round to integer)
    for i = 1:npop
        B(i,1) = round(xmin + (xmax - xmin) * (0.5 + 0.5 * sin(q_population(i,1)))); % Integer x
        B(i,2) = round(ymin + (ymax - ymin) * (0.5 + 0.5 * sin(q_population(i,2)))); % Integer y
        B(i,3) = round(-10 + 20 * (0.5 + 0.5 * sin(q_population(i,3)))); % Continuous r1
        B(i,4) = round(-10 + 20 * (0.5 + 0.5 * sin(q_population(i,4)))); % Continuous r2
    end

    % Final safety check to strictly enforce boundaries
    B(:,1) = max(min(B(:,1), xmax), xmin);
    B(:,2) = max(min(B(:,2), ymax), ymin);
    B(:,3) = max(min(B(:,3), 10), -10);
    B(:,4) = max(min(B(:,4), 10), -10);
end
