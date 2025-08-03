
%{
function uv = QDElocalsearch(uv, m, n)
    perturb = [-10, -10, 10, 5]; % Small changes in all dimensions
    uv = uv + perturb .* randn(size(uv)); % Gaussian perturbation
    uv(:,1) = max(min(uv(:,1), n-150), 150);
    uv(:,2) = max(min(uv(:,2), m-150), 150);
end
%}
%{
function uv = QDElocalsearch(uv, m, n)
    % Small rotation angle (in radians)
    theta = 0.1; % Adjust this value as needed

    % Rotation matrix for 2D (since we are working with x and y coordinates)
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

    % Apply rotation to the x and y coordinates
    uv(1:2) = (R * uv(1:2)')';

    % Add Gaussian perturbation to all dimensions
    perturb = [-10, -10, 10, 5]; % Small changes in all dimensions
    uv = uv + perturb .* randn(size(uv)); % Gaussian perturbation

    % Ensure the values remain within valid bounds
    uv(:,1) = max(min(uv(:,1), n-150), 150);
    uv(:,2) = max(min(uv(:,2), m-150), 150);
end
%}

function uv = QDElocalsearch(uv, m, n, gen, ngen)
    % Adaptive Perturbation Scale for Faster Convergence
    base_perturb = 3; % Initial perturbation
    perturb_scale = base_perturb * (1 - gen/ngen); % Decrease perturbation as generations progress

    % Apply Gaussian Perturbation
    perturb = perturb_scale * randn(size(uv)); 
    uv = uv + perturb;  

    % Ensure the values remain within valid bounds
    uv(1) = max(min(uv(1), n-150), 150);
    uv(2) = max(min(uv(2), m-150), 150);
    uv(3) = max(min(uv(3), 10), -10);
    uv(4) = max(min(uv(4), 10), -10);
end


