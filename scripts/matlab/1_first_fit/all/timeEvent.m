% Define the timeEvent function
function [value,isterminal,direction] = timeEvent(~,~)
    value = toc - 1; % Stop if more than 1 seconds have passed
    isterminal = 1; % Terminate integration
    direction = 0; % Detect both increasing and decreasing events
end
