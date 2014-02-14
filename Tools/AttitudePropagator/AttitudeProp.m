%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacob Cook
% 6/30/2013
% AttitudeProp.m
%
% input:    y0          = 7x1 state vector initial conditions
%           cubeSat     = cubeSat struct containing geometric info
%           Orbit       = Struct containing current orbit pos, vel, sun, B
%                           and time
%
% output:   y           = 7xn Solution to the spacecraft dynamic equations
%           N_distArray = 15xn Array of disturance torques 
%                         [N_aero; N_gravity; N_solar; N_magnet; N_dist]
%
% This function implements a 4th order Runge-Kutta propagator to integrate
% the spacecraft attitude in the orbit determined by the passed orbit 
% structure.  The solution depends on the initial condition y0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y, Outputs] = AttitudeProp(y0, cubeSat, Orbit)

%% book keeping
nSteps = length(Orbit.Time);

% Set initial condition and allocate space for the solution
y_n = y0;
y = zeros(length(y0(:,1)), nSteps);

%% Runge Kutta solver
for idx = 1:nSteps-1
    
    % get current set of orbit parameters
    CurrentOrbitParams = GetCurrentOrbitParams(Orbit, idx);
    
    % save current state
     y(:,idx) = y_n;
    
    % get current time and time step from current set of orbital params
    t_n = CurrentOrbitParams.time;
    h = CurrentOrbitParams.delta_t;
    
    % RK4
    [k1, outputs]= f_spacecraft(t_n, y_n, cubeSat, CurrentOrbitParams);
    k2 = f_spacecraft(t_n+0.5*h, y_n+0.5*h*k1, cubeSat, CurrentOrbitParams);
    k3 = f_spacecraft(t_n+0.5*h, y_n+0.5*h*k2, cubeSat, CurrentOrbitParams);
    k4 = f_spacecraft(t_n+h, y_n+h*k3, cubeSat, CurrentOrbitParams);
    
    % calculate next state
    y_np1 = y_n + h/6*(k1 + 2*k2 + 2*k3 + k4);
    y_n = y_np1;
   
    
    % handle outputs
    if idx == 1
        fnames = fieldnames(outputs);
        for i=1:numel(fnames)
            out = outputs.(fnames{i});
            Outputs.(fnames{i}) = zeros(numel(out), nSteps);
        end
    end
    
    if numel(fnames) > 0
        for i=1:numel(fnames)
            Outputs.(fnames{i})(:,idx) = outputs.(fnames{i});
        end
    end
                 
end

% stuff last one
y(:,end) = y_n;
% smudge this in 
if numel(fnames) > 0
    for i=1:numel(fnames)
        Outputs.(fnames{i})(:,end) = outputs.(fnames{i});
    end
end