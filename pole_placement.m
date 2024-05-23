% Clearing workspace 
clear all; close all; clc;

% Loading the discretized linear system
System;
% The system is stored in the variable sysd
sysd

%% Compensator design via Pole Placement

% Step 1) defining the desired location of the poles 

% Open-loop poles
pole(sysd)

% Step 2) 

