function [M] = delta_t2M(delta_t, mu, a)

n = (mu/a^3)^.5;
M = n*delta_t;