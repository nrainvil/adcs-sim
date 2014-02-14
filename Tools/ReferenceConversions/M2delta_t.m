function [delta_t] = M2delta_t(M, mu, a)

n = (mu/a^3)^.5;
delta_t = M/n;