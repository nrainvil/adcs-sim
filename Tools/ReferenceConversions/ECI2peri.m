function [Q] = ECI2peri(Orbit)

RAAN = Orbit.RAAN*pi/180;
omega = Orbit.ArgOfPerigee*pi/180;
i = Orbit.Inclination*pi/180;

c_i = cos(i);
s_i = sin(i);
c_RAAN = cos(RAAN);
s_RAAN = sin(RAAN);
c_omega = cos(omega);
s_omega = sin(omega);

Q = [-s_RAAN*c_i*s_omega+c_RAAN*c_omega  c_RAAN*c_i*s_omega+s_RAAN*c_omega s_i*s_omega;
     -s_RAAN*c_i*c_omega-c_RAAN*s_omega  c_RAAN*c_i*c_omega-s_RAAN*s_omega s_i*c_omega;
      s_RAAN*s_i                        -c_RAAN*s_i                        c_i];
  