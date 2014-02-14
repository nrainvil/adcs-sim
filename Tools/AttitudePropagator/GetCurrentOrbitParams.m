function orbitParams = GetCurrentOrbitParams(Orbit, idx)

orbitParams.posECEF_km  = Orbit.posECEF_km(:,idx);
orbitParams.posECI_km   = Orbit.posECI_km(:,idx);
orbitParams.velECI_kmps = Orbit.velECI_kmps(:,idx);
orbitParams.sun_ECI     = Orbit.Sun_ECI(:,idx);
orbitParams.B_ECI       = Orbit.B_ECI(:,idx);
orbitParams.time        = Orbit.Time(idx);
orbitParams.delta_t     = Orbit.Time(idx+1) - Orbit.Time(idx);