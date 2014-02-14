function [posECI_km posECEF_km velECI_kmps S_ECI B_ECI]=GenPosVelSunMag(time,tle)
%This function generates spacecraft position, velocity, sun & magnetic field vectors
% INPUTS:
%   time: [1xN double] (UTC) time in matlab datenum format
%   tle: [1x1 structure] contains TLE information. Has fields:
%       header: [char] indentifying info for the TLE. Not used in this function.
%       line1: [1x69 char] line 1 of the TLE
%       line2: [1x69 char] line 2 of the TLE
% OUTPUTS:
%   posECI_km: [3xN double] (km) satellite position in the ECI frame
%   posECEF_km: [3xN double] (km) satellite position in the ECEF frame
%   velECI_km: [3xN double] (km/s) satellite velocity in the ECI frame
%   S_ECI: [3xN double] (unitless) unit vector from earth to sun in ECI frame 
%       (in LEO, this vector is within 0.002° of the satellite to sun vector,
%       smaller than the 0.01° accuracy of the empircal sun position model)
%   B_ECI: [3xN double] (Tesla) earth magnetic flux density at the satellite
%       in the ECI frame, as calculed by the IGRF
%
%   Code written by David Gerhardt
%   Last update 2013.06.03

%% User Input
    if nargin==0
        %(if no user input given to function, use below instead)
        begDateNum = datenum('03 Oct 2012 01:00:00'); %[UTC]
        endDateNum = datenum('03 Oct 2012 02:37:00'); %[UTC]
        timeStep = 1; %[sec]
        time = begDateNum:timeStep/(3600*24):endDateNum;
        tle.header = 'Object J (updated 2012.09.24 12:00:00 from cubesat.org)';
        tle.line1 = '1 90038U 0        12266.69397895 +.00001805 +00000-0 +17857-3 0 00168';
        tle.line2 = '2 90038 064.6725 013.7141 0218862 286.5894 188.9890 14.79123198001287';
    end

%% Determine satellite position vector for each time step
    [posECI_km posECEF_km velECI_kmps S_ECI GST]=OrbitProp(time,tle);
    
%% Use IGRF to Determine ECEF H-field vector vs. time
    B_ECEF=IGRF(posECEF_km); %[Tesla]

%% Convert ECEF to ECI H-field vector
    B_ECI = zeros(3,size(B_ECEF,2));
    for i=1:size(posECEF_km,2)
        R_E2I = [cos(GST(i)) -sin(GST(i)) 0; sin(GST(i)) cos(GST(i)) 0; 0 0 1]; %rotation matrix: ECI to ECEF
        B_ECI(:,i)=R_E2I*B_ECEF(:,i); %[Tesla]
    end
    
%% Determine Vector Magnitudes
    B_mag = sqrt(B_ECI(1,:).^2+B_ECI(2,:).^2+B_ECI(3,:).^2); %[Tesla]
    S_mag = sqrt(S_ECI(1,:).^2+S_ECI(2,:).^2+S_ECI(3,:).^2); %[unitless] 1 in sun, 0 in umbra/penumbra
    pos_mag = sqrt(posECI_km(1,:).^2+posECI_km(2,:).^2+posECI_km(3,:).^2); %[km]
    vel_mag = sqrt(velECI_kmps(1,:).^2+velECI_kmps(2,:).^2+velECI_kmps(3,:).^2); %[km/s]
    
%% BONUS: plots (for sanity checks)
%     figure(1); clf;
%     plot(time,B_mag*1e6);
%     datetick('x','HH:MM','keeplimits');
%     xlabel(['UTC time on ',datestr(begDateNum,'yyyy.mm.dd')]);
%     ylabel('B-field magnitude [\muTesla]');
%     
%     figure(2); clf;
%     plot(time,S_mag); ylim([-0.1 1.1]);
%     datetick('x','HH:MM','keeplimits');
%     set(gca,'YTick',[0 1],'YTickLabel',{'Eclipse','Insolation'});
%     xlabel(['UTC time on ',datestr(begDateNum,'yyyy.mm.dd')]);
%     
%     figure(3); clf;
%     plot(pos_mag-6371,vel_mag);
%     xlabel('Spacecraft altitude [km]');
%     ylabel('Spacecraft velocity magnitude [km/s]');
%     
end
function [posECI_km posECEF_km velECI_kmps posEarth2SunECI_unit GST]=OrbitProp(time,tle)
%% Constants
    % distance from earth to sun
      AU = 149597871; %[km]
    % average sun radius
      Rs_avg = 696342; %[km]
    % average earth radius 
      Re_avg = 6371; %[km]

%% Define directories where extra matlab functions are found
%
% Commented out by Jacob Cook this is not needed when working in my tool
% set
% 
%     if ~isdeployed
%         addpath([pwd,'/coordinate_transformations']);
%         addpath([pwd,'/TimeConversions']);
%         addpath([pwd,'/vallado_tle']);
%     end

%% Determine Satellite Position
    %TLE epoch
    dateNumEpYear = datenum(['0 Jan 20',tle.line1(19:20)]);
    dateNumEpoch = dateNumEpYear + str2double(tle.line1(21:32));
    
    t_min = (time-dateNumEpoch)*24*60; % propagator uses time in minutes since tle epoch
    % Propagate tle
    [posECEF_km posECI_km velECI_kmps] = propagate_tle(tle,t_min);
    
    %determine Greenwich Sidereal Time Vector (to convert ECI <--> ECEF)
    %(based on http://www.astro.umd.edu/~jph/GST_eqn.pdf)
        %GST at start of 2011
        GST_2011=6.6208844; %[hours]
        %day fractions from start of 2012 to each time
        dR2011 = time-datenum([2011 0 0 0 0 0]);
        %whole days from start of 2012 to each time
        dayR2011 = floor(dR2011);
        %hours from start of 2012 to each time
        hrR2011= mod(mod(dR2011,1)*24,24);
    GST = mod(GST_2011+0.0657098244*dayR2011+1.00273791*hrR2011,24)*pi/12; %[rad]
    
%% Compute sun direction in ECI frame
    %convert time vector to julian date
    tJ = juliandate(time);
    %calculate distance vector to sun 
    %  from http://www.dept.aoe.vt.edu/~cdhall/courses/aoe4140/attde.pdf
    T_UT1 = (tJ - 2451545)/36525;
    lambda_M_sun = 280.4606184+36000.77005361*T_UT1;
    M_sun = 357.5277233+35999.05034*T_UT1;
    lambda = lambda_M_sun + 1.914666471*sind(M_sun)+0.918994643*sind(2*M_sun);
    epsilon=23.439291-0.0130042*T_UT1;
    %calculate unit vector from earth (or satellite) to sun
    posEarth2SunECI_unit = [cosd(lambda);
                            sind(lambda).*cosd(epsilon);
                            sind(lambda).*sind(epsilon)];
    %calculate magnitude of earth to sun vector
    magS = AU*(1.000140612-0.016708617*cosd(M_sun)-0.000139589*cosd(2*M_sun));

%% Determine Umbra/Penumbra Eclipse Times
%from http://celestrak.com/columns/v03n01/
    %position vector from satellite to sun
    posEarth2SunECI_km = [magS.*cosd(lambda);
                          magS.*sind(lambda).*cosd(epsilon);
                          magS.*sind(lambda).*sind(epsilon)];
    posSat2SunECI_km = posEarth2SunECI_km+posECI_km; %[km]
    %position vector from satellite to earth
    posSat2EarthECI_km = -posECI_km;
    %calculate distance magnitude of vectors for each timestep
    magE = sqrt(posSat2EarthECI_km(1,:).^2+posSat2EarthECI_km(2,:).^2+posSat2EarthECI_km(3,:).^2); %[km]
    %calculate semidiameters of earth & sun
    thetaE = asind(Re_avg./magE);
    thetaS = asind(Rs_avg./magS);
    %calculate the angle between the center of earth & the sun
    theta = acosd((posSat2EarthECI_km(1,:).*posSat2SunECI_km(1,:)+...
                   posSat2EarthECI_km(2,:).*posSat2SunECI_km(2,:)+...
                   posSat2EarthECI_km(3,:).*posSat2SunECI_km(3,:))./(magE.*magS));
    %calculate the inds which are in an umbral eclipse
    umbraInds = find(theta<thetaE-thetaS);
    %calculate the inds which are in penumbral eclipse
    penumbraInds = intersect(find(abs(thetaE-thetaS)<theta),find(theta<thetaE+thetaS));
    
%% Set sun vector during eclipse to zero
    posEarth2SunECI_unit(:,umbraInds)=0; %#ok<FNDSB>
    posEarth2SunECI_unit(:,penumbraInds)=0;

%     figure(1);
%     plot(tV,theta,tV,thetaE-thetaS,tV(umbraInds),theta(umbraInds),'r.',tV(penumbraInds),theta(penumbraInds),'y.');
%     ylabel('Angle between center of earth & sun [degrees]');
%     datetick('x','keeplimits');
    
end
function B_ECEF=IGRF(posECEF_km)
%INPUTS
%   posECEF_km - [km] satellite position vector in ECEF frame, cartesian
%                coordinates

    %% Constants
    % Mean radius for IGRF (6371.2 km)
    R_mean = 6371.2;

    %% IGRF stuff
    % Load IGRF coefficients
    [G,H] = LoadCoIGRF(2012);
    % max degree of geopotential
    nmax = 10;
    % max order  of geopotential
    mmax = 10;
    % call function to compute schmidt coefficients
    Kschmidt = schmidt(nmax,mmax);
    %define output vector
    B_ECEF = zeros(3,size(posECEF_km,2));
    for i=1:size(posECEF_km,2)
        % call function to compute legendre polynomials
        [A,ctilde,stilde] = recursion(posECEF_km(:,i),nmax,mmax);
        % determine b-field B_E
        bepe = bfield(posECEF_km(:,i),nmax,mmax,Kschmidt,A,ctilde,stilde,G,H,R_mean); %[Tesla]
        B_ECEF(:,i) = bepe; %[Tesla] ECEF frame
    end

end
function [G,H] = LoadCoIGRF(year)
% function [G,H] = igrf_coeffs_gen(year);
% Enter year between 1900 - 2005 in multiple of 5 (e.g. 1915)
% From 2005-2010, enter year in multiple of 1 (e.g. 2006)

% Edits added to by Jacob Cook in order to be able to call the function
% from a different file loation of the GenPosVelSun.m file.
    filename = mfilename('fullpath'); % added by Jacob
    [path] = fileparts(filename); % added by Jacob
    load([path,'/IGRF_Coefficients/g2011.txt']); % changed from pwd to path
    g = g2011; clear g2011
    load([path,'/IGRF_Coefficients/h2011.txt']); % changed from pwd to path
    h = h2011; clear h2011

    G = zeros(14,14);
    H = zeros(14,14);

    year_to_column = [  1900	3
                        1905	4
                        1910	5
                        1915	6
                        1920	7
                        1925	8
                        1930	9
                        1935	10
                        1940	11
                        1945	12
                        1950	13
                        1955	14
                        1960	15
                        1965	16
                        1970	17
                        1975	18
                        1980	19
                        1985	20
                        1990	21
                        1995	22
                        2000	23
                        2005    24
                        2010    25];

    if year > 2010
        column = 25;
    else
        indx = find(year_to_column >= year);
        column = year_to_column(indx(1),2);
    end

    svg = zeros(14,14); % Secular variation past year 2010
    svh = zeros(14,14); % Secular variation past year 2010

    if year > 2010
        column = 25;
        for i = 1:length(g)
            svg(g(i,1)+1,g(i,2)+1) = (year-2010)*g(i,26)/1e9;
        end
        for i = 1:length(h)
            svh(h(i,1)+1,h(i,2)+1) = (year-2010)*h(i,26)/1e9;
        end
    end

    for i = 1:length(g)
        G(g(i,1)+1,g(i,2)+1) = svg(g(i,1)+1,g(i,2)+1) + g(i,column)/1e9;
    end

    for i = 1:length(h)
        H(h(i,1)+1,h(i,2)+1) = svh(h(i,1)+1,h(i,2)+1) + h(i,column)/1e9;
    end
end
function K = schmidt(nmax,mmax)

%+---------------------------------------------------------------------+
%
%     Purpose:
%
%     Compute coefficients that relate Schmidt functions to associated
%     Legendre functions.
%
%+---------------------------------------------------------------------+
%
%     Argument definitions:
%
%     nmax              Maximum degree of contributing spherical harmonics
%
%     mmax              Maximum order of contributing spherical harmonics
%
%     K		            coefficients that relate Schmidt functions to
%	                    associated Legendre functions (Ref. [1]).
%+---------------------------------------------------------------------+
%
% The number 1 is added to degree and order since MATLAB can't have an array
% index of 0.

    % Seed for recursion formulae
    K(2,2) = 1;

    % Recursion formulae
    for n = 1:nmax
        i=n+1;
        for m = 0:n
            j=m+1;
            if m == 0
                % Eq. (3), Ref. [2]
                K(i,j) = 1;
            elseif ((m >= 1) && (n >= (m+1)))
                % Eq. (4), Ref. [2]
                K(i,j) = sqrt((n-m)/(n+m))*K(i-1,j);
            elseif ((m >= 2) && (n >= m))
                % Eq. (5), Ref. [2]
                K(i,j) = K(i,j-1)/sqrt((n+m)*(n-m+1));
            end
        end
    end
end
function [A,ctilde,stilde] = recursion(repe,nmax,mmax)
%+---------------------------------------------------------------------+
%
%     Purpose:
%
%     Recursive calculations of derived Legendre polynomials and other
%     quantities needed for gravitational and magnetic fields.
%
%+---------------------------------------------------------------------+

% The number 1 is added to degree and order since MATLAB can't have an
% array index of 0.

    clear A;
    A=zeros(nmax+3,nmax+3);         % A(n,m) = 0, for m > n

    R_m = sqrt(repe'*repe);
    rhat = repe/R_m;

    u = rhat(3);                    % sin of latitude

    A(1,1)=1;                       % "derived" Legendre polynomials
    A(2,1)=u;
    A(2,2)=1;
    clear ctilde
    clear stilde
    ctilde(1) = 1; ctilde(2) = rhat(1);
    stilde(1) = 0; stilde(2) = rhat(2);

    for n = 2:nmax
        i=n+1;

        % Calculate derived Legendre polynomials and "tilde" letters
        % required for gravitational and magnetic fields.

        % Eq. (4a), Ref. [2]
        A(i,i) = prod(1:2:(2*n - 1));

        % Eq. (4b), Ref. [2]
        A(i,(i-1))= u*A(i,i);

        if n <= mmax
            %   p. 9,     Ref. [1]
            ctilde(i)  = ctilde(2) * ctilde(i-1) - stilde(2) * stilde(i-1);
            stilde(i)  = stilde(2) * ctilde(i-1) + ctilde(2) * stilde(i-1);
        end

        for m = 0:n
            j=m+1;
            if (m < (n-1)) && (m <= (mmax+1))
                %     Eq. I, Table 1, Ref. [2]
                A(i,j)=((2*n - 1)*u*A((i-1),j) - (n+m-1)*A((i-2),j))/(n-m);
            end
        end
    end
end
function bepe = bfield(repe,nmax,mmax,K,A,ctilde,stilde,G,H,R_mean)
%+---------------------------------------------------------------------+
%
%     Purpose:
%
%     Compute magnetic field exerted at a point P.
%
%+---------------------------------------------------------------------+
%
%     Argument definitions:
%
%     repe    (km)      Position vector from Earth's center, E*, to a
%                       point, P, expressed in a basis fixed in the
%                       Earth (ECF): 1 and 2 lie in equatorial plane
%                       with 1 in the plane containing the prime
%                       meridian, in the direction of the north pole.
%
%     nmax              Maximum degree of contributing spherical harmonics
%
%     mmax              Maximum order of contributing spherical harmonics
%
%     K		            coefficients that relate Schmidt functions to
%						associated Legendre functions.
%
%     A                 Derived Legendre polynomials
%
%     ctilde            See pp. 4--9 of Ref. [1]
%
%     stilde            See pp. 4--9 of Ref. [1]
%
%     G, H     Tesla    Schmidt-normalized Gauss coefficients
%
%     R_mean   km       Mean radius for International Geomagnetic
%                       Reference Field (6371.2 km)
%
%     bepe     Tesla    Magnetic field at a point, P, expressed in ECF
%                       basis
%
%+---------------------------------------------------------------------+
%
%     Conversion factors:
%
%       1 Tesla = 1 Weber/(meter-meter) = 1 Newton/(Ampere-meter)
%               = 1e+4 Gauss  =  1e+9 gamma
%
%+=====================================================================+

% The number 1 is added to degree and order since MATLAB can't have an array
% index of 0.

    e1=[1; 0; 0];
    e2=[0; 1; 0];
    e3=[0; 0; 1];

    rmag = sqrt(repe'*repe);
    rhat = repe/rmag;

    u = rhat(3);	% sin of latitude

    bepe = [0; 0; 0];

    % Seed for recursion formulae
    scalar = R_mean*R_mean/(rmag*rmag);
    for n = 1:nmax
        % Recursion formula
        scalar = scalar*R_mean/rmag;
        i=n+1;
        for m = 0:n
            j=m+1;
            if m <= mmax
                ttilde = G(i,j)*ctilde(j) + H(i,j)*stilde(j);
                %     ECF 3 component {Eq. (2), Ref. [2]}
                b3 = -ttilde*A(i,j+1);
                %     rhat component {Eq. (2), Ref. [2]}
                br = ttilde*(u*A(i,j+1) + (n+m+1)*A(i,j));
                %     Contribution of zonal harmonic of degree n to magnetic
                %     field.  {Eq. (2), Ref. [2]}
                bepe = bepe + scalar*K(i,j)*(b3*e3 + br*rhat);
            end
            if ((m > 0) && (m <= mmax))
                %     ECF 1 component {Eq. (2), Ref. [2]}
                b1 = -m*A(i,j)*(G(i,j)*ctilde(j-1) + H(i,j)*stilde(j-1));
                %     ECF 2 component {Eq. (2), Ref. [2]}
                b2 = -m*A(i,j)*(H(i,j)*ctilde(j-1) - G(i,j)*stilde(j-1));
                %     Contribution of tesseral harmonic of degree n and order m to
                %     magnetic field.  {Eq. (2), Ref. [2]}
                bepe = bepe + scalar*K(i,j)*(b1*e1 + b2*e2);
            end
        end
    end
end
function [x_ecf_km x_eci_km v_eci_kmps] = propagate_tle(tle,t_min)
    %% Inputs
    % tle: structure containing strings for each line of the tle
    %      tle.line1, tle.line2
    % t_min: time vector in minutes from epoch

    %% Outputs
    % x_ecf_km: satellite position in ECF (km)
    % x_eci_km: satellite position in ECI (km)

    satrec = twoline2rvMOD(tle.line1,tle.line2);
    jd0 = satrec.jdsatepoch;

    x_ecf_km = zeros(length(t_min),3);
    x_eci_km = zeros(length(t_min),3);
    v_eci_kmps = zeros(length(t_min),3);

    for i = 1:length(t_min)
    %     [satrec, x_ecf_km(i,:)] = spg4_ecf(satrec,t_min(i));
        jd(i) = jd0 + t_min(i)/60/24;
    %     x_eci_km(i,:) = ecf2eci(x_ecf_km(i,:),jd(i));
        [satrec, x_eci_km(i,:), v_eci_kmps(i,:)] = sgp4(satrec,t_min(i));
        x_ecf_km(i,:) = eci2ecf(x_eci_km(i,:),jd(i));
    end
    
    x_ecf_km = x_ecf_km';
    x_eci_km = x_eci_km';
    v_eci_kmps = v_eci_kmps';
    
end