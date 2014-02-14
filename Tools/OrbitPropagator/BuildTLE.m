%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacob Cook
% 10/5/2013
% OE2TLE.m
%
% input:    SatNum = 5 ASCII character string Satellite number
%           Class =  1 ASCII character satellite class
%           IntDsgntr = 8 ASCII char. String [year(2) launch#(3) peice(3)]
%           Epoch = Satellite epoch [datenum]
%           OE = Classical orbital element structure
%           
%
% output:   TLE = Two Line Element 
%
% This function builds the Two Linw Element from the passed variables
% listed in the inputs.  
% ** As of right now the function is only concerned with the sencond line
% of the TLE and the Epoch time int he first TLE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TLE] = BuildTLE(SatNum, Class, IntDsgntr, Epoch, OE)

%%
global MU_EARTH

%% Get Epoch date and put in TLE format
MonthInDays = [0 31 59 90 120 151 181 212 243 273 304 334];
E_year = year(Epoch);
E_month = month(Epoch);
E_day = day(Epoch);
E_hour = hour(Epoch);
E_min = minute(Epoch);
E_sec = second(Epoch);

TLE_year = E_year-floor(E_year/100)*100;
TLE_day = MonthInDays(E_month) + E_day + (E_hour+E_min/60+E_sec/3600)/24;
TLE_Epoch = sprintf('%02d000.00000000',TLE_year);

if TLE_day >= 100
    TLE_Epoch(3:14) = sprintf('%0.8f',TLE_day);
elseif TLE_day >= 10 
    TLE_Epoch(4:14) = sprintf('%0.8f',TLE_day);
else
    TLE_Epoch(5:14) = sprintf('%0.8f',TLE_day);
end


%% Check Length of parameters
if length(SatNum) ~= 5
    disp('Error: Satellite number length needs to be 5 characters')
    return
end

if length(Class) ~= 1
    disp('Error: Satellite class needs to be a single character')
    return
end

if length(IntDsgntr) ~= 8
    disp('Error: International Designator should have length 8')
    return
end

%% TLE Header
TLE.header = 'Object J (updated 2012.09.24 12:00:00 from cubesat.org)';

%% TLE Line 1
% Base Line (We're not concerned with the pertubation stuff right now)
TLE.line1 = '1 90038U 0        12266.69397895 +.00000000 +00000-0 +00000-0 0 00168';
% insert passed paramters
TLE.line1(3:7) = SatNum;
TLE.line1(8) = Class;
TLE.line1(19:32) = TLE_Epoch;

%% TLE line 2
TLE.line2 = '2 90038 000.0000 000.0000 0000000 000.0000 000.0000 00.00000000001287';
TLE.line2(3:7) = SatNum;

% inclination [degrees]
TLE_i = sprintf('%.4f',OE.i*180/pi);
TLE.line2(17-length(TLE_i):16) = TLE_i;

% RAAN [degrees]
TLE_RAAN = sprintf('%.4f',OE.RAAN*180/pi);
TLE.line2(26-length(TLE_RAAN):25) = TLE_RAAN;

% eccentricity
e_str = sprintf('%.7f',OE.e);
TLE_e = e_str(3:9);
TLE.line2(27:33) = TLE_e;

% argument of perigee [degrees]
TLE_omega = sprintf('%.4f',OE.omega*180/pi);
TLE.line2(43-length(TLE_omega):42) = TLE_omega;

% Mean Anomoly [degrees]
E = nu2E(OE.nu, OE.e);
M = E2M(E, OE.e);
TLE_M = sprintf('%03.4f', M*180/pi);
TLE.line2(52-length(TLE_M):51) = TLE_M; 

% Mean Motion [revs/day]
T = 2*pi*(OE.a^3/MU_EARTH)^.5;
N = 24*3600/T;
TLE_N = sprintf('%.8f', N);
TLE.line2(64-length(TLE_N):63) = TLE_N;

% Rev number at Epoch [revs]
Revs = 0;
TLE.line2(64:68) = sprintf('%05d', Revs);
