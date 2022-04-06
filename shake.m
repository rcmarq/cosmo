function shake(element,sample,eff,cycles,R,R_reference,T,deltat)

% SHAKE   Sets up the parameters necessary for performing calculations
% under the COSMO and ERRORCURVESPLIT functions
%
%  SHAKE(element,sample,eff,cycles,R,R_reference,T,deltaT)
%             element -- the element of interest. 
%                sample -- the amount of sample to be analyzed - in
%             nanograms(ng)
%             eff -- 2 x 1 array that contains the ionization efficiency of
%                the two instruments used for the IC and ID measurements,
%                respectively. Range from 0 to 1 
%                e.g. [0.01 0.005]
%             cycles -- 2 x 1 array that contains the no. of cycles
%                for the IC and ID measurements, respectively
%                e.g. [0.01 0.005]
%             R -- 2 x nisos array that contains the resistance values for
%                the IC (1st row) and ID (2nd row) measurements, respectively.
%                By default, all resistors are 10^11 ohms. 
%             R_reference -- reference resistance used for describing beam intensity.
%                The total beam current is intensity/R_reference
%                e.g. the 10 V default beam corresponds to 100 pA with R_reference=1e11 Ohms.
%             T -- absolute temperature of the resistor. Defualt values is
%                300K.  
%             deltat -- integration time for each measuremebnt cycle. In
%                seconds (s)
%      
%                
% All of these parameters (or resulting derived values are stored in the
% global ISODATA variable. 
%
% Example
%   shake('Sr',5.0,[0.01 0.005],[1e11 1e11 1e11 1e11; 1e11 1e11 1e11
%   1e11],1e11)

global ISODATA

% Fundamental constants
elementarycharge = 1.60217646e-19;  % Coulombs
k                = 1.3806504e-23;   % Boltzman's constant (m^2 kg s^-2 K^-1)
Avn              = 6.0221409e+14;   % Avogadro's number (nmol^-1)

nisos = ISODATA.(element).nisos; % number of isotopes for element of interest
errormodel = ISODATA.(element).errormodel; % errormodel as pre-defined by seterrormodel.m

if (nargin<8)
        deltat = 8; % intergration time; in seconds
end

if (nargin<7)
        T = 300; % temperature of amplifier/detector; in Kelvin
end

if (nargin<6)
        R_reference = 1e11; 
end

if (nargin<5)
        R = zeros(2,nisos)+1e11;
end

if (nargin<4)
        cycles = [40 40]; % fix this at some point - should be defined by R and saturation voltage
end

if (nargin<3)
		eff= [0.01 0.01]; % by default, use a different error model for these, as these requires two runs.
end

% Check that the resistor values are set for all isotopes
if size(R,2) ~= nisos
    disp(R)
    disp(size(R,2))
    disp(nisos)
    disp("ERROR: The resistor value array does not match the number of isotopes for the element of interest")
    return
end

% Set the amount of sample in nanograms
MM    = sum(ISODATA.(element).standard.*ISODATA.(element).mass); % molar mass of element
atoms = Avn*sample/MM; % number of atoms introduced

errormodel.V100  = atoms*elementarycharge*R_reference/(deltat); % voltage at ref. resistance at 100% ionization efficiency

% NOTE: The 'intensity' parameter in the 'measured' and 'standard'
% substructures are not yet set here, as it will be determined by the
% splitting of the sample - which is a variable set during the calculation
% of uncertainties by ERRORESTIMATESPLIT.m

% Set all error models to be of type 'fixed-total' 
% and add 'cycles','effIC', and 'effID' parameter
errormodel.measured.type = 'fixed-sample';
errormodel.spike.type    = 'fixed-sample';
errormodel.standard.type = 'fixed-sample';

errormodel.standard.eff  = eff(1); % efficiency of the IC measurement
errormodel.measured.eff  = eff(2); % efficiency of the ID measurement

errormodel.standard.cycles  = cycles(1); % no. of cycles of the IC measurement
errormodel.measured.cycles  = cycles(2); % no. of cycles of the ID measurement

aIC = 4*k*T*R_reference^2./(deltat*R(1,:)); % Johnson-Noise for IC Measurement
aID = 4*k*T*R_reference^2./(deltat*R(2,:)); % Johnson-Noise for ID Measurement

b   = elementarycharge*R_reference/(deltat).*ones(1,nisos); % Counting statistics
c   = zeros(1,nisos);

errormodel.standard.a=aIC;
errormodel.standard.b=b;
errormodel.standard.c=c;

errormodel.measured.a=aID;
errormodel.measured.b=b;
errormodel.measured.c=c;

ISODATA.(element).errormodel = errormodel;

end
