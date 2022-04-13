function [error,ppmperamu]=errorwsplit(element,split,prop,spike,isoinv,INisos,errorratio,alpha,beta)
%ERRORWSPLIT Calculates the error in the natural fractionation factor or a chosen ratio by linear error propagation
%    ERRORWSPLIT(element,split,prop,spike,isoinv,INisos,errorratio,alpha,beta)
%             element -- element used in double spike, e.g. 'Fe'
%             split -- the splitting of the sample, i.e., the fraction of
%                total sample going to the double spike measurement.
%             prop -- proportion of double spike in double spike-sample mix.
%             spike -- the isotopic composition of the spike e.g. [0 0.5 0 0.5]
%                corresponds to a 50-50 mixture of the 2nd and 4th isotopes
%                (56Fe and 58Fe) in the case of Fe.
%             isoinv -- the isotopes used in the inversion, e.g. [54 56 57 58].
%                By default the first 4 isotopes are used.
%             INisos -- theisotopes used for the internal normalization,
%                assumes [n d]. 
%             errorratio -- by default, the error on the natural fractionation
%                factor (known as alpha) is given. Instead, the error on a 
%                particular ratio can be given by setting errorratio. e.g.
%                setting errorratio=[58 56] will give the error on 58Fe/56Fe.
%             alpha, beta -- there is a small dependance of the error on the fractionation
%                factors (instrumental and natural, or alpha and beta). Values of alpha and
%                beta can be set here if desired, although the effect on the optimal spikes
%                is slight unless the fractionations are very large. Default is zero.
%
% Output: error -- the error on the fractionation factor, or the specified ratio. 
%         ppmperamu -- the error converted to an approximate ppm per atomic mass unit
%
% Note that a number of parameters are specified in the global variable ISODATA.
%
% Example
%   error=errorestimate('Fe',0.5,[0 0.5 0 0.5])
%
% See also dsstartup
global ISODATA

if isempty(ISODATA)
	dsstartup;
end
if (nargin<9) || isempty(beta)
	beta=0;
end
if (nargin<8) || isempty(alpha)
	alpha=0;
end
if (nargin<7) || isempty(errorratio)
	errorratio=[];
end
if (nargin<6)
	INisos=[];
end
if (nargin<5) || isempty(isoinv)
	isoinv=[1 2 3 4];
end

rawdata=ISODATA.(element);
spike=spike/sum(spike);

% Splitting the sample and updating the errormodel in the ISODATA structure
cyclesIC=rawdata.errormodel.standard.cycles;
cyclesID=rawdata.errormodel.measured.cycles;

sampleIC=rawdata.errormodel.V100*(1-split)*rawdata.errormodel.standard.eff/cyclesIC; % total voltage of sample per cycle
sampleID=rawdata.errormodel.V100*split*rawdata.errormodel.measured.eff/cyclesID;     % total voltage of sample per cycle

rawdata.errormodel.standard.intensity=sampleIC;
rawdata.errormodel.measured.intensity=sampleID;

% Convert isotope mass numbers to index numbers
errorratio=rawdata.isoindex(errorratio);
isoinv=rawdata.isoindex(isoinv);
INisos=rawdata.isoindex(INisos);

% Check that the denominator for the internal normalization is part of
% the DS inversion isotopes

if ~isempty(INisos) && ~ismember(INisos(2),isoinv)
    disp('Error: the denominator for the internal normalization must be one of the inversion isotopes')
    return
end

if isempty(INisos)
    [m,ix]=max(spike(isoinv));  % find largest spike of those used in inversion
    d=isoinv(ix);
else
    d=rawdata.isoindex(INisos(2)); % the denominator isotope is the one spec. by the user
end

n=isoinv(isoinv~=d);
isoinv=[d n];              % rearrange isoinv so we denominator with largest isotope 

% Calculate ratios
in=calcratioeddata(element,isoinv);       % Calculates ratios for standard, ratio of masses etc
in.AT=spike(in.Ani)./spike(in.di);        % Calculate ratio of spike

% Now calculate sample ratio, lambda etc
in.AN=in.An.*exp(-in.AP*alpha);
lambda=realproptoratioprop([prop 1-prop],[in.AT; in.AN]);
lambda=lambda(1);
z=[lambda alpha beta];
in.AM=lambda.*in.AT +(1-lambda).*in.AN;
in.Am=in.AM.*exp(in.AP*beta);

% Error propagation
measured=ones(1,rawdata.nisos);
measured(in.Ani)=in.Am;
measured=measured/sum(measured);

%isonorm=isoinv;  % normalise so that only the sum of beams used in the inversion is the mean intensity
isonorm=1:rawdata.nisos; % normalise so that the sum of all beams is the mean intensity
if isempty(INisos)
    in.VAn=calcratiocov(rawdata.standard,rawdata.errormodel.standard,in.di,isonorm,0)./cyclesIC; % prop set at 0 since no spike at N meas.
else
    in.VAn=calcratiocovIN(element,rawdata.standard,rawdata.errormodel.standard,INisos)./cyclesIC;
end
in.VAT=calcratiocov(spike,rawdata.errormodel.spike,in.di,isonorm,prop);
in.VAm=calcratiocov(measured,rawdata.errormodel.measured,in.di,isonorm,prop)./cyclesID;
[Vz, VAN]=fcerrorpropagation(z,in.AP,in.An,in.AT,in.Am,in.VAn,in.VAT,in.VAm,in.srat);

% Error to return
if isempty(errorratio)
	error=sqrt(Vz(2,2));    % use alpha
else
	% Now change coordinates to get variance of ratio we're interested in
	newVAN=changedenomcov(in.AN,VAN,in.di,errorratio(2));
	isonums=1:ISODATA.(element).nisos;
	newAni=isonums(isonums~=errorratio(2));  % numerators
	erat=find(errorratio(1)==newAni);        % the ratio that we're interested in
	error=sqrt(newVAN(erat,erat));
end

if (isempty(errorratio))
	ppmperamu=(1e6*error)./mean(ISODATA.(element).mass);
else
	stdratio=ISODATA.(element).standard(errorratio(1))/ISODATA.(element).standard(errorratio(2));
	massdiff=abs(ISODATA.(element).mass(errorratio(1))-ISODATA.(element).mass(errorratio(2)));
	ppmperamu=(1e6*error)./(stdratio*massdiff);
end
