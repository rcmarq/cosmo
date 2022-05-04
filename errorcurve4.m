function errorcurve4(element,split,spike,isoinv,INisos,epsisos,errorratio,alpha,beta,plottype,varargin)
%ERRORCURVE4    A plot of error as a function of splitting between the spiked (ID) and
% unspiked (IC) measurement for a given DS scheme. 
%  ERRORCURVE3(element,type,prop,spike,isoinv,errorratio,alpha,beta,...)
%             element -- element used in double spike, e.g. 'Fe'
%             split -- the splitting of the sample, i.e., the fraction of
%                total sample going to the double spike measurement.
%             spike -- the composition of the double spike as a composition vector e.g. [0 0 0.5 0.5]
%                represents a 50-50 mixture of the third and fourth isotopes (57-58 for Fe).
%             isoinv -- the isotopes used in the inversion, e.g. [54 56 57 58].
%                By default the first four isotopes are chosen.
%             INisos -- the isotopes used for the internal normalization,
%                assumes [n d]. 
%             epsios -- isotopes whose uncertainties in the unspiked
%                measurement is also shown. In epsilon (part per  10^4) units.
%             errorratio -- by default, the error on the natural fractionation
%                factor (known as alpha or alpha) is given. Instead, the
%                error on a particular ratio can be given by setting errorratio. e.g.
%                setting errorratio=[58 56] will give the error on 58Fe/56Fe.
%             alpha, beta -- there is a small dependance of the error on the fractionation
%                factors (natural and instrumental). Values of beta and
%                alpha can be set here if desired, although the effect on the optimal spikes
%                is slight unless the fractionations are very large. Default is zero.
%             plottype -- by default, the error is plotted. By setting this to 'ppmperamu'
%                an estimate of the ppm per amu is plotted instead.
%             ... -- additional arguments are passed to the plot command.
%
% Note that a number of parameters are specified in the global variable ISODATA.
%
% Example
%    errorcurve4('Fe',0.5,[0.5 0.5 0.5 0.5],[54 56 57 58],[54 56])
%
% See also errorwsplit

global ISODATA

if isempty(ISODATA)
	dsstartup;
end
if (nargin<10) || isempty(plottype)
	plottype='default';
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
if (nargin<6) || isempty(epsisos)
	epsisos=[];
end
if (nargin<5) || isempty(INisos)
	INisos=[];
end
if (nargin<4) || isempty(isoinv)
	isoinv=[1 2 3 4];
end

rawdata=ISODATA.(element);
rawspike=rawdata.rawspike;

% check that shake has been ran on the elements of interest

if ~isfield(ISODATA.(element).errormodel,'V100')
    disp(['The shake.m function has not been ran for this element.'])
    return;
end

% Convert isotope mass numbers to index numbers
errorratio=rawdata.isoindex(errorratio);
isoinv=rawdata.isoindex(isoinv);
INisos=rawdata.isoindex(INisos);
epsisos=rawdata.isoindex(epsisos);

pvals=linspace(0.001,0.999,1000);
errvals=zeros(size(pvals));
ppmperamuvals=zeros(size(pvals));

for i=1:length(pvals)
	[errvals(i) ppmperamuvals(i)]=errorwsplit(element,split,pvals(i),spike,isoinv,INisos,errorratio,alpha,beta);
end

if isequal(plottype,'ppmperamu')
	plotvals=ppmperamuvals;
else
	plotvals=errvals;
end

plot(pvals,plotvals,varargin{:},'DisplayName','ID Meas.');
mine=min(plotvals);
xlim([0 1]);
ylim([0 5*mine]);

xlabel(['proportion of double spike in double spike-sample mix']);

if isempty(errorratio)
	ylabel('Error in \alpha (1SE)');
else
	ylabel(['Error in ' rawdata.isolabel{errorratio(1)} '/' rawdata.isolabel{errorratio(2)} ' (1SD)']);
end
title([rawdata.isolabel{isoinv(1)} ', ' rawdata.isolabel{isoinv(2)} ', ' rawdata.isolabel{isoinv(3)} ', ' rawdata.isolabel{isoinv(4)} ' inversion' ])


