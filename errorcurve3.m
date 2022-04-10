function errorcurve3(element,prop,spike,isoinv,INisos,epsisos,errorratio,alpha,beta,plottype,varargin)
%ERRORCURVE3    A plot of error as a function of splitting between the spiked (ID) and
% unspiked (IC) measurement for a given DS scheme. 
%  ERRORCURVE3(element,type,prop,spike,isoinv,errorratio,alpha,beta,...)
%             element -- element used in double spike, e.g. 'Fe'
%             prop -- proportion of double spike in double spike-sample mixture e.g. 0.5.
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
%    errorcurve3('Fe',0.5,[0.5 0.5 0.5 0.5],[54 56 57 58],[54 56])
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
if (nargin<5) || isempty(isoinv)
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

svals=linspace(0.001,0.999,1000);
errvals=zeros(size(svals));
ppmperamuvals=zeros(size(svals));

for i=1:length(svals)
	[errvals(i) ppmperamuvals(i)]=errorwsplit(element,svals(i),prop,spike,isoinv,INisos,errorratio,alpha,beta);
end

if isequal(plottype,'ppmperamu')
	plotvals=ppmperamuvals;
else
	plotvals=errvals;
end

plot(svals,plotvals,varargin{:},'DisplayName','ID Meas.');
mine=min(plotvals);
xlim([0 1]);
ylim([0 5*mine]);

xlabel(['proportion of the sample going to the double spike measurement']);

if isempty(errorratio)
	ylabel('Error in \alpha (1SD)');
else
	ylabel(['Error in ' rawdata.isolabel{errorratio(1)} '/' rawdata.isolabel{errorratio(2)} ' (1SD)']);
end
title([rawdata.isolabel{isoinv(1)} ', ' rawdata.isolabel{isoinv(2)} ', ' rawdata.isolabel{isoinv(3)} ', ' rawdata.isolabel{isoinv(4)} ' inversion' ])

if ~isempty(epsisos)
    
    epserrvals=[];
    cyclesIC=rawdata.errormodel.standard.cycles; % Number of cycles
    stdR=rawdata.standard./rawdata.standard(INisos(2)); % Get standard ratios
    stdR=stdR(1:rawdata.nisos~=INisos(2)); % Exclude identity ratio of norm. isotope
    
    for j=1:length(svals)
        % Update the voltage based o the splitting value
        sampleIC=rawdata.errormodel.V100*(1-svals(j))*rawdata.errormodel.standard.eff/cyclesIC; % total voltage of sample per cycle
        rawdata.errormodel.standard.intensity=sampleIC;
        
        % Solve for the covariance matrix of the standard
        V=calcratiocovIN(element,rawdata.standard,rawdata.errormodel.standard,INisos);
        
        % Converting to uncertainty in epsilon (parts per 10^4) units
        epserrvals(j,:)=sqrt(diag(V))'./stdR.*10000;
    end
    
    WAcolors = [0.7412    0.1216    0.1373
                0.8314    0.5686    0.3490
                0.2667    0.2824    0.4118
                0.4667    0.6353    0.5843
                0.4078    0.6000    0.6471
                0.6549    0.4863    0.4627
                0.8353    0.5255    0.4627
                0.4549    0.0510    0.1020
                0.2784    0.5294    0.6627];
    
    yyaxis right;
    minvals=[]; % minimum uncertainty values
    ci=1;
    
    hold on
    for k=epsisos
        plot(svals,epserrvals(:,k),'--','Color',WAcolors(ci,:),'DisplayName',['\epsilon ' rawdata.isolabel{k}]);
        minvals(end+1)=min(epserrvals(:,k));
        ci=ci+1;
        legend;
    end
    
    set(gca,'YColor','k');
    ylim([0 5*max(minvals)]);
    ylabel('Error in \epsilon (1SD)','Color','k');
    
    hold off
    
end
