function errorsurface(element,spike,isoinv,INisos,errorratio,alpha,beta,resolution,threshold,ncontour,plottype,varargin)
%ERRORSURFACE    A 2D contour plot of error as a function of double spike composition and double spike-sample proportions
%  ERRORSURFACE(element,spike,isoinv,INisos,alpha,beta,resolution,threshold,ncontour,plottype,...)
%             element -- element used in double spike, e.g. 'Fe'
%             spike -- the composition of the double spike as a composition vector e.g. [0.1 0.1 0.4 0.4]
%                represents a ~50-50 mixture of the third and fourth isotopes (57-58 for Fe) w some impurities.
%             isoinv -- the isotopes used in the inversion, e.g. [54 56 57 58].
%                By default the first four isotopes are chosen.
%             INisos -- the isotopes used for the internal normalization,
%                assumes [n d]. 
%             errorratio -- by default, the error on the natural fractionation
%                factor (known as alpha) is given. Instead, the error on a
%                 particular ratio can be given by setting errorratio. e.g.
%                setting errorratio=[58 56] will give the error on 58Fe/56Fe.
%             alpha, beta -- there is a small dependance of the error on the fractionation
%                factors (natural and instrumental). Values of alpha and
%                beta can be set here if desired, although the effect on the optimal spikes
%                is slight unless the fractionations are very large. Default is zero.
%             resolution -- number of grid points in x and y. Default is 100.
%             threshold -- maximum contour to plot, relative to the minimum error.
%                Default is 0.25 i.e. 25% in excess of the minimum.
%             ncontour -- number of countours. Default is 25.
%             plottype -- by default, the error is plotted. By setting this to 'ppmperamu'
%                an estimate of the ppm per amu is plotted instead.
%             ... -- additional arguments are passed to contour command.
%
% Note that a number of parameters are specified in the global variable ISODATA.
%
% Example
%    errorsurface('Ba',[0.1 0.1 0.4 0.4],[54 56 57 58],[54 56])
%
% See also errorwsplit, contour
global ISODATA

% Set some default values
if isempty(ISODATA)
	dsstartup;
end
if (nargin<11) || isempty(plottype)
	plottype='default';
end
if (nargin<10) || isempty(ncontour)
	ncontour=25;
end
if (nargin<9) || isempty(threshold)
	threshold=0.25;
end
if (nargin<8) || isempty(resolution)
	resolution=100;
end
if (nargin<7) || isempty(beta)
	beta=0;
end
if (nargin<6) || isempty(alpha)
	alpha=0;
end
if (nargin<5) || isempty(errorratio)
	errorratio=[];
end
if (nargin<4) || isempty(INisos)
	INisos=[];
end
if (nargin<3) || isempty(isoinv)
	isoinv=[1 2 3 4];
end

rawdata=ISODATA.(element);

% Convert isotope mass numbers to index numbers
errorratio=rawdata.isoindex(errorratio);
isoinv=rawdata.isoindex(isoinv);
INisos=rawdata.isoindex(INisos);

pvals=linspace(0.001,0.999,resolution);
svals=linspace(0.001,0.999,resolution);

[iv jv]=meshgrid(1:length(pvals),1:length(svals));
[errvals ppmperamuvals]=arrayfun(@(i,j) errorwsplit(element,svals(j),pvals(i),spike,isoinv,INisos,errorratio,alpha,beta),iv,jv);

% Getting minimum uncertainties in alpha and ppmperamu
[opterr,opterri]=min(errvals,[],'all','linear');
optppmperamu=min(ppmperamuvals,[],'all');

if strcmp(plottype,'ppmperamu')
	C=contour(pvals,svals,ppmperamuvals,linspace(optppmperamu,(1+threshold)*optppmperamu,ncontour+1),varargin{:});
else
	C=contour(pvals,svals,errvals,linspace(opterr,(1+threshold)*opterr,ncontour+1),varargin{:});
end

xlim([0 1]);
ylim([0 1]);
xlabel('proportion of double spike in double spike-sample mix');
ylabel('proportion of sample going to double spike measurement');

invisostring=[rawdata.isolabel{isoinv(1)} ', ' rawdata.isolabel{isoinv(2)} ', ' rawdata.isolabel{isoinv(3)} ', ' rawdata.isolabel{isoinv(4)} ' inversion'];
if isempty(errorratio)
	title(['Error in \alpha (1SD) with ' invisostring]);
else
	title(['Error in ' rawdata.isolabel{errorratio(1)} '/' rawdata.isolabel{errorratio(2)} ' (1SD) with ' invisostring]);
end

hold on;
plot(pvals(iv(opterri)),svals(jv(opterri)),'rx');
hold off;

