function cosmo(filename,elements)
% COSMO   Generates a version of the cocktail list (Rudge et al., 2009)
% that takes into account the splitting of the sample between spiked and
% unspiked measurements 
%
% Input parameters are similar to the original cocktail list, save for the
% lack of spike 'type' since this only works with 'real' spikes (i.e., all
% isotopes of the element are present in the tracer). 
%
%  COSMO(filename,elements)
%             filename -- file to store output (a CSV file).  Default is 'cosmococktail.csv'.
%             elements -- which elements to include in the cookbook. Specify as a cell array
%                e.g. {'Ca','Fe'}. Default is all possible elements.
%
% This generates an exhaustive list of all possible double spikes for specified elements
% sorted in order of error. To improve efficiency of the optimization, the
% code only tests out (1) four-isotope inversion combinations that contain the
% spiked isotopes, and (2) normalization isotopes whose denominator is part
% of the four inversion isotopes. 
%
% Note that a number of parameters are specified in the global variable ISODATA.
%
% Example
%   cosmo('cosmoBa.csv',{'Ba'})
%
% See also dsstartup
global ISODATA

% default argument
if isempty(ISODATA)
	dsstartup;
end

if (nargin<1) || isempty(filename)
	filename=['cosmo.csv'];
end
if (nargin<2)
	elements=fieldnames(ISODATA);
end

% check that shake has been ran on the elements of interest
errel = {};
for i=1:length(elements)
	if ~isfield(ISODATA.(elements{i}).errormodel,'V100')
        errel{end+1} = elements{i};
    end
end

if ~isempty(errel)
    disp(['Run shake.m for the ff:' errel])
    return
end

disp(['Writing to ' filename]);

title=['Cosmo cocktail list:'];

fwritecell(filename,'%s','w',{title});
fwritecell(filename,'%s','a',{''});

for i=1:length(elements)
	element=elements{i};
	disp(element);
	in=ISODATA.(element);

    [optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu,optsplit,optnorm]=optimalspike(element);
    optisoinv=in.isoindex(optisoinv);
    optisonams=[{in.isoname{optisoinv(:,1)}}' {in.isoname{optisoinv(:,2)}}' {in.isoname{optisoinv(:,3)}}' {in.isoname{optisoinv(:,4)}}'];
    optnormnams=[{in.isoname{optnorm(:,1)}}' {in.isoname{optnorm(:,2)}}'];

    % write output to file
    isohead=strcat(repmat({'iso'},4,1),num2str((1:4)'))';

    spikehead=strcat(repmat({'spike'},in.nspikes,1),num2str((1:in.nspikes)'))';
    
    output=[optisonams num2cell([optspike optspikeprop optprop (1-optprop) optsplit opterr optppmperamu]) optnormnams];
    header=[isohead in.isoname spikehead {'spike'} {'sample'} {'split'} {'error'} {'ppmperamu'} {'normnum'} {'normden'}];
    fwritecell(filename,[repmat('%s,',1,4) repmat('%s,',1,in.nisos) repmat('%s,',1,in.nspikes) '%s,%s,%s,%s,%s,%s,%s'],'a',header);
    fwritecell(filename,[repmat('%s,',1,4) repmat('%f,',1,in.nisos) repmat('%f,',1,in.nspikes) '%f,%f,%f,%f,%f,%s,%s'],'a',output);
    
    fwritecell(filename,'%s','a',{''});
		
end
end

function [optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu,optsplit,optnorm]=optimalspike(element,beta,alpha,errorratio,isospike,isoinv,INisos)
%OPTIMALSPIKE    Finds the best real spike - based on the
%optimalrealspike.m subrouting from Rudge et al. (2009)
%    OPTIMALSPIKE(rawdata,beta,alpha,errorratio,isospike,isoinv)
%             rawdata -- data about a particular element
%             beta -- instrumental fractionation
%             alpha -- natural fractionation
%             errorratio -- the ratio whose error we are targeting
%             isospike -- the isotopes to spike
%             isoinv -- the isotopes used in the inversion

global ISODATA
rawdata=ISODATA.(element);

% Have some default arguments
if (nargin<7) || isempty(INisos)
	INisos=[];
end
if (nargin<6) || isempty(isoinv)
	isoinv=[];
end
if (nargin<5) || isempty(isospike)
	isospike=[];
end
if (nargin<4) || isempty(errorratio)
	errorratio=[];
end
if (nargin<3) || isempty(alpha)
	alpha=0;
end
if (nargin<2) || isempty(beta)
	beta=0;
end

% Convert isotope mass numbers to index numbers
errorratio=rawdata.isoindex(errorratio);
isospike=rawdata.isoindex(isospike);
isoinv=rawdata.isoindex(isoinv);

if (isempty(isoinv))
	isoinv=combnk(1:rawdata.nisos,4);
end

if isempty(isospike)
	isospikev=combnk(1:rawdata.nspikes,2);
else
	isospikev=isospike;
end

if isempty(INisos)
	INisosv=combnk(1:rawdata.nspikes,2);
    INisosv=sortrows([INisosv; INisosv(:,2) INisosv(:,1)],2);
else
	INisosv=INisos;
end

isoinvvals=[];            % inverse values to check
isospikevals=[];          % spike values to check
isonormvals=[];           % internal normalization values to check

for (i=1:size(isoinv,1))
    isospikev = combnk(isoinv(i,:),2); % spiked isotopes are always part of the inversion
    for (j=1:size(isospikev,1))
        INv = INisosv(find(ismember(INisosv(:,2),isoinv(i,:))),:); % internal norm denominator always part of inversion
        isonormvals=[isonormvals; INv];
        isospikevals=[isospikevals; repmat(isospikev(j,:),size(INv,1),1)];
        isoinvvals=[isoinvvals; repmat(isoinv(i,:),size(INv,1),1)];
    end
end


for i=1:size(isoinvvals,1)
	[optspike(i,:),optprop(i,:),opterr(i,:),optspikeprop(i,:),optppmperamu(i,:),optsplit(i,:),optnorm(i,:)]=singleoptimalspike(element,beta,alpha,errorratio,isospikevals(i,:),isoinvvals(i,:),isonormvals(i,:));
end
optisoinv=isoinvvals;

% Sort in ascending order of error
[opterr,ix]=sort(opterr);
optppmperamu=optppmperamu(ix,:);
optspike=optspike(ix,:);
optprop=optprop(ix,:);
optisoinv=optisoinv(ix,:);
optisoinv=rawdata.isonum(optisoinv);
optspikeprop=optspikeprop(ix,:);
optsplit=optsplit(ix,:); % optimal splitting of sample
optnorm=optnorm(ix,:); % optimal normalization isotopes

end

function [optspike,optprop,opterr,optspikeprop,optppmperamu,optsplit,optnorm]=singleoptimalspike(element,beta,alpha,errorratio,isospike,isoinv,INisos)
% Calculate the composition of the optimal double spike given the isotopes used in the inversion
% and of those the isotopes we are spiking

global ISODATA
rawdata=ISODATA.(element);

spikevector1=rawdata.rawspike(isospike(1),:);
spikevector2=rawdata.rawspike(isospike(2),:);

if (verLessThan('optim','4.0'))
	options=optimset('Display','notify','TolX',1e-8,'TolFun',1e-10,'TolCon',1e-6,'LargeScale','off','MaxFunEvals',10000);
else
	options=optimset('Display','notify','TolX',1e-8,'TolFun',1e-10,'TolCon',1e-6,'Algorithm','active-set','MaxFunEvals',10000);
end

tol=1e-5;
lb=[tol;
    tol;
    tol];

ub=[1-tol;
    1-tol;
    1-tol];

y0=[0.5 0.5 0.5]';

% Helpful to rescale the error, to make everything roughly order 1 for the optimiser
initialerror=errorwsplit(element,y0(1),y0(2),y0(3).*spikevector1+(1-y0(3)).*spikevector2,isoinv,INisos,errorratio,beta,alpha);

[y,opterr,exitflag] = fmincon(@(y) errorwsplit(element,y(1),y(2),y(3).*spikevector1 +(1-y(3)).*spikevector2,isoinv,INisos,errorratio,beta,alpha)./initialerror,y0,[],[],[],[],lb,ub,[],options);
[opterr,optppmperamu]=errorwsplit(element,y(1),y(2),y(3).*spikevector1 +(1-y(3)).*spikevector2,isoinv,INisos,errorratio,beta,alpha); 

optsplit=y(1);
optprop=y(2);
optspike=y(3)*spikevector1 +(1-y(3))*spikevector2;
optnorm=INisos;

optspikeprop=zeros(rawdata.nspikes,1);
optspikeprop(isospike(1))=y(3);
optspikeprop(isospike(2))=1-y(3);
end    
