function V=calcratiocovIN(element,composition,errormodel,INisos)
% Calculate the covariance matrix of the ratios with internal normalization
% errormodel is the sub-structure for the element of interest that has 
% INisos are the isotopes used for the internal normalisation [n d]

global ISODATA


% Convert input of normalization isotopes to indeces
INisos=ISODATA.(element).isoindex(INisos); 

% Get rho values for the element of interest
m_i=ISODATA.(element).mass;
m_j=m_i(INisos(1));
m_k=m_i(INisos(2));

rho_i=log(m_i./m_k)./log(m_j/m_k);

% first normalise composition so it is really a composition (unit sum)
composition=composition./sum(composition);

meanbeams=composition.*errormodel.intensity;

covbeams=calcbeamcov(meanbeams,errormodel);
V=covbeamtoratioIN(meanbeams,covbeams,INisos,rho_i);

function beamcov=calcbeamcov(meanbeams,errormodel)
% the beam covariance matrix
beamvar=errormodel.a + meanbeams.*errormodel.b + (meanbeams.^2).*errormodel.c;
beamcov=diag(beamvar);

function V=covbeamtoratioIN(meanbeams,covbeams,INisos,rho)
% converts a covariance matrix for beams to one for ratios
% INisos
% assumes last row and column of M correspond to denominator
di=INisos(2);
isonums=1:length(meanbeams);
ni=isonums(isonums~=di);
[~,ninorm]=ismember(INisos(1),ni); % index of the numerator for the normalization amongst the non-denominator isotopes (ni)

n=meanbeams(ni);
d=meanbeams(di);
r=rho(ni);

M=covbeams([ni di],[ni di]);  % move denominator to end
%A=[diag(repmat(1/d,1,length(n))) -n'./(d^2)];
%A=[(1/d).*eye(length(n)) -n'./(d^2)];
A=[diag((1/d).*ones(1,length(n))) n'./(d^2).*(r-1)'];
A(:,ninorm)=-r.*n./(n(ninorm)*d);
A(ninorm,:)=0;
V=(A*M)*(A');
