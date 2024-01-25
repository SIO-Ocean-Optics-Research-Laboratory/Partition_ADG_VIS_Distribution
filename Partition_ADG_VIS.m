function [lambda_out,adopt,agopt] = Partition_ADG_VIS(lambda_adg,adg,PT)
%Implements the absorption partitioning model to partition the
%non-phytoplankton absorption coefficient, adg, into its two constituents:
%absorption due to depigmented (non-algal) particles, ad, and absorption
%due to chromophoric dissolved organic matter, ag. The model operates in
%the visible (VIS) spectral region covering the light wavelength range from
%400 to 700 nm.
%
%Reference:
%
%Kehrli M. D., Stramski D., Reynolds R. A., Joshi I. D., A model for
%partitioning the non-algal absorption coefficient of seawater in the
%ultraviolet and visible spectral range into the contribution of non-algal
%particulate and dissolved matter. Submitted to Applied Optics January 24, 
%2024 (Referenced in code documentation as KSRJ).
%
%
%Required function inputs: lambda_adg, adg, PT
%   lambda_adg [m-by-1 numeric]: Input values of light wavelength [nm]
%   corresponding to the input values of spectral absorption coefficient
%   adg.
%
%   adg [m-by-1 numeric]: Input values of spectral absorption coefficient
%   adg [m^-1].
%
%   PT [1-by-n numeric]: Desired percentile values of feasible solution
%   pool to return (Optional - if not specified the 10th and 90th
%   percentile values of feasible solutions are returned in addition to the
%   optimal solutions).
%
%Required ancillary data file: VIS_lib.mat
%   VIS_lib.mat: Visible (VIS) spectral shape function library in the
%   spectral range of 400-700 nm required to be in MATLAB path for model to
%   operate. Derived from joint probability distributions of spectral
%   steepness parameters from the development dataset described in the
%   Reference Manuscript.
%
%Outputs: lambda_out, adopt, agopt
%   lambda_out [301-by-1 double]: Light wavelengths [nm] corresponding to
%   output values of spectral absorption coefficients ad and ag from the
%   partitioning model.
%
%   adopt [301-by-(1+length(PT)) double]: Matrix containing the optimal
%   solutions for spectral values of ad [m-1] determined as the median
%   value of all feasible solutions at each wavelength (1st column) and
%   additional desired percentile solutions from the pool of feasible
%   solutions (remaining output column(s)).
%
%   agopt [301-by-(1+length(PT)) double]: Matrix containing the optimal
%   solutions for spectral values of ag [m^-1] determined as the median
%   value of all feasible solutions at each wavelength (1st column) and
%   additional desired percentile solutions from the pool of feasible
%   solutions (remaining output column(s)).
%
%Version 1.0 (v1.0)
%
%Version history: 
%2023-11-09: Revised ADG partitioning model and Matlab version, M. D.
%Kehrli, D. Stramski, R. A. Reynolds, I. D. Joshi.
%2024-01-25: Final revised MATLAB version (v1.0), M. Kehrli, D. Stramski,
%R. A. Reynolds, and I. D. Joshi
%
%Adapted from: GSCM_insitu_adg_final.m (2019-01-28), original ADG
%partitioning model and Matlab version, L. Li, R. A. Reynolds, D. Stramski,
%described in Stramski, Li, Reynolds, 2019, Applied Optics, 58, 3790-3806.
%[doi: 10.1364/AO.58.003790].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input arguments and set defaults for optional input parameters
arguments
    lambda_adg (:,1) {mustBeNumeric}
    adg (:,1) {mustBeNumeric}
    %default to 10th and 90th percentile values if not specified by user
    PT (1,:) {mustBeNumeric} = [10 90]
end

%number of percentile values to include in output
nPT = size(PT,2);

%load library of spectral shape functions, adhat and aghat, which
%characterize the variation in the spectral shapes of ad and ag
%coefficients, respectively, in the VIS spectral range. Note that .mat file
%must be in function path to load to workspace.
shapelib = load('VIS_lib.mat');
%get the row indices where input adg has valid values between 400 and 700
%nm
idx_dg = find(~isnan(adg) & 400 <= lambda_adg & lambda_adg <= 700);

%get the corresponding spectral row indices for the ad and ag spectral
%shapes where adg has valid values
[~, idx_shape] = ismember(lambda_adg(idx_dg), shapelib.lambda);

%calculate the spectral shape function adghat of the input adg spectrum
%(KSRJ Eq. 1)
adghat = numel(idx_dg)*adg/sum(adg(idx_dg),'omitnan'); 

%initialize a counter for looping through ad and ag shape functions
count = 1;
%initialize arrays for holding all speculative solutions, adspec and
%agspec, for ad and ag coefficients
adspec_all = [];
agspec_all = [];

%% Calculation of all speculative solutions
%for each available ad shape function in the library, loop through all ag
%shape functions and build matrix of speculative solutions

%turn off warning regarding rank deficient matrix
warning('off','MATLAB:rankDeficientMatrix'); 

%start of loop for all ad shape functions
for idx_adhat=1:size(shapelib.ad,2)
    %get a single ad shape function and all ag shape functions from the
    %library
    adhat = shapelib.ad(idx_shape,idx_adhat);
    aghat = shapelib.ag(idx_shape, :);
    
    %preallocate array for ag weighting factors
    ag_weight = nan(1,size(aghat,2)); 
    
    %calculate the ag weighting factor for each ag shape function
    for i=1:size(aghat,2)
        %see Eq. 7 in KSRJ
        ag_weight(i) = (aghat(:,i)-adhat(:)) \ (adghat(idx_dg)-adhat(:));
    end
    
    %find and eliminate the speculative solutions for which ag weight is
    %not between 0 and 1
    idx_valid = 0 < ag_weight & ag_weight < 1;
    ag_weight = ag_weight(idx_valid);
    aghat = aghat(:,idx_valid);
    %number of speculative solutions for a given ad shape function and all
    %valid ag shape functions
    numspec = numel(ag_weight); 
    
    %if there are no speculative solutions found, skip to next ad shape
    %function
    if numspec == 0
        continue
    end
    
    %calculate the speculative solutions, adspec, for ad
    ad_weight = 1 - ag_weight;
    %see Eq. 8 in KSRJ
    adspec = (1/numel(idx_dg)) * repmat(ad_weight,numel(idx_shape),1) .* repmat(adhat,1,numspec) * sum(adg(idx_dg),'omitnan');
    
    %calculate the speculative solutions, agspec, for ag
    %see Eq. 9 in KSRJ
    agspec = (1/numel(idx_dg)) * repmat(ag_weight,numel(idx_shape),1) .* aghat * sum(adg(idx_dg),'omitnan');
    
    %append all valid speculative solutions for this specific ad shape
    %function to main output array for all speculative solutions
    adspec_all(:,count:count+numspec-1) = adspec;
    agspec_all(:,count:count+numspec-1) = agspec;
   
    %increment counter and repeat loop for next ad shape function
    count = count + numspec;
end
%reset the warning back to its default
warning('on','MATLAB:rankDeficientMatrix');

%if there are no speculative solutions found for ad and ag, return NaN for
%the optimal and range of feasible solutions of partitioning model
if isempty(adspec_all)
    adopt = nan(numel(lambda_adg),nPT+1);
    agopt = nan(numel(lambda_adg),nPT+1);
    return
end

%% Determination of feasible and optimal solutions for ad and ag coefficients

%calculate residual statistics between measured adg and speculative values
%of adg obtained as a sum of all pairs of ordinary residuals for the
%speculative solutions of ad and ag
Resids = (adspec_all + agspec_all) - repmat(adg(idx_dg), 1, size(adspec_all,2));
%residuals normalized to measured adg
NResids = Resids./repmat(adg(idx_dg),1,size(adspec_all,2));
%sum of normalized squared residuals for wavelengths <=600 nm - see Eq. 10
%in KSRJ
NSSR = sum(NResids(lambda_adg(idx_dg)<=600,:).^2,1,'omitnan');

%indices of feasible solutions 
idx_feas = find(NSSR <= prctile(NSSR,10));

%feasible solutions, adfeas_allf and agfeas_all, for ad and ag coefficients
%selected from the pool of valid speculative solutions
adfeas_allf = adspec_all(:,idx_feas); 
agfeas_allf = agspec_all(:,idx_feas);

%calculation of optimal solutions, adopt and agopt, for ad and ag
%coefficients (the optimal solution is the median value of all feasible
%solutions at each output wavelength) and the percentile values
%characterizing the range of feasible solutions
adopt(:,1) = median(adfeas_allf,2,'omitnan');
agopt(:,1) = median(agfeas_allf,2,'omitnan');
adopt(:,2:nPT+1) = prctile(adfeas_allf, PT, 2);
agopt(:,2:nPT+1) = prctile(agfeas_allf, PT, 2);
lambda_out = shapelib.lambda;
end
