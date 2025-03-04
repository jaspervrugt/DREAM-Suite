% Verifies whether there are any correlations between LM estimates and
% other measured data of samples

load SWIG_1D
% Rename content of file
data_SWIG = SWIG_1D;
% Unpack the data
n_soil = size(data_SWIG,1);

% Load the SWIG database
[NUM,TXT,RAW] = xlsread('SWIG database.xlsx','Metadata','A3:AG5025');
% Load results of LM fitting
load SWIG_646_LM.mat
% Load data of
for soil_type = 1:n_soil
    % Extract sample
    str_s = data_SWIG{soil_type,2};
    % Now extract the right columns
    Dat(soil_type,:) = NUM(str_s,:);
end

for par = 1:3
    for i = 2:size(Dat,2)
        % Determine correlation coefficient
        ii = ~isnan(Dat(:,i));
        R2 = corrcoef(opt_LM_SWIG(ii,par,1),Dat(ii,i));
        switch par
            case 1
                out1(i,1:3) = [ i R2(1,2) sum(ii) ];
            case 2
                out2(i,1:3) = [ i R2(1,2) sum(ii) ];
            case 3
                out3(i,1:3) = [ i R2(1,2) sum(ii) ];
        end
        % %     plot(opt_LM_SWIG(:,3,1),Dat(:,i),'r.');
        % %     i
        % %     pause;
    end
    % Build linear regression model
    switch par
            case 1
                A = out1;
            case 2
                A = out2;
            case 3
                A = out3;
    end
    A(:,2) = abs(A(:,2)); 
    As = -sortrows(-A,2);
    T1 = 6; T2 = 7;%T2 = 5;
    % Pick M predictors with highest R2
    idx_p = As(T1:T2,1); M = numel(idx_p);
    % Now check data they have in common
    idx_data = find(sum(~isnan(Dat(:,idx_p)),2) == M); numel(idx_data)
    % Regression
    Y = opt_LM_SWIG(idx_data,par,2);
    X = [ Dat(idx_data,idx_p) ones(numel(idx_data),1) ]; 
    [B,BINT,R,RINT,STATS] = regress(Y,X);
    STATS
    % the R-square statistic, the F statistic and p value
    % for the full model, and an estimate of the error variance.
end