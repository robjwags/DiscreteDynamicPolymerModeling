function [fitresult, gof] = createFits_kexc_wrtc(x_temp, y, nk, mt)
%CREATEFITS(X_TEMP,Y)
%  Create fits.
%
%  Data for 'bs12' fit:
%      X Input : x_temp
%      Y Output: y
%  Data for 'ms12' fit:
%      X Input : x_temp
%      Y Output: y
%  Data for 'bs18' fit:
%      X Input : x_temp
%      Y Output: y
%  Data for 'ms18' fit:
%      X Input : x_temp
%      Y Output: y
%  Data for 'bs36' fit:
%      X Input : x_temp
%      Y Output: y
%  Data for 'ms36' fit:
%      X Input : x_temp
%      Y Output: y
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 27-Mar-2024 16:32:13

%% Initialization.

% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell( 1, 1 );
gof = struct( 'sse', cell( 1, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

if nk==12 && mt==0
    %% Fit: 'bs12'.
    [xData, yData] = prepareCurveData( x_temp, y );
    
    % Set up fittype and options.
    ft = fittype( 'k*(c*x^(1/3)-1)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.898570798588752 0.73172238565867];
    
    % Fit model to data.
    [fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );
    
    % % Plot fit with data.
    % figure( 'Name', 'bs12' );
    % h = plot( fitresult{1}, xData, yData );
    % legend( h, 'y vs. x_temp', 'bs12', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % % Label axes
    % xlabel( 'x_temp', 'Interpreter', 'none' );
    % ylabel( 'y', 'Interpreter', 'none' );
    % grid on
    
elseif nk==12 && mt==1
    %% Fit: 'ms12'.
    [xData, yData] = prepareCurveData( x_temp, y );
    
    % Set up fittype and options.
    ft = fittype( 'k*(c*x^(1/3)-1)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.451582284453393 0.649115474956452];
    
    % Fit model to data.
    [fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );
    
    % % Plot fit with data.
    % figure( 'Name', 'ms12' );
    % h = plot( fitresult{1}, xData, yData );
    % legend( h, 'y vs. x_temp', 'ms12', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % % Label axes
    % xlabel( 'x_temp', 'Interpreter', 'none' );
    % ylabel( 'y', 'Interpreter', 'none' );
    % grid on
    
elseif nk==18 && mt==0
    %% Fit: 'bs18'.
    [xData, yData] = prepareCurveData( x_temp, y );
    
    % Set up fittype and options.
    ft = fittype( 'k*(c*x^(1/3)-1)', 'independent', 'x', 'dependent', 'y' );
    excludedPoints = excludedata( xData, yData, 'Indices', 2 );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.0845266111088704 0.168990029462704];
    opts.Exclude = excludedPoints;
    
    % Fit model to data.
    [fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );
    
    % % Plot fit with data.
    % figure( 'Name', 'bs18' );
    % h = plot( fitresult{1}, xData, yData, excludedPoints );
    % legend( h, 'y vs. x_temp', 'Excluded y vs. x_temp', 'bs18', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % % Label axes
    % xlabel( 'x_temp', 'Interpreter', 'none' );
    % ylabel( 'y', 'Interpreter', 'none' );
    % grid on
    
elseif nk==18 && mt==1
    %% Fit: 'ms18'.
    [xData, yData] = prepareCurveData( x_temp, y );
    
    % Set up fittype and options.
    ft = fittype( 'k*(c*x^(1/3)-1)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.769850513631213 0.0430238016578078];
    
    % Fit model to data.
    [fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );
    
    % % Plot fit with data.
    % figure( 'Name', 'ms18' );
    % h = plot( fitresult{1}, xData, yData );
    % legend( h, 'y vs. x_temp', 'ms18', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % % Label axes
    % xlabel( 'x_temp', 'Interpreter', 'none' );
    % ylabel( 'y', 'Interpreter', 'none' );
    % grid on
    
elseif nk==36 && mt==0
    %% Fit: 'bs36'.
    [xData, yData] = prepareCurveData( x_temp, y );
    
    % Set up fittype and options.
    ft = fittype( 'k*(c*x^(1/3)-1)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.0133663714160782 0.353158571222071];
    
    % Fit model to data.
    [fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );
    
    % % Plot fit with data.
    % figure( 'Name', 'bs36' );
    % h = plot( fitresult{1}, xData, yData );
    % legend( h, 'y vs. x_temp', 'bs36', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % % Label axes
    % xlabel( 'x_temp', 'Interpreter', 'none' );
    % ylabel( 'y', 'Interpreter', 'none' );
    % grid on
    
elseif nk==36 && mt==1
    %% Fit: 'ms36'.
    [xData, yData] = prepareCurveData( x_temp, y );
    
    % Set up fittype and options.
    ft = fittype( 'k*(c*x^(1/3)-1)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.291192867476437 0.234779913372406];
    
    % Fit model to data.
    [fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );
    
    % % Plot fit with data.
    % figure( 'Name', 'ms36' );
    % h = plot( fitresult{1}, xData, yData );
    % legend( h, 'y vs. x_temp', 'ms36', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % % Label axes
    % xlabel( 'x_temp', 'Interpreter', 'none' );
    % ylabel( 'y', 'Interpreter', 'none' );
    % grid on
end


