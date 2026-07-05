function result = case1_evaluate_scheme(dataset, scheme, options)
%CASE1_EVALUATE_SCHEME Evaluate one calibration scheme for Case1.
%   method = 'spline'  : MATLAB-style not-a-knot cubic spline.
%   method = 'natural' : accepted as an alias for the submit spline evaluator.

if nargin < 3
    options = struct();
end

if ~isfield(options, 'method') || isempty(options.method)
    options.method = 'spline';
end

switch lower(options.method)
    case 'spline'
        result = case1_evaluate_spline(dataset, scheme, options);
    case 'natural'
        result = case1_evaluate_spline(dataset, scheme, options);
    otherwise
        error('Unknown interpolation method: %s', options.method);
end
end

