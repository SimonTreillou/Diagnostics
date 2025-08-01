function [varargout] = merge2netcdf(out,varargin)
% All input files must have the same structures
% Concatenated along the time dimension
    if ~ischar(out)
       error(' out must be a string ')
    end
    varargout = cell(1, nargout);
    info = ncinfo(varargin{1});
    ncwriteschema(out,info);
    disp(class(varargin))
    nccat(out,varargin{:});
end