function [f,str] = constructpolynomialmatrix3d_1(matrixsize,locs,degree,weights,str)

% function [f,str] = constructpolynomialmatrix3d(matrixsize,locs,degree,weights)
%
% <matrixsize> is a 3D matrix size like [100 50 100]
% <locs> is a row or column vector of indices into that matrix size
% <degree> is the maximum polynomial degree desired
% <weights> (optional) is a 1 x N vector of values.
%   if supplied, we automatically weight and sum the basis functions.
%   the point of this input is to avoid having to explicitly create
%   all the basis functions (and therefore we save on memory requirements).
%   default is [] which means do nothing special.
%
% if <weights> is not supplied, then:
%   return <f>, a matrix of dimensions length(<locs>) x N
%   with polynomial basis functions evaluated at <locs> in
%   the columns.  the polynomial basis functions are evaluated
%   over the range [-1,1] which is presumed to correspond to
%   the beginning and ending element along each of the three dimensions.
%   (if a dimension has only one element, the values are all set to 1.)
%   also, return <str>, the algebraic expression that corresponds to
%   the columns of <f>.  'x' refers to the first matrix dimension; 'y'
%   refers to the second matrix dimension; 'z' refers to the third
%   matrix dimension.
%
% if <weights> is supplied, then:
%   return <f>, a vector of dimensions length(<locs>) x 1 with the weighted
%   sum of the polynomial basis functions.  also, return <str>, a cell
%   vector of algebraic expressions describing the various basis functions.
%
% if <str> is supplied, then:
%   return <f>, a vector of dimensions length(<locs>) x 1. will be created
%   acourding to the given str is not the str we be created by the function.
%   note that differnt computer we make a differnt str (in respect  order of parametersather).
%   The input str will control this differences
%
% note that there may be various gain factors on the basis functions
% (e.g. don't assume that they are all unit-length).
%
% also, beware of numerical precision issues for high degrees...
%
% see also constructpolynomialmatrix2d.m.
%
% example:
% [f,str] = constructpolynomialmatrix3d([30 30 30],find(ones(30,30,30)),2);
% str
% f = reshape(f,30,30,30,[]);
% for p=1:size(f,4)
%   drawnow; figure; imagesc(makeimagestack(f(:,:,:,p)));
% end

% input
if ~exist('weights','var') || isempty(weights)
    weights = [];
end

if exist('str','var') && ~isempty(str)
    
    % prep the linear coordinates
    [x,y,z] = ind2sub(matrixsize,locs(:));
    if matrixsize(1)~=1
        x = normalizerange(x,-1,1,1,matrixsize(1));
    end
    if matrixsize(2)~=1
        y = normalizerange(y,-1,1,1,matrixsize(2));
    end
    if matrixsize(3)~=1
        z = normalizerange(z,-1,1,1,matrixsize(3));
    end
    
    % handle regular case
    if ~isempty(weights)
        error ('not support weights polynomyals and input str')
    else
        % do it
        f = eval(str);
        
    end
else
    % prep
    x = sym('x');
    y = sym('y');
    z = sym('z');
    
    % do the algebra
    str = char(expand((x+y+z+1)^degree));
    
    % sort the stuff in between + signs to try to ensure consistent ordering!!!
    str = sort(strsplit(str,'+'));
    str = cat(1,str,repmat({'+'},[1 length(str)]));
    str = cat(2,str{:});
    str = str(1:end-1);
    
    % add a little padding so the 1 step below will work for degree 0
    str = [' ' str ' '];
    
    % change * to .*
    old = pwd;  % THIS IS A TOTAL HACK TO AVOID FUNCTION NAME CONFLICTS!
    cd(fullfile(matlabroot,'toolbox','matlab','funfun'));
    str = vectorize(str);
    cd(old);
    
    % remove +
    str(str=='+') = ' ';
    
    % change 1 to ones(size(x),1)
    str0 = strrep(str,' 1 ',' ones(size(x)) '); assert(length(str0) ~= length(str));
    str = str0;
    
    % prep the linear coordinates
    [x,y,z] = ind2sub(matrixsize,locs(:));
    if matrixsize(1)~=1
        x = normalizerange(x,-1,1,1,matrixsize(1));
    end
    if matrixsize(2)~=1
        y = normalizerange(y,-1,1,1,matrixsize(2));
    end
    if matrixsize(3)~=1
        z = normalizerange(z,-1,1,1,matrixsize(3));
    end
    
    % handle regular case
    if isempty(weights)
        
        % add brackets
        str = [ '[' str ']' ];
        
        % do it
        f = eval(str);
        
        % handle special case
    else
        
        % divide them up
        str = strsplit(str,' ');
        
        % throw away empty elements
        str = str(cellfun(@(x) ~isempty(x),str));
        
        % weight and sum
        f = 0;
        for p=1:length(str)
            f = f + eval(str{p})*weights(p);
        end
        
    end
    
end







%   switch degrees(p)
%   case 0
%     f = [f ones(size(x))];
%   case 1
%     f = [f x y z];
%   case 2
%     f = [f x.^2 y.^2 z.^2 x.*y x.*z y.*z];
%   case 3
%     f = [f x.^3 y.^3 z.^3 x.^2.*y x.*y.^2 x.^2.*z x.*z.^2 y.^2.*z y.*z.^2 x.*y.*z];
%   case 4
%     f = [f x.^4 y.^4 z.^4 x.^3.*y x.*y.^3 x.^3.*z x.*z.^3 y.^3.*z y.*z.^3 x.^2.*y.^2 x.^2.*z.^2 y.^2.*z.^2 x.^2.*y.*z x.*y.^2.*z x.*y.*z.^2];
%   otherwise
%     die;
%   end

% FOURIER STUF
% % <maxcpfov> is the maximum allowable cycles per FOV (measured
% %   along the first dimension)
%
% % all basis functions
% % have unit length and have mean 0 (except for the DC basis
% % function).
%
%         % not handling nyquist
%         % get rid of bogus basis functions/ DC special
%         % normalize basis functions
%
% %
% % % do it
% % f = [];
% % for a=0:maxcpfov
% %   for b=0:maxcpfov
% %     for c=0:maxcpfov
% %       cpfov = sqrt(a.^2 + b.^2 + c.^2);
% %       if cpfov <= maxcpfov
% %
% %         tt = a*xx + b*yy + c*zz;
% %         if a==0 & b==0 & c==0
% %           f(:,end+1) = cos(-2*pi*tt) / sqrt(2);
% %         else
% %           f(:,end+1) = cos(-2*pi*tt);
% %           f(:,end+1) = sin(-2*pi*tt);
% %         end
% %
% %       end
% %     end
% %   end
% % end
