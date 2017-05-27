function cscheme(varargin)
if nargin>0
    scheme = varargin{1};
else
    scheme = jet;
end
set(groot,'DefaultFigureColormap',scheme);
set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75;...
    .75 0 .75; .75 .75 0; .25 .25 .25]);
close
