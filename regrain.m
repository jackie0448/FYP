%
%  regraining post-process to match colour of IR and gradient of I0
%
%   IRR = degrain(I_original, I_graded, [smoothness]);
%
%   the smoothness (default=1, smoothness>=0) sets the fidelity of the 
%   original gradient field. e.g. smoothness = 0 implies IRR = I_graded.
%
%  (c) F. Pitie 2007
%
%  see reference:
%     Automated colour grading using colour distribution transfer. (2007)
%     Computer Vision and Image Understanding.
%
%  note: this implementation follows a simple top-down approach
%        with jacobi iterations.
%
%IR is the image after the n dims
function IRR = regrain(I0, IR, varargin)
%{
numvarargs = length(varargin);
if numvarargs > 1
    error('regrain:TooManyInputs', ...
        'requires at most 1 optional input');
end

optargs = {1};
optargs(1:numvarargs) = varargin;
[smoothness] = optargs{:};
%}
% IRR is a temp var to store the original image
IRR = I0;
[IRR] = regrain_rec(IRR, I0, IR, [4 16 32 64 64 64 ], 1, 0);

end


function [IRR] = regrain_rec(IRR, I0, IR, nbits, smoothness, level)

hres = size(I0,2); %give the number of columns/ aka in the horizontal direction
vres = size(I0,1); %give the number of rows
vres2 = ceil(vres/2);
hres2 = ceil(hres/2);

if (length(nbits) > 1 && hres2 > 30 && vres2 > 30)
    I02 = imresize(I0, [vres2 hres2], 'bilinear');
    IR2 = imresize(IR, [vres2 hres2], 'bilinear');
    IRR2 = imresize(IRR, [vres2 hres2], 'bilinear');
    [IRR2] = regrain_rec(IRR2, I02, IR2, nbits(2:end), smoothness, level+1);
    IRR = imresize(IRR2, [vres hres], 'bilinear');
end

IRR = solve(IRR, I0, IR, nbits(2), smoothness, level);

end

function IRR = solve(IRR, I0, IR, nbits, smoothness, level)

hres = size(I0,2); % gives number of cols
vres = size(I0,1); % gives the number of rows
K = size(I0,3);

%5 sets of xy value bc its the 4 neighbouring pixels oc the main one
y0 = 1:vres;%x,y
y1 = [1 1:vres-1]; %ok
y2 = [2:vres vres];
y3 = 1:vres; %ok
y4 = [1:vres-1 vres-1] ;

x0 = 1:hres; %x,y
x1 = [1: hres]; %ok
x2 = 1:hres;
x3 = [1 1:hres-1]; %ok
x4 = 1:hres;

G0 = I0;
%all the row, from 
%G0y = (G0(:,[2:end end], :) - G0(:,[1 1:end-1], :)); % just taking a variations of the column and subtracts
%G0x = (G0([2:end end], :, :) - G0([1 1:end-1], :, :));
[G0x, G0y]=  gradient(G0); %ok
dI = sqrt(sum(G0x.^2 + G0y.^2, 3)); %ok 

%h = 2^(-level);
psi = min(256*dI/5, 1); %based on formula in paper
phi = 30./(1 + 10*dI); %based on formula in paper, a weight field, so it has a range

a1 = (phi(y1,x1) + phi(y0, x0))/2;
a2 = (phi(y2,x2) + phi(y0, x0))/2;
a3 = (phi(y3,x3) + phi(y0, x0))/2;
a4 = (phi(y4,x4) + phi(y0, x0))/2;
%a5= 0.5.* (phi(y0,x0)+ phi(y1,x1) + phi(y2,x2)+ phi(y3,x3)+ phi(y4,x4))+ psi;

rho = 1/5;
%this is jacobi
while (nbits~= 0)
    den =  psi + a1 + a2 + a3 + a4;
    %repmat is to ensure the size of the matrix is the same in 3 dim, so 
    %repmat(a1, [1 1 k]) means you are not changing the size of the img,
    %but duplicating it in 3 dim
    num =  repmat(psi, [1 1 K]).*IR   ...
        + repmat(a1, [1 1 K]).* (IRR(y1,x1,:) - I0(y1,x1,:) + I0) ...
        + repmat(a2, [1 1 K]).*(IRR(y2,x2,:) - I0(y2,x2,:) + I0) ...
        + repmat(a3, [1 1 K]).*(IRR(y3,x3,:) - I0(y3,x3,:) + I0) ...
        + repmat(a4, [1 1 K]).*(IRR(y4,x4,:) - I0(y4,x4,:) + I0);
    %{    
    + repmat(a5, [1 1 K]).*(IRR(y0,x0,:) - I0(y0,x0,:));
        %}
    
    IRR = num./repmat(den + eps, [1 1 K]) .* (1-rho) + rho.*IRR;
    nbits= nbits-1;
end

end
