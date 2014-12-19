function w = compute_kappa(lr_only,nit,N,sres,btype,b_dur,iwi,lr_type,amp_pot,amp_dep,tau_pot,tau_dep,v,vsd,fr,tau_r,tau_d,zscale,w_type,sd)
%
% [w1 wo] = conv_w(lr_only,nit,N,OO,olag,L,llag,btype,b_dur,iwi,lr_type,amp_pot,amp_dep,tau_pot,tau_dep,v,vsd,fr,tau_r,tau_d,zscale,w_type,sd)
%
% Inputs:
% lr_only - 1: return STDP functions and FTs only; 0: run numerical integration
%     nit - number of iterations
%       N - number of presynaptic cells
%   btype - burst type, e.g. 'boxcar', 'rand_spks' etc. (see line-419 for full list)
%   b_dur - burst duration (s)
%     iwi - inter wave interval (s)
% lr_type - learning rule type (sym, asym)
% amp_pot - amplitude of potentiation
% amp_dep - amplitude of depression
% tau_pot - time scale for potentiation window
% tau_dep - time scale for depression window
%       v - wave-speed
%     vsd - SD in wave-speed
%      fr - presynaptic firing rate
%   tau_r - EPSP rise time
%   tau_d - EPSP decay time
%  zscale - learning rate
%  w_type - initial weights type (gaussian, uniform, box, random)
%      sd - width of inital weights (for gaussian & box) (units of cell spacings)
%
% Outputs:
% w - final synaptic weights

% Position of cells in space (0 corresponds to middle of input layer)
x = 0:sres:floor(N/2)*sres;
n = length(x);

% Compute kappa
z = get_lr(x,n,sres,btype,b_dur,iwi,lr_type,amp_pot,amp_dep,tau_pot,tau_dep,v,vsd,fr,tau_r,tau_d);
if lr_only % return kappa, comprising functions and FTs
  w = z;
  return;
end;
z1 = z.z1; z2 = z.z2; dv = z.dv;

% Initialise synaptic weights
if strcmp(w_type,'gaussian')==1
  w = gaussian(2*n-1,sd)';
  w = w / max(w);
  w = circshift(w,[0 -(n-1)]);
elseif strcmp(w_type,'uniform')==1
  w = 0.5 * ones(1,2*n-1) - 0.01 * randn(1,2*n-1); % ON
elseif strcmp(w_type,'box')==1  
  w = zeros(1,2*n-1);
  rad=sd;
  del = 0;
  w(ceil((2*n-1)/2)-sd+del:ceil((2*n-1)/2)+sd+del) =  1 - 0.01*rand(1,2*sd+1);
  w = circshift(w,[0 -(n-1)]);
elseif strcmp(w_type,'rand')==1  
  w = rand(1,2*n-1);
  w = circshift(w,[0 -(n-1)]);
end;

% Limit wynaptic weights to the range [0,1]
w(w<0) = 0;
w(w>1) = 1;

fw = fft(w); % FT of the weights

% % % % Arbour function
% A = zeros(size(w));
% A = gaussian(2*n-1,sd)'; % Gaussian arbour
% A = A / max(A(:));
% % % 
% A(ceil(N/2)-sd*1.5:ceil(N/2)+sd*1.5) = 1; % Flat arbour
% A = circshift(A,[0 -(n-1)]);

% Run numerical integration
figure;
for j=1:nit  
  % Recompute kappa if wave speed is variable
  if vsd>0
    z = get_lr(x,n,sres,btype,b_dur,iwi,lr_type,amp_pot,amp_dep,tau_pot,tau_dep,v,vsd,fr,tau_r,tau_d);
    z1 = z.z1; z2 = z.z2; dv = z.dv;
  end;
  
  % Convolve weights with kappa, reversing the wave direction after each iteration
  if mod(j,2)==1
    dw = zscale * fw .* z1;
  elseif mod(j,2)==0
    dw = zscale * fw .* z2;
  end;
  fw = fw + dw;
  w = ifft(fw);
  
  % Apply bounds to weights (without arbour)
  w(w<0) = 0; 
  w(w>1) = 1;
  
  % Apply bounds to weights (with arbour)
%   w(w>A)=A(w>A);
%   w(w<0)=0;
  
  % FT of weights after bounds are applied
  fw = fft(w); 

  % For viewing synaptic weight evolution
  if mod(j,101)==1
    plot(circshift(w,[0 n-1]),'b');
    hold on;
    plot(circshift(ifft(dw),[0 n-1])/zscale,'g');    
    set(gca,'xlim',[1 2*n-1],'ylim',[-0.2 1.4]);
    hold off;
    pause(0.1);
  end;
end;

% Return final synaptic weights
w = circshift(w,[0 n-1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute kappa
function z = get_lr(x,n,sres,btype,b_dur,iwi,lr_type,amp_pot,amp_dep,tau_pot,tau_dep,v,vsd,fr,tau_r,tau_d)

dv = vsd * randn;
v = v + dv;
if v<0.05 % Ensure wave speed is not too slow
  v = 0.05;
end;

% Compute STDP rule
if strcmp(lr_type,'sym')
  lr = amp_pot*exp(-((x/v)/tau_pot).^2/2) - amp_dep*exp(-((x/v)/tau_dep).^2/2);
  lr = [fliplr(lr(2:end)) lr];
  lr = circshift(lr,[0 (n-1)]);
elseif strcmp(lr_type,'asym')
  lr1 = amp_pot*exp(-(x/v)/tau_pot);
  lr1(1) = 0;
  lr2 = amp_dep*exp(-(x/v)/tau_dep);
  lr = [fliplr(lr1) -lr2(2:end)];
  lr = circshift(lr,[0 n]);
else
  error('??? CONV_W: lr_type not recognised\n');
end;

% Compute alpha
a = zeros(1,2*n-1);
if strcmp(btype,'boxcar') % Boxcar function
  a(x/v<b_dur) = fr * sres / v;
elseif strcmp(btype,'reg_spks') % Burst of regular spike times
  a(1:1000:end) = 1;
elseif strcmp(btype,'rand_spks') % Burst of random spike times
  na = sum(x/v<b_dur);
  a(1:na) = rand(1,na);
  pr = 1 - (fr * b_dur / na); 
  if pr<0 pr = 0; end;
  a(a>pr) = 1; a(a<=pr) = 0;
%   if sum(a>0)==0 a(:) = 1e-10; end;
  a = a*fr*sres/v;
elseif strcmp(btype,'sine') % Sine wave
  a = sin((0:sres:(2*n-2)*sres)/b_dur/v*2*pi);
  a = (a-min(a(:)))*fr/1000;
elseif strcmp(btype,'biphasic') % biphasic LGN stimulus response
  a = (x/v/0.025).^3.*exp(-x/v/0.025).*(1/factorial(3) - ((x/v/0.025).^2)/factorial(5)).*exp(-x/v.^2/2/0.25^2);
  a = [a zeros(1,size(x,2)-1)];
elseif strcmp(btype,'alpha') % Alpha function
  a = ((0:sres:(n-1)*2*sres)/v) .* exp(-((0:sres:(n-1)*2*sres)/v)/b_dur); 
  a = a / max(a); 
  a=a*fr*sres/v;
end;

% Inter-wave interval
if ((2*n-1)*sres/(v*iwi)-1)>=1
  b = zeros(size(a));
  for j=0:((2*n-1)*sres/(v*iwi)-1)
    b = b + circshift(a,[0 floor(j*v*iwi/sres)]);
  end;
  a = b;
end;

% Compute alpha for other wave direction
a2           = fliplr(a);
a2           = circshift(a2,[0 1]);

% Compute epsilon
e = exp(-((0:sres:(2*n-2)*sres)/v)/tau_d) - exp(-((0:sres:(2*n-2)*sres)/v)/tau_r);
e = e / sum(e);

% Compute epsilon for other wave direction
e2 = fliplr(e);
e2 = circshift(e2,[0 1]);

% FT of STDP rule
flr1  = fft(lr);
flr2 = fft(circshift(fliplr(lr),[0 1]));
% FT of alpha
fa   = fft(a);
fa2  = fft(a2);
% FT of epsilon
fe   = fft(e);
fe2  = fft(e2);

% FT of kappa
z1 = (flr1 .* fa .* fa2 .* fe);

% Kappa in real space
zz1 = ifft(z1);
zz1 = circshift(zz1,[0 n-1]);

% FT of kappa for other wave direction
z2 = (flr2 .* fa .* fa2 .* fe2);

% Kappa in real space for other wave direction
zz2 = ifft(z2);
zz2 = circshift(zz2,[0 n-1]);


% Normalise kappa so that changes in parameters do not greatly change speed of convergence
nrml = sqrt(sum(abs(z1).^2)); % To normalise z1 power.
if nrml>0
  z1 = z1 / nrml;
  z2 = z2 / nrml;
end;

% For output if lr_only==1
z.z1 = z1; z.z2 = z2; z.flr1 = flr1; z.flr2 = flr2; z.fa = fa; z.fa2 = fa2;
z.fe = fe; z.fe2 = fe2; z.zz1 = zz1; z.zz2 = zz2; z.dv = dv; z.lr = lr;
