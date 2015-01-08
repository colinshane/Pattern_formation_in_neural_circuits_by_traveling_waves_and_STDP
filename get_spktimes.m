function times = get_spktimes(filename,varargin)
%
% times = get_spktimes(filename)
%
% Load all spike times from the .nd file called filename into a cell array.
% Assumes that the input layer is square. If it is not square, see optional
% inputs below. To specifiy the time period over which spike times should
% be read in, additional input arguments are needed.
%
% Inputs:
%  filename - .nd file name, including path
% 
%    Additional but optional inputs:
%   
%     * times = get_spktimes(filename,[NY NX 1]);
%    Specifies the size of the input layer, with NY units along the vertical
%    axis and NX units along the horizontal axis.
%     * times = get_spktimes(filename,[NY NX 1 NT]);
%    Load the first NT milliseconds. The input layer size must be specified in
%    order to specify NT.
%     * times = get_spktimes(filename,[NY NX 1 NT NS]);
%    Spike times are loaded from NS milliseconds into the simulation until NT
%    milliseconds into the simulation.
%
% Outputs
%  times - NY-by-NX cell array. Each cell contains a N-by-1 vector of N
%          spike times (in units of ms).


fid = fopen(filename);
fread(fid,1,'*int'); % Each .nd file begins with 1
nc = fread(fid,1,'*int'); %# chars next
fread(fid,nc,'*char'); % class
fread(fid,3,'*int'); % nconst, nvar, ntable
ntrials = fread(fid,1,'*int'); %ntrials
tn = 0;
flag = -1; % when max spike time > NT, stop loading spike times
j = 0;
while flag<0
  j = j + 1; % trial number
  %     if mod(j,20)==1
  %       fprintf('Loading spike times from trial: %d of %d\n',j,ntrials);
  %     end;
  fread(fid,1,'*int'); % tcode
  tref = fread(fid,1,'*int'); % tref
  fread(fid,1,'*int'); % nparam
  nunit = fread(fid,1,'*int'); % # units (should be a square)
  if numel(varargin)==0
    NY = sqrt(single(nunit));
    NX = NY;
    if mod(NY,1)~=0
      fclose(fid);
      error('*** GET_SPKTIMES: network is not square: size must be defined');
    end;
  elseif numel(varargin{1})<=2
    fclose(fid);
    error('*** GET_SPKTIMES: define size of network in 3 dimensions, [NY x NX x NPOP]');
  elseif numel(varargin{1})>2
    NY = varargin{1}(1);
    NX = varargin{1}(2);
    NZ = varargin{1}(3);
  end;
  if numel(varargin{1})>3
    NT = varargin{1}(4); % maximum spike time allowed (raw time units)
  else
    NT = inf;
  end;
  if numel(varargin{1})>4
    NS = varargin{1}(5); % Use spike times after this time (raw time units)
  else
    NS = 0;
  end;
  if nunit~=(NY*NX*NZ)
    fclose(fid);
    error('*** GET_SPKTIMES: nunit does not match grid size');
  end;
  if numel(varargin)>1
    fclose(fid);
    error('*** GET_SPKTIMES: define network/time dimensions in a vector: "[...]"');
  end;
  if j==1
    times = cell(NY,NX,NZ);
    nspk = cell(NY,NX,NZ);
    
    for k=1:NY
      for l=1:NX
        for m=1:NZ
          nspk{k,l,m}=0;
        end;
      end;
    end;
  end;
  if feof(fid)
    fprintf('Warning: reached end of file without exceeding \n');
    fprintf('designated max spike time. max(t) = %d\n',mt);
    %       fclose(fid);
    break;
  end;
  for l=1:NX
    for k=1:NY
      for m=1:NZ
        fread(fid,1,'*int'); % rtype
        nc = fread(fid,1,'*int'); % # characters in name
        fread(fid,nc,'*char'); % name
        fread(fid,1,'*int'); % rcode
        samp = fread(fid,1,'*float'); % sampling rate (bins/sec)
        t0 = fread(fid,1,'*int'); % t0 - trial start time (# bins)
        tn = fread(fid,1,'*int'); % tn - trial duration (# bins)
        n = fread(fid,1,'*int'); % n - number of spike times in trial
        nspk{k,l,m}=[nspk{k,l,m};n];
        if n>0
          t = tref + t0 + fread(fid,n,'*int'); % spike times
          if numel(t)>0
            mt = max(t);
          end;
          if max(t)>NT
            flag = 1;
          end;
          t = t(t>NS);
          times{k,l,m} = [times{k,l,m}; t(t<=NT)];
          if (sum(t > tref+t0+tn+(NZ-1)*5000)>0)
            fclose(fid);
            error('*** GET_SPKTIMES: spike time too large');
          end;
        end;
      end;
    end;
  end;
end;
fclose(fid);
% Clean out repeated spikes
for k=1:NY
  for l=1:NX
    for m=1:NZ
      d = diff(times{k,l,m});
      times{k,l,m}(d==0)=[];
      %       times{k,l}(times{k,l}>NT)=[];
    end;
  end;
end;
if NZ>1
  for j=2:NZ
    varargout{j-1} = times(:,:,j);
  end;
  times(:,:,2:end) = [];
  varargout{NZ} = nspk;
else
  varargout{1} = nspk;
end;
