function [y, unstable, a, u, Pf, n, w,eta_0,eta_U,eta_H] = MISO_sim(x,c,Bff,Bfb,Qff,Qfb,Optff,Optfb,seed,varargin)
%MISO_sim simulates a MISO simulation given input spike timing,
%coefficients, basis functions
%
%Input Arguments:
%   x, T x N binary input spiking vector
%   c, either normalized or not:
%       If normalized, contains fields (sig,c.k, c.h) 
%       If un-normalized, contains fields (c_0,c.k, c.h) 
%   Bff, Bfb: feedforward and feedback basis function objects
%   Qff, Qfb, feedforward and feedback model order
%   Optff, Optfb, options for creating design matrix, default values should
%       be empty,
%   seed, rng seed for reproducability
%   Name Value Pairs:
%       'unstable_a_thresh', default value is 1, used for detection of
%           instability
%       'unstable_timeout', default value is empty, if non-empty, and the
%           system becomes unstable, output spiking will be prevented for
%           unstable_timeout time bins to regain stability
%
%Output Arguments:
%   y, T x 1 vector of output sike timing
%   unstable, 1 or 0, whether or not the simulation was unstable
%   a, V_a*c.h  (feedback model component, affected by normalization)
%   u, V_u*c.k (feedforward model component, affected by normalization)
%   Pf, T x 1 firing probability at each time step
%   n, T x 1, noise values used in simulation
%   w, T x 1, combined feedforward and feedback,  Pf = g^-1(w), (NOT
%       affected by normalization) (w = eta_0+ eta_U+eta_H ).  In GLM
%       literature w is equivalent to eta
%   eta_0, T x 1, offset component (NOT affected by normalization)
%   eta_U, T x 1, feedforward component (NOT affected by normalization)
%   eta_H, T x 1, feedback component (NOT affected by normalization)


[unstable_a_thresh,unstable_timeout,link,g] = process_options(varargin,'unstable_a_thresh',1,'unstable_timeout',[],'link','probit','g',[]);
unstable=false;
%normalized=1 is the standard case when eta(t) = (-1+V_u(t)*c.k+V_a(t)*c.h)/c.sig;
%normalized=0 is the case when eta(t) = c.c0+V_u(t)*c.k+V_a(t)*c.h;
if isfield(c,'sig')
    normalized=1;
else
    normalized=0;
end
%set seed
rng(seed);

%create U
T = size(x,1);
y = zeros(T,1);  %y vector creation
y_f = zeros(T,1);  %y vector creation
tstart = tic;
[V, ~, ff_inds, ~] = createDesign(x,y,Bff,Bfb,Qff,Qfb,Optff,Optfb);
tend = toc(tstart);
disp(['Initial design matrix made in ' num2str(tend) ' s'])
V_u = V(:,ff_inds);

% u = V_u*c.k;
if ~isempty(V_u) % Modified by Xiwei, Dec.14, 2021
    u = V_u*c.k;
else
    u = zeros(T, 1);
end
if ~isempty(g)
    u = u.*g';
end    

%initialize a and feedback
a = zeros(T,1);
B_h = Bfb.calcB();
L_h = Bfb.getN();  %this is the number of basis functions

%initialize other intermediate vectors
Pf = nan(T,1);
w = nan(T,1);
n = rand(T,1);
eta_0= nan(T,1);  %these are just for debugging for now
eta_U= nan(T,1);    %these are just for debugging for now
eta_H= nan(T,1); %these are just for debugging for now
%loop through time
tstart = tic;
t=1;
while t<T % Modified by Xiwei, previous :while t<=T
    %for t=1:T
    if normalized
        if numel(c.sig)==1  %this happens in the stationary case when there is only one baseline firing rate
            w(t) = (-1+u(t)+a(t))/c.sig;  %w is not officially w from journal papers, w here is equivalent to eta in GLM lit
            eta_0(t) = -1/c.sig; %these are just for debugging for now
            eta_U(t) = u(t)/c.sig; %these are just for debugging for now
            eta_H(t) = a(t)/c.sig; %these are just for debugging for now
        else
            w(t) = (-1+u(t)+a(t))/c.sig(t);
            eta_0(t) = -1/c.sig(t); %these are just for debugging for now
            eta_U(t) = u(t)/c.sig(t); %these are just for debugging for now
            eta_H(t) = a(t)/c.sig(t); %these are just for debugging for now
        end
    else
        if numel(c.c_0)==1
            w(t) = c.c_0+u(t)+a(t);  %this is the case when c0 is stationary
            eta_0(t) = c.c_0; %these are just for debugging for now
            eta_U(t) = u(t); %these are just for debugging for now
            eta_H(t) = a(t); %these are just for debugging for now
        else
            w(t) = c.c_0(t)+u(t)+a(t);  %this is the case when c0 is NONstationary
            eta_0(t) = c.c_0(t); %these are just for debugging for now
            eta_U(t) = u(t); %these are just for debugging for now
            eta_H(t) = a(t); %these are just for debugging for now
        end
    end
    switch link
        case 'probit'
            Pf(t) = normcdf(w(t));
        case 'logit'
            Pf(t) = exp(w(t))./(1+exp(w(t)));
    end
    if n(t) <Pf(t)
        y(t) = 1;
        y_f(t+1) = 1;
        if Qfb~=0
            V_a = X2V(y_f, Bfb, Qfb,Optfb);
            a = V_a*c.h;
            
            if a(t)>unstable_a_thresh   %this stops the simulation once the feedback component reaches above a certain threshold at a = 1, Pf should equal 0.5
                unstable=true;
                if isempty(unstable_timeout)
                    break
                else
                    t=t+unstable_timeout;
                end
            end
        end
    end
    if mod(t,50000)==0
        tend = toc(tstart);
        disp([ 'sim t=' num2str(t) ', Elapsed time is ' num2str(tend) ' s'])
        tstart =tic;
    end
    t=t+1;
end
toc(tstart)