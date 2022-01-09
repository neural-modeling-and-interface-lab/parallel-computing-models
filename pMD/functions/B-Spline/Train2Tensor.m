function P = Train2Tensor(X, E, L)
% X: spike trains with 2ms resolution
% E: event timings in sec
% L: number of bins in pattern tensor

[Len, N] = size(X);
NumTrial = length(E);

P = zeros(NumTrial,L+1,N); % spike tensor % NumTrial - # of trials(Sample_response); L+1 - Length of bin window; N - # of Channels(CA3+CA1)

for i = 1:NumTrial
    t = round(E(i)*500);
    tmin = t - L/2;
    tmax = t + L/2;
    a = X(tmin:tmax,:);
    P(i,:,:) = a;
end
