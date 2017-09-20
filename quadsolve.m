% function [P1, P2 ] = quadsolve(resp1, resp2, ntrials)
%
% This function carries out the quadratic analysis from Dugue et al (2015).
% resp1 and resp2 are vectors of boolean responses to the probes.
% P1 and P2 are the probabilty of being correct on the 'attended' and
% 'unattended' sides respectively.
%
% JCM update 1: Don't assign the greater value to P1, just assign the first.
% This should reduce the bias.

function [P1, P2] = quadsolve(resp1, resp2, ntrials)

Pboth = sum(resp1 == 1 & resp2 == 1)/ntrials;
Pnone = sum(resp1 == 0 & resp2 == 0)/ntrials;

b = 1 + Pboth - Pnone;
c = Pboth;

det = b^2 - 4*c;

tmp = [(b + sign(det)*sqrt(abs(det)))/2,(b - sign(det)*sqrt(abs(det)))/2];
P1 = tmp(1);
P2 = tmp(2);

%% Update 1
% tmp = [(b + sign(det)*sqrt(abs(det)))/2,(b - sign(det)*sqrt(abs(det)))/2];
% 
% % Set so that P1 is always the more attended
% P1 = max(tmp);
% P2 = min(tmp);