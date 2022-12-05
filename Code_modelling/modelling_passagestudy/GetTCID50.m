function TCID50 = GetTCID50(W, params, f)

% relationship between TCID50 titer (per ml) and the mean number of PFU/ml:
% TCID50*0.70 = PFU/ml, where the 0.70 is more exactly given by: -log(0.5)
% Going the other way:
TCID50 = f*W;  

% add measurement noise once tranformed data to log2 scale:
% log10_TCID50 = log10(TCID50);
% log10_TCID50_measured = normrnd(log10_TCID50,params.TCID50_noise);
% TCID50_measured = 10.^log10_TCID50_measured;

% figure;
%subplot(2,1,1); plot(log2_TCID50, 'k'); hold on; plot(log2_TCID50_measured, 'r');
%subplot(2,1,2); semilogy(TCID50, 'k'); hold on; semilogy(TCID50_measured, 'r'); 

% replace data points that fall below the detection limit to NaN:
% locs = find(TCID50 <= params.LOD_TCID50);
% TCID50_measured(locs) = NaN;
% TCID50(locs) = params.LOD_TCID50/2;