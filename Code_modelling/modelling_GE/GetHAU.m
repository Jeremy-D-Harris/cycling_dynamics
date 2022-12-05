function HAU = GetHAU(W, D, params)

% HAU is total number of virions divided by the particles factor:
HAU = (W+D)/params.HAU_particles_factor;

% add measurement noise once tranformed data to log2 scale:
% log2_HAU = log2(HAU);
% log2_HAU_measured = normrnd(log2_HAU, params.HAU_noise);
% HAU_measured = 2.^log2_HAU_measured;

% figure;
%subplot(2,1,1); plot(log2_HAU, 'k'); hold on; plot(log2_HAU_measured, 'r');
%subplot(2,1,2); semilogy(HAU, 'k'); hold on; semilogy(HAU_measured, 'r'); 

% replace data points that fall below the detection limit to NaN:
locs = find(HAU <= params.LOD_HAU);
% HAU_measured(locs) = NaN;
HAU(locs) = params.LOD_HAU/2;