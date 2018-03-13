%Generates Figures 3 and 4
clear;
eps = 0.10;
tau = 0.1;
d = 20;
lambda = 0.2;
subsampSize = 100;

%Process mappings from color names to numbers
rawColor = importdata('./Novembre_etal_2008_misc-master/files/colors.txt');
keys = rawColor.rowheaders;
vals = mat2cell(rawColor.data, ones(37,1));
colorMap = containers.Map(keys, vals);

%Process data points
A = importdata('./Novembre_etal_2008_misc-master/files/fixed.eigs');
eigs = importdata('./Novembre_etal_2008_misc-master/files/POPRES_08_24_01.EuroThinFinal.LD_0.8.exLD.out0-PCA.eval');
eigs = eigs(1:d);
data = zeros(length(A), d);
for i = 1:length(A)
    v = str2num(A{i,1}(1:265));
    data(i,:) = v(3:end).*eigs';
end

%Convert data point color names to numbers
colorsFid = fopen('./Novembre_etal_2008_misc-master/files/POPRESID_Color.txt');
colorsStrings = textscan(colorsFid, '%d %s');
colorsStrings = colorsStrings{2};
dataColors = zeros(length(colorsStrings),3);
for i = 1:length(colorsStrings)
    dataColors(i, :) = colorMap(colorsStrings{i})/255;
end

%Generate noisy points
N = length(A)/(1-eps);
noise = (1/24)*[randi([0, 2], round(eps*N), d/2) randi([2, 3], round(eps*N), d/2)];

%Generate noisy point colors
noiseColors = zeros(round(eps*N), 3);

%Combine data and noise,
randRot = orth(randn(d, d));
data = data * randRot;
noise = noise * randRot;
D = [data; noise];
C = [dataColors; noiseColors];
subsampD = datasample(D, subsampSize);

%Generate original data's projection
[dataU, ~, ~] = svd(cov(data));

%Generate noised data's projection
[noisedU, ~, ~] = svd(cov(D));

%Generate filter output and projection
[filterM, filtered, filteredMetadata] = filterGaussianCovTuned(D, C, eps, tau, 0);
[filterU, ~, ~] = svd(filterM);

%Generate filter output and projection (no noise)
[filterNoNoiseM, filteredNoNoise, filterNoNoiseMetadata] = filterGaussianCovTuned(data, dataColors, eps, tau, 0);
[filterNoNoiseU, ~, ~] = svd(filterNoNoiseM);

%Generate LRV projection
[~, lrvM, ~] = agnosticCovarianceGeneral(D, eps);
[lrvU, ~, ~] = svd(lrvM);

%Generate Pruning projection
[pruneM, ~, ~] = pruneGaussianCov(D, C, tau, 0);
[pruneU, ~, ~] = svd(pruneM);

%Generate RANSAC projection
[ransacM] = ransacMVE(D, eps, ceil(d / eps), 1000);
[ransacU, ~, ~] = svd(ransacM);

%Generate CLMW SDP projection
cvx_begin quiet
    variable clmwsdpL(subsampSize, d) 
    variable S(subsampSize, d) 
    dual variable Q
    minimize norm_nuc(clmwsdpL) + lambda * sum(sum(abs(S)))
    subject to 
        clmwsdpL + S == subsampD;
cvx_end
clmwsdpM = clmwsdpL' * clmwsdpL;
[clmwsdpU, ~, ~] = svd(clmwsdpM);

%Generate CLMW AM projection
[clmwamL, ~] = ADPCP(D, lambda);
clmwamM = clmwamL' * clmwamL;
[clmwamU, ~, ~] = svd(clmwamM);

%Generate XCS SDP projection
cvx_begin quiet
    variable xcsL(subsampSize, d) 
    variable S2(subsampSize, d) 
    dual variable Q
    minimize norm_nuc(xcsL) + lambda * sum(norms(S2, 2, 2))
    subject to 
        xcsL + S2 == subsampD;
cvx_end
xcsM = xcsL' * xcsL;
[xcsU, ~, ~] = svd(xcsM);

%Generate plot for original data
fig = figure(1);
clf;
scatter(data*dataU(:,1), data*dataU(:,2), [], dataColors);
hold on;
scatter(noise*dataU(:,1), noise*dataU(:,2), [], noiseColors);
axis([-0.25 0.35 -0.15 0.2])
title('Original Data')
set(fig, 'Position', [100, 100, 640, 480]);

%Generate plot for noised data
fig = figure(2);
clf
scatter(data*noisedU(:,1), data*noisedU(:,2), [], dataColors);
hold on;
scatter(noise*noisedU(:,1), noise*noisedU(:,2), [], noiseColors);
axis([-0.25 0.35 -0.15 0.2])
title('Empirical Noised Singular Vectors')
set(fig, 'Position', [100, 100, 640, 480]);

%Generate plot for filter's output
fig = figure(3);
clf
scatter(filtered*filterU(:,1),filtered*filterU(:,2), [], filteredMetadata);
axis([-0.25 0.35 -0.15 0.2])
title('Filter Output');
set(fig, 'Position', [100, 100, 640, 480]);

%Generate plot for filter's output without noise
fig = figure(4);
clf
scatter(filteredNoNoise*filterNoNoiseU(:,1), filteredNoNoise*filterNoNoiseU(:,2), [], filterNoNoiseMetadata);
axis([-0.25 0.35 -0.15 0.2])
title('Filter Output (no noise)');
set(fig, 'Position', [100, 100, 640, 480]);

%Generate plot for the filter's projection
fig = figure(5);
clf
scatter(data*filterU(:,1), data*filterU(:,2), [], dataColors);
hold on;
scatter(noise*filterU(:,1), noise*filterU(:,2), [], noiseColors);
axis([-0.25 0.35 -0.15 0.2])
title('Filter Projection')
set(fig, 'Position', [100, 100, 640, 480]);

%Generate plot for LRV's projection
fig = figure(6);
clf
scatter(data*lrvU(:,1), data*lrvU(:,2), [], dataColors);
hold on;
scatter(noise*lrvU(:,1), noise*lrvU(:,2), [], noiseColors);
axis([-0.25 0.35 -0.15 0.2])
title('LRV Projection')
set(fig, 'Position', [100, 100, 640, 480]);

%Generate plot for pruned data's projection
fig = figure(7);
clf
scatter(data*pruneU(:,1), data*pruneU(:,2), [], dataColors);
hold on;
scatter(noise*pruneU(:,1), noise*pruneU(:,2), [], noiseColors);
axis([-0.25 0.35 -0.15 0.2])
title('Pruning Projection')
set(fig, 'Position', [100, 100, 640, 480]);

%Generate plot for RANSAC's projection
fig = figure(8);
clf
scatter(data*ransacU(:,1), data*ransacU(:,2), [], dataColors);
hold on;
scatter(noise*ransacU(:,1), noise*ransacU(:,2), [], noiseColors);
axis([-0.25 0.35 -0.15 0.2])
title('RANSAC Projection')
set(fig, 'Position', [100, 100, 640, 480]);

%Generate plot for CLMW SDP projection
fig = figure(9);
clf
scatter(data*clmwsdpU(:,1), data*clmwsdpU(:,2), [], dataColors);
hold on;
scatter(noise*clmwsdpU(:,1), noise*clmwsdpU(:,2), [], noiseColors);
axis([-0.25 0.35 -0.15 0.2])
title('CLMW SDP Projection')
set(fig, 'Position', [100, 100, 640, 480]);

%Generate plot for CLMW AM projection
fig = figure(10);
clf
scatter(data*clmwamU(:,1), data*clmwamU(:,2), [], dataColors);
hold on;
scatter(noise*clmwamU(:,1), noise*clmwamU(:,2), [], noiseColors);
axis([-0.25 0.35 -0.15 0.2])
title('CLMW ADMM Projection')
set(fig, 'Position', [100, 100, 640, 480]);

%Generate plot for XCS projection
fig = figure(11);
clf
scatter(data*xcsU(:,1), data*xcsU(:,2), [], dataColors);
hold on;
scatter(noise*xcsU(:,1), noise*xcsU(:,2), [], noiseColors);
axis([-0.25 0.35 -0.15 0.2])
title('XCS Projection')
set(fig, 'Position', [100, 100, 640, 480]);