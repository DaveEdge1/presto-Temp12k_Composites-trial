function run_scc(fieldsMat, sandboxDir, seed)
% Run the published SCC_GMST_122719.m verbatim in a sandbox dir.
% fieldsMat : TS_fields.mat produced by ndjson_to_tsmat.py (per-field cells)
% sandboxDir: output dir; TS.mat is assembled here, inputs copied in,
%             script executed with cwd = sandbox (its save calls land here)
% seed      : rng seed (MATLAB batch sessions otherwise all start identical)
repoT12k = fullfile(getenv('HOME'), 'GitHub', 'Temperature12k');
plotcode = fullfile(repoT12k, 'ScientificDataDescriptor', 'PlottingCode');
sccdir   = fullfile(repoT12k, 'ScientificDataAnalysis', 'SCC');
mlscripts = fullfile(getenv('HOME'), 'Dropbox', 'ml_scripts');

if ~exist(sandboxDir, 'dir'); mkdir(sandboxDir); end

% helpers: PlottingCode first (correct find_nearest/gridMat/bin_x/lldistkm),
% ml_scripts as fallback for anything else
addpath(plotcode);
addpath(mlscripts);

% ---- assemble TS.mat (struct array from per-field cell arrays) ----
S = load(fieldsMat);
TS = struct( ...
  'dataSetName', S.f_dataSetName, ...
  'geo_latitude', S.f_geo_latitude, ...
  'geo_meanLat', S.f_geo_meanLat, ...
  'geo_meanLon', S.f_geo_meanLon, ...
  'paleoData_units', S.f_paleoData_units, ...
  'paleoData_inCompilation', S.f_paleoData_inCompilation, ...
  'interpretation1_seasonalityGeneral', S.f_interpretation1_seasonalityGeneral, ...
  'paleoData_temperature12kUncertainty', S.f_paleoData_temperature12kUncertainty, ...
  'age', S.f_age, ...
  'paleoData_values', S.f_paleoData_values); %#ok<NASGU>
save(fullfile(sandboxDir, 'TS.mat'), 'TS', '-v7');

% ---- copy the script's load-from-cwd inputs ----
copyfile(fullfile(plotcode, 'grid.mat'), sandboxDir);
copyfile(fullfile(plotcode, 'PAGES_multiMethodMeadian.txt'), sandboxDir);

% ---- run verbatim script in sandbox ----
oldd = cd(sandboxDir);
cleaner = onCleanup(@() cd(oldd));
rng(seed);   % [environmental patch] batch MATLAB otherwise always seeds identically
try
  run(fullfile(sccdir, 'SCC_GMST_122719.m'));
catch err
  % plotting tail may fail headless; outputs are saved before the plots
  fprintf('[run_scc] script ended with error (often just plotting): %s\n', err.message);
end
if exist(fullfile(sandboxDir, 'globalComp.mat'), 'file')
  fprintf('[run_scc] OK: globalComp.mat written in %s\n', sandboxDir);
else
  error('[run_scc] FAILED: globalComp.mat not produced');
end
end
