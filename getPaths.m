function paths = getPaths() 

paths.project     = fileparts(which('getPaths.m'));
paths.figures     = fullfile(paths.project, 'Figures'); 
paths.stats       = fullfile(paths.project, 'Stats'); 
paths.data        = fullfile(paths.project, 'Data');
paths.results     = fullfile(paths.project, 'Results');

if ~isdir(paths.figures), mkdir(paths.figures); end
if ~isdir(paths.stats),   mkdir(paths.stats); end
if ~isdir(paths.results),   mkdir(paths.results); end