function [y, fs, info] = straight_process(varargin)

%STRAIGHT_PROCESS - Cached STRAIGHT process of files
%   [Y, FS] = STRAIGHT_PROCESS(FILENAME, DF0)
%       Processes FILENAME with STRAIGHTV40_006b to apply DF0, which
%       represents:
%       - If DF0 is a number, it is treated as the distance in semitones from the
%         geometric mean F0.
%       - If DF0 is a string of the form '100.1Hz', it is treated as the new
%         absolue mean F0 in Hertz.
%       - If DF0 is a string starting with '+' or '-' and ending with 'Hz',
%         it is treated as a linear offset in Hertz.
%
%   [Y, FS] = STRAIGHT_PROCESS(FILENAME, DF0, DVTL)
%       DVTL is specified as a shift in semitones. ATTENTION: this is VTL
%       change, i.e. if DVTL is positive, formant frequencies are shifted
%       down.
%
%   [Y, FS] = STRAIGHT_PROCESS(FILENAME, DF0, DVTL, DURATION)
%       DURATION specifies the duration modificator:
%       - If DURATION is a number, it specifies the new duration in
%         seconds.
%       - If DURATION is a string of the form '*1.5', the duration is
%         multiplied by the number following '*'.
%       - If DURATION is a string of the form '+1.5s' or '-1.5s', the duration is
%         extended by the given number of seconds.
%       - If DURATION is empty ([] or '') or NaN, the duration is unchanged.
%
%   [Y, FS] = STRAIGHT_PROCESS(..., METHOD)
%       METHOD is a string that is one of the following: 
%       - 'straight': then STRAIGHTV40_006b is used (Matlab).
%       - 'tandem': then baseTandemSTRAIGHTV009x is used (Matlab).
%       - 'straightlib': then straightLib 1.1.4 is used (binary).
%       - 'world': Morise's WORLD vocoder (Matlab).
%
%   [Y, FS] = STRAIGHT_PROCESS(..., PARAMS)
%       Where PARAMS is a struct() giving extra parameters:
%       - cache_folder: the path where the cached files are stored
%       - cache_format: 'flac' (default) or 'wav'
%       - straight_path, tandem_path, straightlib_path: the paths to the
%         various STRAIGHT libraries.
%       - lower_f0_bound, upper_f0_bound: these are the min and max F0
%         values used during the F0 extraction process.
%
%   [Y, FS] = STRAIGHT_PROCESS(X, FS, ...)
%       Same as above but passing the signal directly instead of a file
%       name (no cache is created).
%
%   Important note on clipping: The returned signals are scaled such that
%   their RMS is identical to the original sound. That may mean that their
%   peak amplitude is larger than 1, which would result in clipping. The
%   cache files are stored in a way that prevents clipping. It is the
%   user's job to make sure that clipping does not occur at playback.

%--------------------------------------------------------------------------
% Authors:
%   2018-05-17: Etienne Gaudrain <etienne.gaudrain@cnrs.fr>
%
% Copyright CNRS (FR), UMCG (NL) and the author(s).
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------

% TODO:
%   - Find a way to get analysis only, with possibility to get F0 only

args = struct();
args.do_cache = true;

%-- Parsing arguments

if nargin<1
    error('Try harder: You need at least 1 argument...');
end

% Signal or filename
if ischar(varargin{1})
    args.filename = varargin{1};
    nargin_offset = 1;
elseif isnumeric(varargin{1})
    if numel(varargin{1})<2048
        error('The first argument must be either a filename or a vector containing more than 2048 elements');
    end
    args.x = varargin{1};
    if nargin<2 || ~isnumeric(varargin{2})
        error('If the first argument is a vector, the second argument has to be the sampling frequency');
    end
    args.fs = varargin{2};
    args.do_cache = false;
    nargin_offset = 2;
else
    error('The first argument is either a filename or a vector. %s given.', class(varargin{1}));
end

% dF0
nargin_offset = nargin_offset+1;
if nargin<nargin_offset
    args.dF0_type = 'semitones';
    args.dF0 = 0;
elseif isnumeric(varargin{nargin_offset})
    args.dF0_type = 'semitones';
    args.dF0 = varargin{nargin_offset};
elseif ischar(varargin{nargin_offset})
    dF0 = varargin{nargin_offset};
    if endswith(dF0, 'Hz')
        dF0 = dF0(1:end-2);
        dF0_num = str2double(dF0);
        if isnan(dF0_num)
            error('If DF0 is a string ending with "Hz", the part preceding the unit must be numeric (''%s'' given).', dF0);
        end
        if startswith(dF0, '+') || startswith(dF0, '-')
            args.dF0_type = 'linear_shift';
            args.dF0 = dF0_num;
        else
            args.dF0_type = 'absolute';
            args.dF0 = dF0_num;
        end
    else
        error('DF0 must be either numeric or a string ending with "Hz" (''%s'' given)', dF0);
    end
else
    error('The type of DF0 is not recognised (%s).', class(varargin{nargin_offset}));
end


% dVTL
nargin_offset = nargin_offset+1;
if nargin<nargin_offset
    args.dVTL = 0;
elseif isnumeric(varargin{nargin_offset})
    args.dVTL = varargin{nargin_offset};
else
    error('DVTL must be numeric (%s given).', class(varargin{nargin_offset}));
end

% duration
nargin_offset = nargin_offset+1;
if nargin<nargin_offset || isempty(varargin{nargin_offset}) || any(isnan(varargin{nargin_offset}))
    args.duration = 1.0;
    args.duration_type = 'factor';
elseif isnumeric(varargin{nargin_offset})
    args.duration = varargin{nargin_offset};
    args.duration_type = 'absolute';
elseif ischar(varargin{nargin_offset})
    dur = varargin{nargin_offset};
    if endswith(dur, 's') && (startswith(dur, '+') || startswith(dur, '-'))
        args.duration = str2double(dur(1:end-1));
        args.duration_type = 'offset';
    elseif startswith(dur, '*')
        args.duration = str2double(dur(2:end));
        args.duration_type = 'factor';
    else
        error('If DURATION is a string, it must either start with ''*'' or end with ''s''.');
    end
    if isnan(args.duration)
        error('DURATION could not be converted into a number ("%s" given).', dur);
    end
else
    error('DURATION must be numeric or string (%s given).', class(varargin{nargin_offset}));
end

% method & extra parameters
args.params = struct();
args.params.cache_folder = './straight_process_cache';
args.params.cache_format = 'flac';
args.params.straight_path = './straight';
args.params.tandem_path   = './tandem';
args.params.straightlib_path = './straightlib';
args.params.world_path = './world';
for i=1:2
    nargin_offset = nargin_offset+1;
    if nargin<nargin_offset
        if ~isfield(args, 'method')
            args.method = 'straight';
        end
    else
        if ischar(varargin{nargin_offset})
            args.method = lower(varargin{nargin_offset});
            switch args.method
                case {'straight', 'straightlib', 'tandem', 'world'}
                otherwise
                    error('METHOD is not valid ("%s" given).', args.method);
            end
        elseif isstruct(varargin{nargin_offset})
            params = varargin{nargin_offset};
            f = fieldnames(params);
            for k=1:length(f)
                args.params.(f{k}) = params.(f{k});
            end
        else
            error('The last argument must be a string (METHOD) or a struct (PARAMS). (%s given)', class(varargin{nargin_offset}));
        end
    end
end

if strcmp(args.method, 'straightlib') && ~strcmp(args.dF0_type, 'semitones')
    error('DF0 can only be a factor (semitones) when straightlib is selected.');
end
if strcmp(args.method, 'straightlib') && ( (args.duration ~= 1.0) || ~strcmp(args.duration_type, 'factor'))
    error('DURATION cannot be anything else than 1.0 when straightlib is selected.');
end

%--- Check if the file exists in cache (if do_cache==1)

info = struct();

[cache_filename_snd, cache_filename_mat] = make_fname(args);

if args.do_cache && exist(cache_filename_snd, 'file')
    info.analysis_time = 0;
    tic();
    [y, fs] = audioread(cache_filename_snd);
    snd_inf = audioinfo(cache_filename_snd);
    info.scaling = NaN;
    if isfield(snd_inf, 'Comment')
        scaling = sscanf(snd_inf.Comment, 'scaling:%f');
        if isempty(scaling)
            warning('Scaling couldn''t be properly retrieved from cached file %s. A scaling of 1 was applied, but that may mean that the RMS is not correct.', cache_filename_snd);
        else
            info.scaling = scaling;
        end
    end
    if info.scaling~=1 && ~isnan(info.scaling)
        y = y*info.scaling;
    end
    info.source = 'cache';
    info.source_file = cache_filename_snd;
    info.synthesis_time = toc();
    return
end

%-- Dispatch based on method

args.cache_filename_snd = cache_filename_snd;
args.cache_filename_mat = cache_filename_mat;

if isfield(args, 'filename')
    [args.x, args.fs] = audioread(args.filename);
end

switch args.method
    case 'straight'
        [y, fs, info] = process_straight(args);
    case 'straightlib'
        [y, fs, info] = process_straightlib(args);
    case 'tandem'
        [y, fs, info] = process_tandem(args);
    case 'world'
        [y, fs, info] = process_world(args);
end

m = max(abs(y));
if m>.98
    scaling = .98/m;
else
    scaling = 1;
end

y = y*scaling;
info.scaling = scaling;

if args.do_cache
    audiowrite(args.cache_filename_snd, y, fs, 'Comment', sprintf('scaling:%10f', 1/scaling));
    % Note: the scaling we write in the file is the scaling that will have to
    % be applied to restore the proper RMS
end

%============================================================================================================================

function [f0, ser, tar] = apply_args(args, f0, duration)

% Applies the modifications specified in ARGS to F0 and calculates the SER,
% and the timeAxisRatio.

%-- F0
switch args.dF0_type
    case {'semitones', 's'}
        f0 = f0 * 2.^(args.dF0/12);
    case {'linear_shift', 'l'}
        f0(f0~=0) = f0(f0~=0) + args.dF0;
    case {'absolute', 'a'}
        mf0 = exp(mean(log(f0(f0~=0))));
        f0 = f0 / mf0 * args.dF0;
    otherwise
        error('dF0_type not recognised ("%s"),', args.dF0_type);
end

%-- SER
ser = 2.^(-args.dVTL/12);

%-- TAR
switch args.duration_type
    case {'factor', 'f'}
        tar = args.duration;
    case {'absolute', 'a'}
        tar = args.duration / duration;
    case {'offset', 'o'}
        tar = args.duration + duration;
    otherwise
        error('duration_type not recognised ("%s"),', args.dF0_type);
end

%============================================================================================================================

function [y, fs, info] = process_straight(args)

% Effect the transformation with STRAIGHT

addpath(args.params.straight_path);

info = struct();
tic();

if args.do_cache && exist(args.cache_filename_mat, 'file')
    load(args.cache_filename_mat);
    info.mat_filename = args.cache_filename_mat;
    info.mat_source = 'cache';
    info.analysis_time = toc();
else
    fs = args.fs;
    prms = struct();
    if isfield(args.params, 'lower_f0_bound')
        prms.F0searchLowerBound = args.params.lower_f0_bound;
    end
    if isfield(args.params, 'upper_f0_bound')
        prms.F0searchUpperBound = args.params.upper_f0_bound;
    end
    [f0, ap] = exstraightsource(args.x, fs, prms);
    spenv = exstraightspec(args.x, f0, fs);
    rms_x = rms(args.x);
    duration = length(args.x)/fs;
    if args.do_cache
        filename = args.filename;
        if ~exist(args.params.cache_folder, 'dir')
            mkdir(args.params.cache_folder);
        end
        save(args.cache_filename_mat, 'fs', 'spenv', 'f0', 'ap', 'rms_x', 'filename', 'duration');
    end
    info.mat_source = 'generated';
    info.analysis_time = toc();
end

info.f0 = f0;
info.f0_temporal_positions = (0:length(f0)-1)/1e3;

tic();

[f0, ser, tar] = apply_args(args, f0, duration);

p = struct();
p.frequencyAxisMappingTable = ser;
p.timeAxisStretchingFactor = tar;

y = exstraightsynth(f0, spenv, ap, fs, p);
y = y / rms(y) * rms_x;

info.synthesis_time = toc();

rmpath(args.params.straight_path);

%============================================================================================================================

function [y, fs, info] = process_tandem(args)

% Effect the transformation with TANDEM

addpath(args.params.tandem_path);

info = struct();
tic();

if args.do_cache && exist(args.cache_filename_mat, 'file')
    load(args.cache_filename_mat);
    info.mat_filename = args.cache_filename_mat;
    info.mat_source = 'cache';
    info.analysis_time = toc();
else
    fs = args.fs;
    
    % Source information extraction
    optP = struct();
    optP.channelsPerOctave = 3;

    if isfield(args.params, 'lower_f0_bound')
        optP.f0floor = args.params.lower_f0_bound;
    end
    if isfield(args.params, 'upper_f0_bound')
        optP.f0ceil = args.params.upper_f0_bound;
    end
    
    r  = exF0candidatesTSTRAIGHTGB(args.x, fs, optP);
    rc = autoF0Tracking(r, args.x);
    q = aperiodicityRatio(args.x, rc, 1);
    
    % Filter information extraction
    prmIn = struct();
    prmIn.exponentControl = 0.25;
    prmIn.compensationCoefficient = -0.3;
    f = exSpectrumTSTRAIGHTGB(args.x, fs, rc,prmIn);
    f.frequencies = (0:size(f.spectrogramSTRAIGHT, 1)-1) / (size(f.spectrogramSTRAIGHT, 1)-1) * f.samplingFrequency/2;
    
    rms_x = rms(args.x);
    duration = length(args.x)/fs;
    if args.do_cache
        filename = args.filename;
        if ~exist(args.params.cache_folder, 'dir')
            mkdir(args.params.cache_folder);
        end
        save(args.cache_filename_mat, 'fs', 'q', 'f', 'rms_x', 'filename', 'duration');
    end
    info.mat_source = 'generated';
    info.analysis_time = toc();
end

info.f0 = q.f0;
info.f0_temporal_positions = q.temporalPositions;

tic();

[q.f0, ser, tar] = apply_args(args, q.f0, duration);

% TODO: check if AP needs to be shifted as well

f.spectrogramSTRAIGHT = interp2(f.temporalPositions, f.frequencies', f.spectrogramSTRAIGHT, f.temporalPositions, f.frequencies' / ser, 'cubic', 0);
f.temporalPositions = f.temporalPositions * tar;
q.temporalPositions = q.temporalPositions * tar;

s = exTandemSTRAIGHTsynthNx(q, f);
y = s.synthesisOut;
y = y / rms(y) * rms_x;

info.synthesis_time = toc();

rmpath(args.params.tandem_path);


%============================================================================================================================

function [y, fs, info] = process_straightlib(args)

% Effect the transformation with StraightLib

info = struct();

cmd = fullfile(args.params.straightlib_path, 'straight');
if ispc()
    cmd = [cmd, '.exe'];
end

if args.do_cache
    info.analysis_time = 0;
    if ~exist(args.cache_filename_mat, 'file')
        tic();
        cmd_args = sprintf('-ana -f0infofile "%s" -apfile "%s" "%s" "%s"', args.cache_filename_mat, args.cache_filename_mat, args.filename, args.cache_filename_mat);
        [s, ~] = system([cmd, ' ', cmd_args]);
        if s~=0
            error('There was an error during the execution of straightlib...');
        end
        info.analysis_time = toc();
    end

    tic();
    [~, ser, tar] = apply_args(args, 1, 1);
    
    cmd_args = sprintf('-syn -pc %f -fc %f "%s" "%s"', 2.^(args.dF0/12), ser, args.cache_filename_mat, args.cache_filename_snd);
    
    if isfield(args.params, 'lower_f0_bound')
        cmd_args = ['-lf0 ', num2str(args.params.lower_f0_bound), ' ', cmd_args];
    end
    if isfield(args.params, 'upper_f0_bound')
        cmd_args = ['-uf0 ', num2str(args.params.upper_f0_bound), ' ', cmd_args];
    end
    
    [s, ~] = system([cmd, ' ', cmd_args]);
    if s~=0
        error('There was an error during the execution of straightlib...');
    end
    info.synthesis_time = toc();
end

[y, fs] = audioread(args.cache_filename_snd);

%============================================================================================================================

function [y, fs, info] = process_world(args)

% Effect the transformation with TANDEM

addpath(args.params.world_path);

info = struct();
tic();

if args.do_cache && exist(args.cache_filename_mat, 'file')
    load(args.cache_filename_mat);
    info.mat_filename = args.cache_filename_mat;
    info.mat_source = 'cache';
    info.analysis_time = toc();
else
    fs = args.fs;
    
    opt = struct();
    if isfield(args.params, 'lower_f0_bound')
        opt.f0_floor = args.params.lower_f0_bound;
    end
    if isfield(args.params, 'upper_f0_bound')
        opt.f0_ceil = args.params.upper_f0_bound;
    end
    
    % Source information extraction
    f0p = Harvest(args.x, fs, opt);
    spe = CheapTrick(args.x, fs, f0p);
    spe.frequencies =  (0:size(spe.spectrogram, 1)-1) / (size(spe.spectrogram, 1)-1) * spe.fs/2;
    src = D4C(args.x, fs, f0p);
    
    rms_x = rms(args.x);
    duration = length(args.x)/fs;
    if args.do_cache
        filename = args.filename;
        if ~exist(args.params.cache_folder, 'dir')
            mkdir(args.params.cache_folder);
        end
        save(args.cache_filename_mat, 'fs', 'spe', 'src', 'rms_x', 'filename', 'duration');
    end
    info.mat_source = 'generated';
    info.analysis_time = toc();
end

info.f0 = src.f0;
info.f0_temporal_positions = src.temporal_positions;

tic();

[src.f0, ser, tar] = apply_args(args, src.f0, duration);

spe.spectrogram = interp2(spe.temporal_positions, spe.frequencies', spe.spectrogram, spe.temporal_positions, spe.frequencies' / ser, 'cubic', 0);
src.aperiodicity = interp2(src.temporal_positions, spe.frequencies', src.aperiodicity, src.temporal_positions, spe.frequencies' / ser, 'cubic', 0);
spe.temporal_positions = spe.temporal_positions * tar;
src.temporal_positions = src.temporal_positions * tar;

y = Synthesis(src, spe);
y = y / rms(y) * rms_x;

info.synthesis_time = toc();

rmpath(args.params.world_path);

%============================================================================================================================
% Some utility functions
%-----------------------------------------------------------------------------------------

function [snd_fname, mat_fname] = make_fname(args)

% Computes the name of the cached file.

if isfield(args, 'filename')
    [p, name, ext] = fileparts(args.filename);
    ext = lower(ext);
    m = md5(p);
    s = '';
    s = [s, sprintf('f%.2f%s', args.dF0, args.dF0_type(1))];
    s = [s, sprintf('_v%.2f', args.dVTL)];
    s = [s, sprintf('_d%.2f%s', args.duration, args.duration_type(1))];
    s = [s, sprintf('.%s', args.method)];
    snd_fname = fullfile(args.params.cache_folder, sprintf('%s._.%s%s._.%s.%s', m(1:8), name, ext, s, args.params.cache_format));
    if strcmp(args.method, 'straightlib')
        ext_dat = '.stf';
    else
        ext_dat = '.mat';
    end
    mat_fname = fullfile(args.params.cache_folder, sprintf('%s._.%s%s._.%s%s', m(1:8), name, ext, args.method, ext_dat));
else
    snd_fname = '';
    mat_fname = '';
end

%-----------------------------------------------------------------------------------------

function b = startswith(str, start)

%STARTSWITH(STR, START) is True if STR starts with START
%   If STR is shorter than START, it returns False. This function has been
%   implemented in newer versions of Matlab.

%--------------------------------------------------------------------------
% Etienne Gaudrain <etienne.gaudrain@crns.fr> - 2017-10-05
% CNRS UMR5292, CRNL, Lyon, FR - RuG, UMCG, KNO, Groningen, NL
%--------------------------------------------------------------------------

if length(str)<length(start)
    b = false;
elseif isempty(start)
    b = true;
else
    b = strcmp(str(1:length(start)), start);
end

%-----------------------------------------------------------------------------------------

function b = endswith(str, sub)

%ENDSWITH(STR, SUB) is True if STR ends with SUB
%   If STR is shorter than SUB, it returns False. This function has been
%   implemented in newer versions of Matlab.

%--------------------------------------------------------------------------
% Etienne Gaudrain <etienne.gaudrain@crns.fr> - 2018-04-30
% CNRS UMR5292, CRNL, Lyon, FR - RuG, UMCG, KNO, Groningen, NL
%--------------------------------------------------------------------------

if length(str)<length(sub)
    b = false;
elseif isempty(sub)
    b = true;
else
    b = strcmp(str(end-length(sub)+1:end), sub);
end

%-----------------------------------------------------------------------------------------
function hash = md5(msg)

md = java.security.MessageDigest.getInstance('MD5');
hash = typecast(md.digest(uint8(msg)), 'uint8');
hash = lower(reshape(dec2hex(hash)', 1, []));
