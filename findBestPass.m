function pass = findBestPass(cfgBase, station, opt)
%==========================================================================
% findBestPass  Search a propagation window for visible passes over a ground
%               station and return the one reaching the highest elevation.
%
%   Given:
%     cfgBase  generateTrack config WITHOUT time fields (mode + tle/elements)
%     station  struct .latDeg .lonDeg .altKm
%     opt      struct (optional):
%                .searchHours  window length from epoch  (default 12)
%                .coarseDtSec  coarse search step        (default 15)
%                .fineDtSec    output sampling step       (default 1)
%                .maskDeg      horizon elevation mask     (default 10)
%
%   Returned struct pass:
%     .found        logical
%     .tSec         1xM  fine sample times over the best pass [s from epoch]
%     .tStartSec .tEndSec   pass bounds (mask crossings) [s]
%     .peakAltDeg   maximum elevation reached [deg]
%     .durationSec  pass duration [s]
%==========================================================================

    if nargin < 3 || isempty(opt), opt = struct(); end
    if ~isfield(opt,'searchHours'), opt.searchHours = 12;  end
    if ~isfield(opt,'coarseDtSec'), opt.coarseDtSec = 15;  end
    if ~isfield(opt,'fineDtSec'),   opt.fineDtSec   = 1;   end
    if ~isfield(opt,'maskDeg'),     opt.maskDeg     = 10;  end

    if ~isfield(station,'rEcef') || isempty(station.rEcef)
        station.rEcef = geodetic2ecef(station.latDeg, station.lonDeg, station.altKm);
    end

    % ---- coarse elevation profile ------------------------------------
    cfg = cfgBase;
    cfg.tSpanSec = [0, opt.searchHours*3600];
    cfg.dtSec    = opt.coarseDtSec;
    trk = generateTrack(cfg);

    alt = zeros(1, numel(trk.tSec));
    for k = 1:numel(trk.tSec)
        [~, a] = computeAzAlt(trk.rEci(:,k), gmstRad(trk.jd(k)), station);
        alt(k) = rad2deg(a);
    end

    % ---- find contiguous above-mask segments -------------------------
    above = alt > opt.maskDeg;
    d = diff([false, above, false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;

    pass.found = false;
    if isempty(starts)
        pass.tSec = []; pass.peakAltDeg = max(alt); pass.durationSec = 0;
        pass.tStartSec = NaN; pass.tEndSec = NaN;
        return;
    end

    % ---- pick the highest-culmination segment ------------------------
    peak = zeros(1, numel(starts));
    for j = 1:numel(starts)
        peak(j) = max(alt(starts(j):ends(j)));
    end
    [peakAlt, jb] = max(peak);

    tStart = trk.tSec(starts(jb));
    tEnd   = trk.tSec(ends(jb));
    % pad by one coarse step so the mask crossings are captured
    tStart = max(0, tStart - opt.coarseDtSec);
    tEnd   = tEnd + opt.coarseDtSec;

    pass.found       = true;
    pass.tSec        = tStart:opt.fineDtSec:tEnd;
    pass.tStartSec   = tStart;
    pass.tEndSec     = tEnd;
    pass.peakAltDeg  = peakAlt;
    pass.durationSec = tEnd - tStart;
end
