function result = classifyMultiChannelArtifacts(dataA, dataB, varargin)
% Notes from Mai-Anh: 
% There are 3 possible artifacts present in the data:
% 1) random spikes upward or downward: sharp, extreme value
% 2) in UG recordings (before 2024/03/14), there was an issues 
% where the IR light source kept overheating and as a result would blink. 
% Usually, this light source has no effect on the fluorescence data 
% collected, since it just adds a fixed baseline shift in the fluorescence 
% collected from the 570-excited (red) sensor. However, the blinking 
% results in a periodic downward artifact (in the red channel).
% 3) motion artifact (probably z-motion), seems to be present on a
% mouse-by-mouse basis. Sometimes there's a motion artifact that causes a 
% large-ish upward or downward artifact. This artifact is usually present 
% across multiple ROIs and in both channels (green and red, i.e., ACh and 
% DA), and can be bidirectional, i.e., some ROIs have a positive-going 
% artifact, while others have a negative-going artifact. This seems more
% present in the AD cohort
%
% The function below was developed iteratively by Mai-Anh Vu with 
% Anthropic's Claude Sonnet 5, July 2026.
% updated Aug 22, 2026
% 
% 
%
%
%CLASSIFYMULTICHANNELARTIFACTS Detect and classify artifacts into three
%mechanism-based types, using both simultaneously-recorded channels.
%Types are NOT mutually exclusive -- see NO MUTUAL EXCLUSION BETWEEN
%TYPES below.
%
%   result = classifyMultiChannelArtifacts(dataA, dataB) where dataA
%   and dataB are [nFrames x nROI] matrices from the two simultaneously-
%   recorded channels (e.g. DA and ACh), frame-aligned (small
%   alternation offset, e.g. ~28ms at 36Hz/18Hz, is tolerated via
%   'CrossChannelMargin'), with matching column order (column r in dataA
%   and column r in dataB must be the same physical fiber -- confirmed
%   with the user).
%
% ARTIFACT TAXONOMY (as described by the user, validated where possible
% against real recordings -- see project history):
%   TYPE 1 - random large spike/dip, OR a stuck/clipped hardware value.
%            Rare (<=2/session for the spike/dip case), very large
%            amplitude, and -- per the user -- a SPIKE: a SHARP rise and
%            fall, not a gradual one (this is what distinguishes it from
%            type 3, which can also be large but is less sharp). A stuck
%            value is treated as the deterministic extreme case of this
%            same bucket and is NOT subject to the sharpness constraint
%            (see 'Type1MinSlope'). The amplitude/slope route is
%            evaluated PER-ROI, with NO cross-ROI support whatsoever --
%            per the user, a given ROI's type-1 label must depend only on
%            that ROI's own z-trace, never on what any other ROI (or the
%            population as a whole) is doing. The stuck-value route is
%            the one exception, since a clipped/stuck value can
%            legitimately be a shared hardware event (e.g. camera-wide
%            clipping) -- that is a different mechanism, not the
%            amplitude/slope criterion the "no cross-ROI support" request
%            was about.
%            UNVALIDATED: testing 'Type1MinSlope' against every confirmed
%            case in hand (true signal, type 2, type 3) showed
%            overlapping sharpness values (roughly 3-4.5 in z/frame) --
%            we do not yet have a confirmed type-1 example to know
%            whether a real spike sits clearly above that range. Treat
%            the default as a placeholder pending a real example.
%   TYPE 2 - blinking light artifact. Subtle amplitude, based on cross-ROI
%            coincidence (a population-level event, same underlying
%            detector as type 3) and a known/expected sign (see
%            'Type2ExpectedSign', so a candidate with the wrong direction
%            doesn't qualify). Does NOT require the OTHER channel to
%            corroborate -- that is what distinguishes it from type 3.
%   TYPE 3 - motion artifact. Per the user, three properties: (1) occurs
%            across SOME PROPORTION of ROIs (see 'Type3MinROIFraction');
%            (2) CAN be sign-heterogeneous; (3) WITHIN a given ROI, both
%            channels should show it at the same timepoints. Point (3) is
%            enforced STRICTLY here: a candidate event must be
%            independently detected as a population-level event by BOTH
%            channels' own criteria (not just "does the other channel's
%            raw trace twitch somewhere"), overlapping in time. The final
%            window and ROI set is the INTERSECTION of what both channels
%            independently found, and the resulting mask is IDENTICAL for
%            both channels. This was added after the user reported that,
%            for one real session, requiring only one channel's own
%            detection (checked against the other channel's raw z-scores)
%            let a noisier channel's false positives through: that
%            channel independently flagged 22 candidates where the user
%            judged most to be true signal, while the paired channel
%            (which the user judged accurate) flagged only 6 -- all 6 of
%            which had at least one overlapping candidate in the noisier
%            channel, while only 7 of the noisy channel's 22 did in
%            reverse. Requiring mutual detection recovers all 6 genuine
%            events while dropping 15 of the 21 likely false ones.
%
% NO MUTUAL EXCLUSION BETWEEN TYPES
%   Per the user, artifact types are allowed to overlap freely -- a given
%   [frame, ROI] cell can be type1 AND type2 AND type3 simultaneously if
%   it independently meets each type's own criteria. Nothing in this
%   pipeline "consumes" or removes a candidate from consideration by one
%   type just because it already matched another. There are exactly two
%   narrower exceptions, and neither is really about type-exclusivity:
%     - A STUCK-dominated population event is excluded from type2/type3
%       evaluation entirely, because a clipped/stuck value isn't real
%       signal to evaluate against anything (consistent with stuck
%       samples being zeroed out of z-scores elsewhere in this pipeline).
%     - The type-2 SENSITIVE pass (see below) skips a candidate only if
%       it's already FULLY covered by type1Mask -- not because type1 and
%       type2 can't coexist in principle, but because a cell already
%       established as an extreme spike or a stuck value isn't a sensible
%       candidate to ALSO relabel as a subtle blinking-light event.
%
% CLASSIFICATION LOGIC (order of evaluation; NOT a priority/exclusion
% order -- see NO MUTUAL EXCLUSION above)
%   1. TYPE 1, in two parts, run in this order:
%      1a. Per-ROI, no cross-ROI support at all (classifyIsolatedType1):
%          for each ROI independently, if that ROI's own peak |z|
%          (excluding stuck samples) exceeds Type1ZThresh AND its own max
%          single-frame |z| change exceeds Type1MinSlope, OR that ROI's
%          own stuck-value overlap exceeds StuckOverlapThresh. Runs
%          first, but does NOT remove anything from consideration by
%          type 2/type 3 below -- see NO MUTUAL EXCLUSION.
%      1b. Population-gated, stuck-value only (classifyType1): for each
%          population-level event, if its AGGREGATE stuck-value overlap
%          (across all ROIs "involved" in that event) exceeds
%          StuckOverlapThresh, mark type 1 and exclude it entirely from
%          type2/type3 evaluation (a data-quality exclusion, not a
%          type-exclusivity one -- see NO MUTUAL EXCLUSION). Otherwise,
%          the event is added, unmodified, to the "remaining" list used
%          by steps 2 and 3 below.
%      .type1Reason records which route fired: 'stuck', 'amplitude', or
%      'stuck+amplitude' (both conditions true for the same ROI) from 1a;
%      'stuck' only from 1b -- amplitude alone never appears from 1b,
%      since 1b no longer evaluates amplitude/slope at all.
%   2. TYPE 3 (joint, both channels together) - for every "remaining"
%      event from 1b (in EITHER channel; no exclusion based on step 1
%      results): find events independently detected in BOTH channels
%      with overlapping time windows. Take the intersection of their
%      onset/offset and involvedROIs. Classify type 3 if the intersected
%      ROI set is at least Type3MinROIFraction of all ROIs AND the
%      intersected duration is at least Type3MinDuration. Applied
%      identically to both channels. The "remaining" lists are NOT
%      filtered afterward -- every event stays available for step 3.
%   3. TYPE 2 (per channel, independently) - for the SAME "remaining"
%      list used in step 2 (unfiltered): type 2 if the event's dominant
%      sign matches Type2ExpectedSign. Classified regardless of
%      Type2Likely -- see DETECTION vs. AUTO-CORRECTION below. A given
%      event can be labeled type2 here even if it was ALSO labeled type3
%      in step 2 above.
%   4. UNCLASSIFIED - determined PER CELL, after steps 2 and 3, as
%      whatever in the "remaining" events is covered by NEITHER
%      type2Mask NOR type3Mask (not "whatever failed type 2", since an
%      event can fail type 2 but still be fully covered by type 3, or
%      vice versa, in which case it is not unclassified).
%   Separately, any stuck-value sample that never formed a full
%   population event is still folded into type1Mask directly.
%   Separately, a SECOND, more sensitive per-channel pass
%   (Type2ZThresh/Type2PopFracThresh, classifySensitivePass) feeds
%   ADDITIONAL type-2-only candidates not caught by the main pass -- see
%   'source' field in the output events ('main', 'sensitive', or
%   'isolated' for step 1a). This runs under the same condition as step 3
%   above: only if Type2Likely is true for that channel -- see DETECTION
%   vs. AUTO-CORRECTION below.
%
% DETECTION vs. AUTO-CORRECTION
%   Type 1 and type 3 are always evaluated and reported (in .events and
%   the per-type masks: .type1Mask, .type3Mask) regardless of
%   'Type3Likely' -- that flag only controls whether .type3Mask gets
%   folded into .autoRemovedMask (and therefore .cleaned); .type1Mask
%   always contributes.
%   Type 2 is different, per the user (reversing an earlier "always
%   evaluate" design): 'Type2Likely' gates BOTH evaluation and
%   auto-correction together. If a channel is not flagged Type2Likely,
%   this pipeline does not attempt to find type-2 events on it at all --
%   .type2Mask stays entirely empty for that channel, and .events will
%   contain no type-2 entries for it. If Type2Likely is true, type 2 is
%   evaluated AND automatically folded into .autoRemovedMask/.cleaned --
%   there is no way to detect-but-not-auto-correct type 2, unlike type 3.
%
% Name-Value params, grouped by what they affect:
%
%   -- General --
%   'ChannelNames'         {default {'A','B'}} 2-element cell, for labeling
%   'Handling'             'nan' (default) | 'interp' | 'mask'
%
%   -- Shared detection infrastructure (the main population-event pass --
%      feeds the population-gated stuck check within type 1, and type 3;
%      also reused, with different threshold values, by the type-2
%      sensitive pass -- see TYPE 2 group below. NOT used by type 1's
%      amplitude/slope criterion, which is evaluated per-ROI only --
%      see TYPE 1 group) --
%   'PopZThresh'            (default 3) per-ROI threshold for the main
%                           population-event pass
%   'PopFracThresh'         (default 0.6) fraction-of-ROIs threshold for
%                           the main pass
%   'SignHetMinCount'       (default 2)
%   'ParticipationZ'        (default 2)
%   'EventGapFrames'        (default 3)
%   'BoundaryMargin'        (default 3) also reused directly by the
%                           isolated type-1 pass, to pad each candidate
%                           spike cluster -- the only shared-infrastructure
%                           parameter type 1 touches at all
%   'BoundaryFracThresh'    (default 0.3) see detectPopulationEvents
%
%   -- Type 1 (random large spike/dip, or stuck/clipped value). Per the
%      user, evaluated purely per-ROI: a given ROI's amplitude/slope
%      label never depends on any other ROI's data. (Population-level
%      events can still overlap in time with a type-1 detection, but
%      only the stuck-value route -- a genuinely shared hardware
%      mechanism, e.g. camera-wide clipping -- uses the aggregate
%      population-event window; see classifyType1 vs.
%      classifyIsolatedType1.) --
%   'Type1ZThresh'         (default 8) per-ROI |z| above which a sample is
%                           a type-1 candidate by amplitude. Used ONLY by
%                           the per-ROI isolated pass.
%   'Type1MinSlope'        (default 8, z-units/frame) required alongside
%                           Type1ZThresh, same per-ROI pass -- see
%                           UNVALIDATED note above.
%   'StuckOverlapThresh'   (default 0.2) fraction of a window's samples
%                           that must overlap detectStuckValue output to
%                           classify it type 1 via the stuck route. Used
%                           by BOTH the per-ROI isolated pass AND the
%                           population-gated pass (the latter is the one
%                           legitimate cross-ROI mechanism left in type 1,
%                           since a stuck/clipped value can be a genuine
%                           shared hardware event).
%
%   -- Type 2 (blinking light) --
%   'Type2Likely'          (default false) whether to evaluate AND
%                           auto-correct type-2 events for that channel
%                           (see DETECTION vs. AUTO-CORRECTION above --
%                           unlike Type3Likely, this gates detection
%                           itself, not just auto-correction: if false,
%                           type 2 is not attempted for that channel at
%                           all). Accepts either a scalar (applies to
%                           both channels) or a 2-element vector [A, B]
%                           matched to 'ChannelNames'.
%   'Type2ExpectedSign'    (default -1) known direction of the type-2
%                           artifact: -1 for a dip (validated case: an
%                           IR light source dropping out), +1 for a peak,
%                           or 0 to accept either. Also accepts a
%                           2-element [A, B] vector.
%   'Type2ZThresh'          (default 1.5) per-ROI threshold for a SECOND,
%                           more sensitive population-event pass, used
%                           only to find type-2 candidates (never type1/
%                           type3), and only run at all for a channel
%                           where Type2Likely is true. See rationale in
%                           project history: several confirmed type-2
%                           events involve only ~35-40% of ROIs at modest
%                           amplitude, below what the main pass can
%                           detect without also flooding it with false
%                           positives.
%   'Type2PopFracThresh'    (default 0.55) fraction threshold for the
%                           sensitive pass. Close to observed baseline
%                           noise at this ZThresh -- expect some
%                           unconfirmed candidates; see 'source' field.
%   'Type2MinDuration'      (default 8 frames) sensitive-pass candidates
%                           shorter than this are dropped.
%   'Type2MaxDuration'      (default 50 frames) sensitive-pass candidates
%                           longer than this are dropped (likely
%                           sustained real drift, not a discrete event).
%
%   -- Type 3 (motion) --
%   'Type3Likely'          (default true) whether to auto-correct type-3
%                           events (see DETECTION vs. AUTO-CORRECTION
%                           above -- type 3 is still detected and
%                           reported either way; the joint matching step
%                           always runs). Per the user, this artifact
%                           turned out to be MOUSE-SPECIFIC -- some mice
%                           show it, others don't -- so you may want
%                           false for a mouse known not to show motion
%                           artifacts, while still being able to inspect
%                           .type3Mask/.events to confirm nothing was
%                           missed. A single session-level flag, not
%                           per-channel like Type2Likely, since type 3 is
%                           inherently a joint/cross-channel decision.
%   'CrossChannelMargin'   (default 2) extra frames of tolerance when
%                           checking whether two channels' independently-
%                           detected event windows overlap in time, to
%                           allow for the alternating-exposure offset
%   'Type3MinDuration'     (default 3 frames) the INTERSECTED window must
%                           be at least this long to be classified type 3.
%   'Type3MinROIFraction'  (default 0.15) minimum fraction of ALL ROIs
%                           that must survive the intersection (of both
%                           channels' independently-detected involvedROIs)
%                           for the event to be confirmed as type 3,
%                           rather than one or two ROIs coincidentally
%                           overlapping. UNVALIDATED, see prior note.
%
% OUTPUT (struct with fields .A and .B, one per channel, using
% 'ChannelNames' if provided)
%   each channel struct contains:
%     .events        struct array; each has onset, offset, type
%                     ('type1','type2','type3','unclassified'),
%                     type1Reason ('stuck','amplitude','stuck+amplitude',
%                     or '' for non-type1 events), source ('main',
%                     'sensitive', or 'isolated' -- see classifyType1 vs.
%                     classifyIsolatedType1 vs. classifySensitivePass),
%                     peakAmplitude, maxSlope, coreDuration, stuckOverlap,
%                     involvedROIs
%     .stuckMask, .type1Mask, .type2Mask, .type3Mask, .unclassifiedMask
%                     [nFrames x nROI] logical (type3Mask is IDENTICAL
%                     between the two channels by construction). type1Mask
%                     and type3Mask always reflect everything DETECTED,
%                     regardless of Type3Likely. type2Mask does NOT: it
%                     stays entirely empty for a channel unless
%                     Type2Likely is true for it -- see DETECTION vs.
%                     AUTO-CORRECTION above. Where populated, masks are
%                     NOT mutually exclusive: the same cell can be true in
%                     more than one mask -- see NO MUTUAL EXCLUSION
%                     BETWEEN TYPES above. This is expected, not a bug.
%     .autoRemovedMask, .cleaned  -- type1Mask always included; type2Mask/
%                     type3Mask included only if Type2Likely/Type3Likely
%                     is true for this channel

ip = inputParser;
% -- General --
ip.addParameter('ChannelNames', {'A','B'});
ip.addParameter('Handling', 'nan');
% -- Shared detection infrastructure --
ip.addParameter('PopZThresh', 3);
ip.addParameter('PopFracThresh', 0.6);
ip.addParameter('SignHetMinCount', 2);
ip.addParameter('ParticipationZ', 2);
ip.addParameter('EventGapFrames', 3);
ip.addParameter('BoundaryMargin', 3);
ip.addParameter('BoundaryFracThresh', 0.3);
% -- Type 1 --
ip.addParameter('Type1ZThresh', 8);
ip.addParameter('Type1MinSlope', 8);
ip.addParameter('StuckOverlapThresh', 0.2);
% -- Type 2 --
ip.addParameter('Type2Likely', false);
ip.addParameter('Type2ExpectedSign', -1);
ip.addParameter('Type2ZThresh', 1.5);
ip.addParameter('Type2PopFracThresh', 0.55);
ip.addParameter('Type2MinDuration', 8);
ip.addParameter('Type2MaxDuration', 50);
% -- Type 3 --
ip.addParameter('Type3Likely', true);
ip.addParameter('CrossChannelMargin', 2);
ip.addParameter('Type3MinDuration', 3);
ip.addParameter('Type3MinROIFraction', 0.15);
ip.parse(varargin{:});
prm = ip.Results;

% ---- resolve Type2Likely / Type2ExpectedSign per channel ----
if isscalar(prm.Type2Likely)
    type2LikelyA = prm.Type2Likely; type2LikelyB = prm.Type2Likely;
else
    type2LikelyA = prm.Type2Likely(1); type2LikelyB = prm.Type2Likely(2);
end
if isscalar(prm.Type2ExpectedSign)
    type2SignA = prm.Type2ExpectedSign; type2SignB = prm.Type2ExpectedSign;
else
    type2SignA = prm.Type2ExpectedSign(1); type2SignB = prm.Type2ExpectedSign(2);
end
prmA = prm; prmA.Type2Likely = type2LikelyA; prmA.Type2ExpectedSign = type2SignA;
prmB = prm; prmB.Type2Likely = type2LikelyB; prmB.Type2ExpectedSign = type2SignB;

% ---- validate that per-ROI joint matching is well-defined ----
if size(dataA,2) ~= size(dataB,2)
    error('classifyMultiChannelArtifacts:roiMismatch', ...
        ['dataA and dataB have different numbers of ROIs (%d vs %d). ' ...
         'Joint type-3 matching requires both channels to index the ' ...
         'same fibers in the same order.'], size(dataA,2), size(dataB,2));
end
nROI = size(dataA,2);

% ---- shared-mechanism detectors, independent per channel ----
stuckA = detectStuckValue(dataA);
stuckB = detectStuckValue(dataB);
popA = detectPopulationEvents(dataA, 'ZThresh', prm.PopZThresh, ...
    'PopFracThresh', prm.PopFracThresh, 'SignHetMinCount', prm.SignHetMinCount, ...
    'ParticipationZ', prm.ParticipationZ, 'EventGapFrames', prm.EventGapFrames, ...
    'BoundaryMargin', prm.BoundaryMargin, 'BoundaryFracThresh', prm.BoundaryFracThresh);
popB = detectPopulationEvents(dataB, 'ZThresh', prm.PopZThresh, ...
    'PopFracThresh', prm.PopFracThresh, 'SignHetMinCount', prm.SignHetMinCount, ...
    'ParticipationZ', prm.ParticipationZ, 'EventGapFrames', prm.EventGapFrames, ...
    'BoundaryMargin', prm.BoundaryMargin, 'BoundaryFracThresh', prm.BoundaryFracThresh);
% Sensitive pass is only needed if this channel is flagged Type2Likely --
% skip the (otherwise wasted) computation entirely when it isn't, using
% an empty result with the same structure detectPopulationEvents returns.
if type2LikelyA
    sensA = detectPopulationEvents(dataA, 'ZThresh', prm.Type2ZThresh, ...
        'PopFracThresh', prm.Type2PopFracThresh, 'SignHetMinCount', nROI+1, ...
        'ParticipationZ', prm.ParticipationZ, 'EventGapFrames', prm.EventGapFrames, ...
        'BoundaryMargin', prm.BoundaryMargin, 'BoundaryFracThresh', prm.BoundaryFracThresh);
else
    sensA.events = struct('onset', {}, 'offset', {}, 'coreOnset', {}, 'coreOffset', {}, ...
        'type', {}, 'involvedROIs', {}, 'peakFraction', {}, 'peakPos', {}, 'peakNeg', {});
    sensA.z = zeros(size(dataA));
    sensA.fracInvolved = zeros(size(dataA,1), 1);
end
if type2LikelyB
    sensB = detectPopulationEvents(dataB, 'ZThresh', prm.Type2ZThresh, ...
        'PopFracThresh', prm.Type2PopFracThresh, 'SignHetMinCount', nROI+1, ...
        'ParticipationZ', prm.ParticipationZ, 'EventGapFrames', prm.EventGapFrames, ...
        'BoundaryMargin', prm.BoundaryMargin, 'BoundaryFracThresh', prm.BoundaryFracThresh);
else
    sensB.events = struct('onset', {}, 'offset', {}, 'coreOnset', {}, 'coreOffset', {}, ...
        'type', {}, 'involvedROIs', {}, 'peakFraction', {}, 'peakPos', {}, 'peakNeg', {});
    sensB.z = zeros(size(dataB));
    sensB.fracInvolved = zeros(size(dataB,1), 1);
end

[nFrames, ~] = size(dataA);
emptyEvents = struct('onset', {}, 'offset', {}, 'type', {}, 'type1Reason', {}, ...
    'peakAmplitude', {}, 'maxSlope', {}, 'coreDuration', {}, 'stuckOverlap', {}, ...
    'involvedROIs', {}, 'source', {});

stateA = initChannelState(nFrames, nROI, emptyEvents);
stateB = initChannelState(nFrames, nROI, emptyEvents);

% ---- STEP 1: ISOLATED type 1 -- per-ROI amplitude/slope, NO cross-ROI
% support at all ----
% Runs FIRST, deliberately. Per the user: type 1 should be judged purely
% on each ROI's own z amplitude and slope, with no dependence on any
% other ROI whatsoever. This scans each ROI's own z-trace directly,
% independent of any population-level criterion. (This also catches a
% real bug found via a real example: a genuine, extreme (z~21), sharp
% spike confined to a single ROI, out of many, never registered as a
% type-1 candidate under the old population-gated-only design, because
% it never crossed detectPopulationEvents' population-level gate --
% requiring either >=60% of ROIs together or >=2 ROIs on each sign --
% regardless of amplitude.)
stateA = classifyIsolatedType1(popA, stuckA, prmA, stateA);
stateB = classifyIsolatedType1(popB, stuckB, prmB, stateB);

% ---- STEP 1b: population-gated STUCK check only (NOT amplitude/slope) ----
% This step no longer makes any type-1 decision based on amplitude or
% slope -- that would reintroduce cross-ROI support (the whole event's
% classification would be driven by whichever ROI happened to have the
% most extreme value, not the ROI actually being labeled), which is
% exactly what STEP 1 above is meant to avoid. This step exists only to:
% (a) catch population-level events that are stuck-dominated (a stuck/
% clipped value CAN legitimately be a shared hardware event, e.g. camera-
% wide clipping -- that is a different mechanism than the per-ROI
% amplitude/slope criterion, and is not what the user's "no cross-ROI
% support" request was about); and (b) build the "remaining" event list
% for type-3/type-2 evaluation, EXCLUDING any ROI with ANY overlap
% already claimed by STEP 1's isolated scan, so a sample legitimately
% labeled type 1 there can't also be reconsidered for type 3/type 2.
[stateA, remA] = classifyType1(popA, stuckA, prmA, stateA);
[stateB, remB] = classifyType1(popB, stuckB, prmB, stateB);

% ---- STEP 2: type 3 (joint, mutual detection required) ----
% Always runs -- Type3Likely no longer gates detection (per the user:
% all error types should always be evaluated). It only gates whether
% .type3Mask contributes to .autoRemovedMask, in finalizeChannel below.
% Per the user, artifact types are allowed to overlap (no mutual
% exclusion) -- remA/remB are NOT filtered/consumed after this step, so
% the SAME events remain fully available for STEP 3's independent type-2
% evaluation below.
for i = 1:numel(remA)
    eA = remA(i);
    matchIdx = find(arrayfun(@(b) overlaps(eA, remB(b), prm.CrossChannelMargin), 1:numel(remB)));
    if isempty(matchIdx)
        continue;
    end
    combOnset = min(arrayfun(@(b) remB(b).onset, matchIdx));
    combOffset = max(arrayfun(@(b) remB(b).offset, matchIdx));
    combROIs = unique([remB(matchIdx).involvedROIs]);

    intOnset = max(eA.onset, combOnset);
    intOffset = min(eA.offset, combOffset);
    intROIs = intersect(eA.involvedROIs, combROIs);

    if intOffset < intOnset || isempty(intROIs)
        continue;
    end
    intDuration = intOffset - intOnset + 1;
    if intDuration < prm.Type3MinDuration || numel(intROIs)/nROI < prm.Type3MinROIFraction
        continue;
    end

    % Pad the final intersected window by CrossChannelMargin. The two
    % channels' independently hysteresis-extended boundaries don't always
    % align exactly (validated case: the known contaminated sample fell
    % one frame outside the raw intersection, since channel B's own
    % boundary recovered one frame earlier than channel A's) -- this
    % margin compensates without needing a separate parameter.
    padOnset = max(1, intOnset - prm.CrossChannelMargin);
    padOffset = min(nFrames, intOffset + prm.CrossChannelMargin);

    stateA.type3Mask(padOnset:padOffset, intROIs) = true;
    stateB.type3Mask(padOnset:padOffset, intROIs) = true;
    ev = struct('onset', padOnset, 'offset', padOffset, 'type', 'type3', ...
        'type1Reason', '', 'peakAmplitude', NaN, 'maxSlope', NaN, ...
        'coreDuration', intDuration, 'stuckOverlap', NaN, ...
        'involvedROIs', intROIs, 'source', 'main');
    stateA.events(end+1) = ev; %#ok<AGROW>
    stateB.events(end+1) = ev; %#ok<AGROW>
end

% ---- STEP 3: type 2 / unclassified fallback (per channel, independent) ----
% Runs on the FULL remA/remB (not filtered by STEP 2 above -- see note
% there), independently checking type-2 sign criteria. A cell only ends
% up unclassified if it matches NEITHER type 2 NOR type 3.
stateA = classifyType2Fallback(remA, prmA, stateA);
stateB = classifyType2Fallback(remB, prmB, stateB);

% ---- fold stuck-value into type1 (covers standalone stuck samples too) ----
stateA.type1Mask = stateA.type1Mask | stuckA;
stateB.type1Mask = stateB.type1Mask | stuckB;

% ---- STEP 4: sensitive second pass, type-2-only, per channel ----
stateA = classifySensitivePass(sensA, stuckA, prmA, stateA);
stateB = classifySensitivePass(sensB, stuckB, prmB, stateB);

chanA = finalizeChannel(dataA, stateA, stuckA, prmA);
chanB = finalizeChannel(dataB, stateB, stuckB, prmB);

result.(prm.ChannelNames{1}) = chanA;
result.(prm.ChannelNames{2}) = chanB;

end

% ===================================================================
function state = initChannelState(nFrames, nROI, emptyEvents)
state.type1Mask = false(nFrames, nROI);
state.type2Mask = false(nFrames, nROI);
state.type3Mask = false(nFrames, nROI);
state.unclassifiedMask = false(nFrames, nROI);
state.events = emptyEvents;
end

% ===================================================================
function tf = overlaps(eA, eB, margin)
tf = ~(eB.offset + margin < eA.onset || eB.onset - margin > eA.offset);
end

% ===================================================================
function [state, remaining] = classifyType1(popResult, stuckMask, prm, state)
%CLASSIFYTYPE1 Per-channel: population-gated STUCK check only. Amplitude/
%slope is handled entirely by classifyIsolatedType1 (per-ROI, no cross-
%ROI support) -- see rationale where this is called. Builds "remaining"
%for type-3/type-2 evaluation. Per the user, artifact types are allowed
%to overlap (no mutual exclusion), so ROIs already labeled type1 by the
%isolated pass are NOT excluded here -- a genuine amplitude/slope spike
%in one ROI can still also be part of a genuine cross-ROI/cross-channel
%event and get labeled type3 and/or type2 too. The one exception is a
%STUCK-dominated event (checked below): that's excluded not to enforce
%type exclusivity, but because a clipped/stuck value is not real signal
%at all, so there's nothing meaningful to evaluate it against for type2/
%type3 (consistent with how stuck samples are zeroed out of z-scores
%elsewhere in this pipeline).
remaining = struct('onset', {}, 'offset', {}, 'involvedROIs', {}, ...
    'peakPos', {}, 'peakNeg', {});

for i = 1:numel(popResult.events)
    ev = popResult.events(i);
    onset = ev.onset; offset = ev.offset; involved = ev.involvedROIs;
    if isempty(involved)
        continue;
    end

    window = stuckMask(onset:offset, involved);
    stuckOverlap = sum(window(:)) / numel(window);
    isStuck = stuckOverlap >= prm.StuckOverlapThresh;

    if isStuck
        state.type1Mask(onset:offset, involved) = true;
        state.events(end+1) = struct('onset', onset, 'offset', offset, 'type', 'type1', ...
            'type1Reason', 'stuck', 'peakAmplitude', NaN, 'maxSlope', NaN, ...
            'coreDuration', ev.coreOffset-ev.coreOnset+1, 'stuckOverlap', stuckOverlap, ...
            'involvedROIs', involved, 'source', 'main'); %#ok<AGROW>
        continue;
    end

    remaining(end+1) = struct('onset', onset, 'offset', offset, ...
        'involvedROIs', involved, 'peakPos', ev.peakPos, 'peakNeg', ev.peakNeg); %#ok<AGROW>
end
end

% ===================================================================
function state = classifyIsolatedType1(popResult, stuckMask, prm, state)
%CLASSIFYISOLATEDTYPE1 Per-ROI, independent of any population-level
%criterion: scan each ROI's own z-trace directly for sharp/large spikes
%that STEP 1 would miss because they never crossed the population gate
%(PopFracThresh/SignHetMinCount) in the first place -- i.e. artifacts
%confined to a single ROI. See rationale where this is called.
[nFrames, nROI] = size(stuckMask);
gap = 2; % small gap tolerance for grouping nearby candidate frames

for c = 1:nROI
    zc = popResult.z(:,c);
    zc(stuckMask(:,c)) = 0; % stuck contamination handled separately
    candidate = abs(zc) > prm.Type1ZThresh;
    if ~any(candidate)
        continue;
    end
    idx = find(candidate);
    splits = find(diff(idx) > gap + 1);
    starts = [1; splits + 1];
    ends = [splits; numel(idx)];
    for k = 1:numel(starts)
        clusterIdx = idx(starts(k):ends(k));
        onset = max(1, clusterIdx(1) - prm.BoundaryMargin);
        offset = min(nFrames, clusterIdx(end) + prm.BoundaryMargin);

        if all(state.type1Mask(onset:offset, c))
            continue; % already fully covered, nothing new here
        end

        window = stuckMask(onset:offset, c);
        stuckOverlap = sum(window) / numel(window);
        isStuck = stuckOverlap >= prm.StuckOverlapThresh;

        zWindow = popResult.z(onset:offset, c);
        zWindow(stuckMask(onset:offset, c)) = 0;
        peakAmp = max(abs(zWindow));

        ctxLo = max(1, onset - 1); ctxHi = min(nFrames, offset + 1);
        zCtx = popResult.z(ctxLo:ctxHi, c);
        zCtx(stuckMask(ctxLo:ctxHi, c)) = 0;
        maxSlope = max(abs(diff(zCtx)));

        isLargeAmplitude = (peakAmp > prm.Type1ZThresh) && (maxSlope > prm.Type1MinSlope);
        if ~(isStuck || isLargeAmplitude)
            continue;
        end
        if isStuck && isLargeAmplitude
            type1Reason = 'stuck+amplitude';
        elseif isStuck
            type1Reason = 'stuck';
        else
            type1Reason = 'amplitude';
        end

        state.type1Mask(onset:offset, c) = true;
        state.events(end+1) = struct('onset', onset, 'offset', offset, 'type', 'type1', ...
            'type1Reason', type1Reason, 'peakAmplitude', peakAmp, 'maxSlope', maxSlope, ...
            'coreDuration', clusterIdx(end)-clusterIdx(1)+1, 'stuckOverlap', stuckOverlap, ...
            'involvedROIs', c, 'source', 'isolated'); %#ok<AGROW>
    end
end
end

% ===================================================================
function state = classifyType2Fallback(remaining, prm, state)
%CLASSIFYTYPE2FALLBACK Per-channel: for every remaining (non-type1,
%non-stuck) event, check the type-2 sign criterion ONLY if Type2Likely is
%true for this channel. Per the user (reversing an earlier "always
%evaluate" design): type 2 should be both evaluated AND auto-included
%based on the same Type2Likely flag -- if this channel isn't flagged
%likely, type2Mask is never populated by this function at all, and every
%cell not covered by type3 (already evaluated earlier, in the joint
%matching step) is simply unclassified. When Type2Likely IS true, a cell
%is unclassified only if it matches NEITHER type 2 NOR type 3 --
%determined per-cell, not per-event, since type3Mask's window/ROI set
%can differ from this event's own. (Type 3 is unaffected by this change
%-- it is still always evaluated regardless of Type3Likely; only type 2
%reverted to flag-gated evaluation.)
for i = 1:numel(remaining)
    ev = remaining(i);
    onset = ev.onset; offset = ev.offset; involved = ev.involvedROIs;

    if prm.Type2Likely && signMatchesType2(ev.peakPos, ev.peakNeg, prm.Type2ExpectedSign)
        state.type2Mask(onset:offset, involved) = true;
        state.events(end+1) = struct('onset', onset, 'offset', offset, 'type', 'type2', ...
            'type1Reason', '', 'peakAmplitude', NaN, 'maxSlope', NaN, ...
            'coreDuration', offset-onset+1, 'stuckOverlap', NaN, ...
            'involvedROIs', involved, 'source', 'main'); %#ok<AGROW>
    end

    classified = state.type2Mask(onset:offset, involved) | state.type3Mask(onset:offset, involved);
    if ~all(classified(:))
        newUnclassified = ~classified;
        subMask = state.unclassifiedMask(onset:offset, involved);
        subMask(newUnclassified) = true;
        state.unclassifiedMask(onset:offset, involved) = subMask;
        state.events(end+1) = struct('onset', onset, 'offset', offset, 'type', 'unclassified', ...
            'type1Reason', '', 'peakAmplitude', NaN, 'maxSlope', NaN, ...
            'coreDuration', offset-onset+1, 'stuckOverlap', NaN, ...
            'involvedROIs', involved, 'source', 'main'); %#ok<AGROW>
    end
end
end

% ===================================================================
function state = classifySensitivePass(sensResult, stuckMask, prm, state)
%CLASSIFYSENSITIVEPASS Second, more sensitive population-event pass for
%subtle type-2 candidates the main pass misses. Runs ONLY if Type2Likely
%is true for this channel (per the user, reversing an earlier "always
%evaluate" design -- type 2 is now fully flag-gated, both evaluation and
%auto-correction). Per the user, type 2 is allowed to overlap type 3 (no
%mutual exclusion), so only type1Mask is checked here to skip -- not
%because of exclusivity, but because a cell already established as type
%1 (stuck, or an extreme per-ROI amplitude/slope spike) isn't a sensible
%candidate to ALSO relabel as a subtle blinking-light event; those are
%different, largely incompatible signal signatures, not just "different
%types allowed to coexist".
for i = 1:numel(sensResult.events)
    ev = sensResult.events(i);
    onset = ev.onset; offset = ev.offset; involved = ev.involvedROIs;
    if isempty(involved) || ~prm.Type2Likely
        continue;
    end
    existingCoverage = state.type1Mask(onset:offset, involved);
    if all(existingCoverage(:))
        continue;
    end
    coreDuration = ev.coreOffset - ev.coreOnset + 1;
    if coreDuration < prm.Type2MinDuration || coreDuration > prm.Type2MaxDuration
        continue;
    end
    if ~signMatchesType2(ev.peakPos, ev.peakNeg, prm.Type2ExpectedSign)
        continue;
    end
    zWindow = sensResult.z(onset:offset, :);
    zWindow(stuckMask(onset:offset, :)) = 0;
    peakAmp = max(max(abs(zWindow)));

    % Mark type2Mask only for cells NOT already covered by type1Mask --
    % not the whole block -- consistent with the "type1 is a different,
    % incompatible mechanism" rationale above (a partially-type1-covered
    % event should still get its non-type1 portion labeled type2).
    notType1 = ~state.type1Mask(onset:offset, involved);
    subMask = state.type2Mask(onset:offset, involved);
    subMask(notType1) = true;
    state.type2Mask(onset:offset, involved) = subMask;
    state.events(end+1) = struct('onset', onset, 'offset', offset, 'type', 'type2', ...
        'type1Reason', '', 'peakAmplitude', peakAmp, 'maxSlope', NaN, ...
        'coreDuration', coreDuration, 'stuckOverlap', NaN, ...
        'involvedROIs', involved, 'source', 'sensitive'); %#ok<AGROW>
end
end

% ===================================================================
function chan = finalizeChannel(data, state, stuckMask, prm)
[nFrames, nROI] = size(data);
% type1 is always auto-corrected (never optional); type2/type3 are
% detected/reported unconditionally (see classifyType2Fallback,
% classifySensitivePass, and the joint type-3 matching step above) but
% only fold into autoRemovedMask -- and therefore .cleaned -- if their
% respective Likely flag is true for this channel.
autoRemovedMask = state.type1Mask;
if prm.Type3Likely
    autoRemovedMask = autoRemovedMask | state.type3Mask;
end
if prm.Type2Likely
    autoRemovedMask = autoRemovedMask | state.type2Mask;
end

cleaned = data;
switch lower(prm.Handling)
    case 'nan'
        cleaned(autoRemovedMask) = NaN;
    case 'interp'
        cleaned(autoRemovedMask) = NaN;
        for c = 1:nROI
            col = cleaned(:,c);
            goodIdx = find(~isnan(col));
            if numel(goodIdx) >= 2
                cleaned(:,c) = interp1(goodIdx, col(goodIdx), (1:nFrames)', 'linear', 'extrap');
            end
        end
    case 'mask'
        % leave untouched
    otherwise
        error('classifyMultiChannelArtifacts:badHandling', ...
            'Handling must be ''nan'', ''interp'', or ''mask''.');
end

chan.events = state.events;
chan.stuckMask = stuckMask; % raw stuck-value flags (subset of type1Mask)
chan.type1Mask = state.type1Mask;
chan.type2Mask = state.type2Mask;
chan.type3Mask = state.type3Mask;
chan.unclassifiedMask = state.unclassifiedMask;
chan.autoRemovedMask = autoRemovedMask;
chan.cleaned = cleaned;
end

% ===================================================================
function tf = signMatchesType2(peakPos, peakNeg, expectedSign)
%SIGNMATCHESTYPE2 Check whether an event's dominant sign matches the
%known/expected sign for type-2 (blinking-light) artifacts.
if expectedSign == 0
    tf = true;
elseif expectedSign < 0
    tf = peakNeg > peakPos;
else
    tf = peakPos > peakNeg;
end
end



function result = detectPopulationEvents(data, varargin)
%DETECTPOPULATIONEVENTS Detect frame-level shared artifacts across many
%simultaneously-imaged ROIs, and share timing across participating ROIs.
%
%   result = detectPopulationEvents(data) where data is [nFrames x nROI].
%
% RATIONALE (validated on real multi-fiber data; see project history)
%   Two independent, complementary signatures distinguish shared
%   hardware/optical artifacts from genuine (possibly correlated)
%   biological signal, when many fibers are imaged in the same frame:
%
%   1) FRACTION INVOLVED: hardware-wide events (camera clipping, a light
%      source dropping out) tend to involve a large fraction of ALL
%      simultaneously-imaged ROIs at once (observed: 74-100% in
%      confirmed artifacts) whereas even strongly-correlated real
%      activity in a spatial sub-population stayed <=45% in the one
%      session tested. There is no theoretical guarantee this gap holds
%      in general -- recalibrate PopFracThresh against your own
%      confirmed-clean references before trusting it on new data.
%
%   2) SIGN HETEROGENEITY: a frame where a meaningful number of ROIs
%      cross threshold in BOTH directions simultaneously (some positive,
%      some negative) is a distinctive signature of motion-type
%      artifacts (which can affect different fibers/regions with
%      different signs), and was NOT observed in any confirmed
%      real-signal event tested (which were uniform-sign). This check
%      catches artifacts that don't reach a high involved-fraction.
%
%   BOUNDARY SHARING: once a shared event is confirmed, per-ROI
%   hysteresis-based boundaries (e.g. from remove_artifact) can disagree
%   substantially due to per-ROI noise (observed spread: 14-32 frames'
%   duration, ~9-frame edge jitter, for the same underlying event). This
%   function instead defines ONE shared onset/offset per event from the
%   population-level signal, and applies it to every ROI that shows
%   at least weak involvement (a permissive per-ROI threshold) within
%   that window -- avoiding both under- and over-extension driven by any
%   single noisy trace. ROIs with no involvement at all during the
%   window are left unmasked, since shared events do not always affect
%   every fiber (observed: "several but not all" in one confirmed case).
%
% Name-Value params:
%   'ZThresh'          per-ROI robust z-score crossing threshold (default 3)
%   'PopFracThresh'     fraction of ROIs required for a "fraction" event
%                        (default 0.6 -- calibrate against your own data;
%                        see rationale above)
%   'SignHetMinCount'   min ROIs required on EACH sign to call a frame
%                        sign-heterogeneous (default 2)
%   'ParticipationZ'    permissive per-ROI threshold used only to decide
%                        which ROIs share in an already-confirmed event's
%                        boundary (default 2, i.e. looser than ZThresh)
%   'BoundaryFracThresh' permissive population-fraction threshold used to
%                        extend the event boundary outward from the
%                        strict core, frame by frame, for as long as the
%                        fraction stays at or above it (default 0.3).
%                        This is a hysteresis step -- PopFracThresh finds
%                        the event, BoundaryFracThresh finds its true
%                        extent -- and matters most for events with a
%                        gradual, staggered recovery across ROIs, where a
%                        fixed margin alone would cut the boundary short
%                        for the slower-recovering ROIs (validated
%                        against real data: one event's fraction declined
%                        from 0.87 to 0.09 over ~10 frames; a fixed
%                        3-frame margin left ~7 frames of real
%                        contamination unmasked for the slowest ROI).
%   'EventGapFrames'    max gap (frames) to merge nearby core-event
%                        frames into one event (default 3)
%   'BoundaryMargin'    additional small fixed buffer (frames) applied on
%                        each side AFTER the BoundaryFracThresh hysteresis
%                        extension (default 3)
%
% OUTPUT (struct)
%   .events         struct array, one per detected event, with fields:
%                     .onset, .offset (final boundary: hysteresis-
%                       extended core, plus BoundaryMargin)
%                     .coreOnset, .coreOffset (the hysteresis-extended
%                       core, WITHOUT the extra BoundaryMargin padding --
%                       this is the best estimate of the event's true
%                       span, and what duration-based classification
%                       downstream should use, since .onset/.offset
%                       include the margin buffer and would overstate it)
%                     .type ('fraction','signhet','both')
%                     .involvedROIs (indices of ROIs sharing this boundary)
%                     .peakFraction, .peakPos, .peakNeg
%   .eventMask      [nFrames x nROI] logical, TRUE where a ROI should be
%                     masked/corrected due to a shared event
%   .z              [nFrames x nROI] robust z-scores (for reuse/plotting)

ip = inputParser;
ip.addParameter('ZThresh', 3);
ip.addParameter('PopFracThresh', 0.6);
ip.addParameter('SignHetMinCount', 2);
ip.addParameter('ParticipationZ', 2);
ip.addParameter('BoundaryFracThresh', 0.3);
ip.addParameter('EventGapFrames', 3);
ip.addParameter('BoundaryMargin', 3);
ip.parse(varargin{:});
prm = ip.Results;

[nFrames, nROI] = size(data);

% ---- robust per-ROI z-scores ----
z = zeros(nFrames, nROI);
for c = 1:nROI
    col = data(:,c);
    med = median(col, 'omitnan');
    mad = 1.4826 * median(abs(col - med), 'omitnan');
    if mad == 0
        continue;
    end
    z(:,c) = (col - med) / mad;
end

% ---- per-frame population features ----
posCount = sum(z > prm.ZThresh, 2);
negCount = sum(z < -prm.ZThresh, 2);
fracInvolved = (posCount + negCount) / nROI;

isFractionEvent = fracInvolved >= prm.PopFracThresh;
isSignHetEvent  = (posCount >= prm.SignHetMinCount) & (negCount >= prm.SignHetMinCount);
isCoreEvent = isFractionEvent | isSignHetEvent;

% ---- group core-event frames into events ----
events = struct('onset', {}, 'offset', {}, 'coreOnset', {}, 'coreOffset', {}, ...
    'type', {}, 'involvedROIs', {}, 'peakFraction', {}, 'peakPos', {}, 'peakNeg', {});
eventMask = false(nFrames, nROI);

idx = find(isCoreEvent);
if ~isempty(idx)
    gaps = find(diff(idx) > prm.EventGapFrames + 1);
    starts = [1; gaps + 1];
    ends = [gaps; numel(idx)];
    for k = 1:numel(starts)
        coreIdx = idx(starts(k):ends(k));
        coreStart = coreIdx(1);
        coreEnd = coreIdx(end);

        % classify type
        isFrac = any(isFractionEvent(coreIdx));
        isSign = any(isSignHetEvent(coreIdx));
        if isFrac && isSign
            evType = 'both';
        elseif isFrac
            evType = 'fraction';
        else
            evType = 'signhet';
        end

        % Boundary: extend outward from the strict core (coreStart/coreEnd)
        % as long as the population fraction stays above a PERMISSIVE
        % threshold (BoundaryFracThresh), THEN pad by a small fixed
        % BoundaryMargin. This is a hysteresis approach -- strict
        % threshold to detect the event, permissive threshold to find its
        % true extent -- the same principle as the user's own two-
        % threshold artifact-removal method, applied here to the
        % population fraction rather than a single trace.
        %
        % Rationale: a fixed margin alone under-covers events with a
        % gradual, staggered recovery across ROIs (observed: one event's
        % fraction-involved declined from 0.87 to 0.09 over ~10 frames as
        % different ROIs recovered at slightly different times -- a fixed
        % 3-frame margin stopped well before the last ROIs actually
        % returned to baseline, leaving several frames of real
        % contamination unmasked for the slower-recovering ROIs).
        hystStart = coreStart;
        while hystStart > 1 && fracInvolved(hystStart - 1) >= prm.BoundaryFracThresh
            hystStart = hystStart - 1;
        end
        hystEnd = coreEnd;
        while hystEnd < nFrames && fracInvolved(hystEnd + 1) >= prm.BoundaryFracThresh
            hystEnd = hystEnd + 1;
        end

        onset = max(1, hystStart - prm.BoundaryMargin);
        offset = min(nFrames, hystEnd + prm.BoundaryMargin);

        % which ROIs participate: any |z| > ParticipationZ within window
        windowZ = z(onset:offset, :);
        involvedROIs = find(any(abs(windowZ) > prm.ParticipationZ, 1));

        if ~isempty(involvedROIs)
            eventMask(onset:offset, involvedROIs) = true;
        end

        events(end+1) = struct( ...
            'onset', onset, 'offset', offset, ...
            'coreOnset', hystStart, 'coreOffset', hystEnd, ...
            'type', evType, ...
            'involvedROIs', involvedROIs, ...
            'peakFraction', max(fracInvolved(coreIdx)), ...
            'peakPos', max(posCount(coreIdx)), ...
            'peakNeg', max(negCount(coreIdx))); %#ok<AGROW>
    end
end

result.events = events;
result.eventMask = eventMask;
result.z = z;
result.fracInvolved = fracInvolved;
result.posCount = posCount;
result.negCount = negCount;

end

function stuckMask = detectStuckValue(data, varargin)
%DETECTSTUCKVALUE Flag exact-repeated, far-from-baseline values.
%
%   stuckMask = detectStuckValue(data) where data is
%   [nFrames x nROI]. Returns a logical mask the same size as data.
%
%   Rationale (validated against real data, see project history):
%   camera/hardware faults (e.g. saturation at the bright end, a
%   dark-current/black-level floor, or a software error code written in
%   place of a real sample) all produce the same signature: a
%   bit-identical repeated value far outside a channel's normal range.
%   Real continuous fluorescence noise essentially never repeats a value
%   exactly, regardless of which direction (high or low) the fault pins
%   to. This check is channel-agnostic and direction-agnostic, and had
%   zero false positives across every confirmed-clean file it was tested
%   against.
%
% Name-Value params:
%   'MinRepeats'  (default 3)
%   'MADThresh'   (default 6)

ip = inputParser;
ip.addParameter('MinRepeats', 3);
ip.addParameter('MADThresh', 6);
ip.parse(varargin{:});
prm = ip.Results;

[nFrames, nROI] = size(data);
stuckMask = false(nFrames, nROI);

for c = 1:nROI
    col = data(:,c);
    med = median(col, 'omitnan');
    mad = 1.4826 * median(abs(col - med), 'omitnan');
    if mad == 0
        continue;
    end
    [uvals, ~, ic] = unique(col);
    counts = accumarray(ic, 1);
    repeatCandidates = uvals(counts >= prm.MinRepeats);
    for i = 1:numel(repeatCandidates)
        v = repeatCandidates(i);
        if abs(v - med) / mad > prm.MADThresh
            stuckMask(col == v, c) = true;
        end
    end
end

end
