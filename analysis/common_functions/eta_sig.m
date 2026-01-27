% Mai-Anh Vu, updated 2025/05/09
function eta_struct = eta_sig(eta_struct,varargin)

%%% optional inputs
ip = inputParser;
ip.addParameter('sig_idx',[])
ip.addParameter('sig_min',3)
ip.addParameter('sig_alpha',1)
ip.addParameter('sig_sem',0)
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

%%% load eta_struct if necessary
if ischar(eta_struct)
    if ~endsWith(eta_struct,'.mat')
        roi = [eta_struct '.mat'];
    end
    eta_struct = load(eta_struct);
end

%%% assignment
if isempty(sig_idx)
    sig_idx = [1 size(eta_struct.activity,1)];
end

%%% significance


% recalculate meanSemSig

this_sem = eta_struct.std/sqrt(size(eta_struct.activity,3));
if sig_alpha == 1
    if ~isfield(eta_struct.bootstrap,'meansU1')
        eta_struct.bootstrap.meansL1 = quantile(eta_struct.bootstrap.means,.05,3);
        eta_struct.bootstrap.meansU1 = quantile(eta_struct.bootstrap.means,.995,3);
    end
    this_upper = eta_struct.bootstrap.meansU1;
    this_lower = eta_struct.bootstrap.meansL1;
else % if sig_alpha == 5
    this_upper = eta_struct.bootstrap.meansU;
    this_lower = eta_struct.bootstrap.meansL;
end
eta_struct.bootstrap.meanSig = zeros(size(eta_struct.mean));
eta_struct.bootstrap.meanSig((eta_struct.mean+sig_sem*this_sem) < this_lower) = -1;
eta_struct.bootstrap.meanSig((eta_struct.mean-sig_sem*this_sem) > this_upper) = 1;

% use this to calculate sig
sig_cols = zeros(size(eta_struct.activity,2),2);
sig = eta_struct.bootstrap.meanSig(sig_idx(1):sig_idx(2),:);
for i = 1:size(sig,2)
    if ~isempty(strfind(sig(:,i)',ones(1,sig_min)))
        sig_cols(i,1) = 1;
    end
    if ~isempty(strfind(sig(:,i)',-1*ones(1,sig_min)))
        sig_cols(i,2) = -1;
    end
end

eta_struct.sig_cols = sig_cols;
eta_struct.sig_idx = sig_idx;