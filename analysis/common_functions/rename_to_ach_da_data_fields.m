% rename the fields in the control data for ease of code re-use
% AChMut -> ACh; tdTomato or DAMut -> DA
function data = rename_to_ach_da_data_fields(data)

% extract matching fieldnames
all_fields = fieldnames(data);
roi_fields = all_fields(~contains(all_fields,{'behav','idx','task'}));
ach_str = roi_fields{contains(roi_fields,'ACh')};
da_str = roi_fields{~contains(roi_fields,'ACh')};

% ROI fields
data.ACh = data.(ach_str);
data = rmfield(data,ach_str);
data.DA = data.(da_str);
data = rmfield(data,da_str);

% behav fields
data.behav_ACh = data.(['behav_' ach_str]);
data = rmfield(data,['behav_' ach_str]);
data.behav_DA = data.(['behav_' da_str]);
data = rmfield(data,['behav_' da_str]);