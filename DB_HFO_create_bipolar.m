function data_bip = DB_HFO_create_bipolar(data)

% bipolar transformation of the data
% input: data 
% output: data_bip, with the bipolar data in .x and the lables in .label
 
data_bip = data;
data_bip.label = data.label_bipolar;
data_bip = rmfield(data_bip, {'x','label_bipolar'});
for cc = 1:length(data.BipChOrder)
    
    bippy1 = data.BipChOrder(1,cc); % top line 
    bippy2 = data.BipChOrder(2,cc); % bottom line
    data_bip.x(cc,:) = data.x(bippy1,:) - data.x(bippy2,:);
 end
 
  