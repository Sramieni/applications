function data = load(fullPath)
% Load the object
tmp = load(fullPath,'-mat');
data = tmp.obj;
disp('Timing: Timing object read from file!');
end

