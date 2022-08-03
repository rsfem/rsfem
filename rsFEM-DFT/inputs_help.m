% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function inputs_help()
fileID = fopen('help.txt','r');
formatSpec = '%c';
S = fscanf(fileID,formatSpec);
C = strsplit(S, '#');
for i=1:length(C)
fprintf('%c', C{i});
end
fprintf('\n');
fclose(fileID);
