function vcprogress(n)
% vcprogress
% 
% Track progress using local txt file.
%
% .............................................................................
%
% Created: December 1, 2014
%
% Copyright 2014 by Vasco Curdia

%% ----------------------------------------------------------------------------

if nargin==1
    f = fopen('progress.txt', 'w');
    fprintf(f, '%d\n', n);
    fclose(f);
else
    f = fopen('progress.txt', 'a');
    fprintf(f, '1\n');
    fclose(f);
    f = fopen('progress.txt', 'r');
    progress = fscanf(f, '%d');
    fclose(f);
    percent = (length(progress)-1)/progress(1)*100;
    if ~mod(percent,10)
        fprintf('Progress: %3.0f%%\n',percent);
    end
end

%% ----------------------------------------------------------------------------

