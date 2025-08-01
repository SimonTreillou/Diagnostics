% Open the file for reading
fileID = fopen('biblio.bib', 'r');
if fileID == -1
    error('File cannot be opened. Check if the file exists and you have permission to read it.');
end

% Read the file line by line
fileContent = {};
line = fgetl(fileID);
while ischar(line)
    % Check if the line contains the word 'abstract'
    if isempty(strfind(line, 'abstract'))
        % If not, add it to the fileContent
        fileContent{end+1} = line; %#ok<SAGROW>
    end
    line = fgetl(fileID);
end
fclose(fileID);

% Open the file for writing (this will overwrite the existing file)
fileID = fopen('biblio.bib', 'w');
if fileID == -1
    error('File cannot be opened for writing. Check if you have permission to write to it.');
end

% Write the modified content back to the file
for i = 1:length(fileContent)
    fprintf(fileID, '%s\n', fileContent{i});
end
fclose(fileID);

disp('Lines containing "abstract" have been removed.');