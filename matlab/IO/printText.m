%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prints out in the MATLAB console text and variables
% INPUT PARAMETERS:
% - Boolean value to activate/deactivate the print
% - Text to print
% - Array containing printable variables for the given text

function printText(isOn, text, variables)

    if(isOn)
        fprintf(text, variables);
    end

end
