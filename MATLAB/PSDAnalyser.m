% HR 05/03/22
% Particle size distribution analyser
% Parse, process and save PSD data from Mastersizer file

classdef PSDAnalyser
   properties
      PropName
   end
   methods
      function obj = MyClass(arg1)
         obj.PropName = arg1;
      end 
   end 
end % End of classdef

function myUtilityFcn
   fprintf('Hi')
end