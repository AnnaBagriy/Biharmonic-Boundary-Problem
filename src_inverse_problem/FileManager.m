classdef FileManager < handle
    
    methods
       
        function WriteRightVectorToFile(obj, file, y)
         vid = fopen(file, 'wt');

         fprintf(vid, '%d\n', y);

         fclose(vid);
      end
      
      function WriteMatrixToFile(obj, file, A)
         mid = fopen(file, 'wt');

         for ii = 1:size(A, 1)
            fprintf(mid, '%g\t', A(ii, :));
            fprintf(mid, '\n');
         end

         fclose(mid);
      end
      
      function y = ReadRightVectorFromFile(obj, file)
         y = dlmread(file);
      end
        
    end
    
end