
classdef DirectProblemSolver < handle
    
    properties (Constant)
        vector_file = 'y_vector.txt';
        matrix_file = 'A_matrix.txt';
        f_file = 'f_vector.txt';
    end
    
   properties (SetAccess = public)
       
      m {mustBeNumeric}
      
      % Utility coeffs
      A0
      A1
      A2
      
      % Constant for Mx operator
      nu = 0.5

      % Functions on boundaries
      g
      q
      
      % Boundaries
       % Г1
      x_vector_1
      y_vector_1
       % Г2
      x_vector_2
      y_vector_2
   end
   
   properties (SetAccess = private)
      h % Step
      s
      
      A % Matrix
      y % Right vector
      x % Solution
      
      f % Fuction on Г2 to find
      
      % Constants
      a0
      a1
      a2
      
      % Densities
      fi_func_1
      fi_func_2
      
      psi_func_1
      psi_func_2
      
      
      x1
      y1
      
      x2
      y2
      
      n1_on_1
      n2_on_1
      
      n1_on_2
      n2_on_2
        
      x_derivative_1
      y_derivative_1
      
      x_derivative_2
      y_derivative_2
      
      x_second_derivative_1
      y_second_derivative_1
      
      x_second_derivative_2
      y_second_derivative_2
      
   end
   
   properties (Access = private)
      q_m
       
      g_k
      q_k
      
      calc
      file_manager
   end
   
   methods (Access = public)
      
      % Initialization 
      function obj = DirectProblemSolver(m, g, q, A0, A1, A2, nu)
         obj.m = m;
         
         obj.h = pi / m;
         obj.s = (0:2 * m - 1) * obj.h;
         
         obj.g = g;
         obj.q = q;
         
         obj.A0 = A0;
         obj.A1 = A1;
         obj.A2 = A2;
         
         obj.nu = nu;
         
         obj.file_manager = FileManager();
         obj.calc = Calculator(obj.s);
      end
      
      function SetBoundary1(obj, x, y)
          obj.x_vector_1 = x;
          obj.y_vector_1 = y;

          obj.x1 = x(obj.s);
          obj.y1 = y(obj.s);
          
          obj.calc.SetBoundary1(obj.x1, obj.y1);
          
          obj.SetNormal1();
      end

      function r_q_m = SetApproximationBoundary1(obj, x, y, r_q_m, q_m)
          r_q_m = r_q_m + q_m;

          obj.x1 = r_q_m .* x(obj.s);
          obj.y1 = r_q_m .* y(obj.s);
          
          obj.calc.SetBoundary1(obj.x1, obj.y1);
      end

      function SetBoundary2(obj, x, y)
          obj.x_vector_2 = x;
          obj.y_vector_2 = y;
          
          obj.x2 = x(obj.s);
          obj.y2 = y(obj.s);
          
          obj.g_k = obj.g(obj.x2, obj.y2);
          obj.q_k = obj.q(obj.x2, obj.y2);
          
          obj.calc.SetBoundary2(obj.x2, obj.y2);
          
          obj.SetNormal2();
      end

      function x = SolveSystem(obj)
         obj.InitializeMatrix();
         obj.InitializeRightVector();
         
         %obj.x = transpose(obj.y \ obj.A);
         obj.x = obj.A \ obj.y;
         x = obj.x;

      end
      
      function [fi_func_1, fi_func_2, psi_func_1, psi_func_2, a0, a1, a2] = FindDensitiesAndConstants(obj)
         obj.SolveSystem();

         obj.a0 = obj.x(8 * obj.m + 1);
         obj.a1 = obj.x(8 * obj.m + 2);
         obj.a2 = obj.x(8 * obj.m + 3);

         obj.fi_func_1 = zeros(2 * obj.m, 1);
         obj.fi_func_2 = zeros(2 * obj.m, 1);

         obj.psi_func_1 = zeros(2 * obj.m, 1);
         obj.psi_func_2 = zeros(2 * obj.m, 1);

         for ii = 1:2 * obj.m 
             obj.fi_func_1(ii) = obj.x(ii);
             obj.fi_func_2(ii) = obj.x(ii + 2 * obj.m);

             obj.psi_func_1(ii) = obj.x(ii + 4 * obj.m);
             obj.psi_func_2(ii) = obj.x(ii + 6 * obj.m);
         end
         
         a0 = obj.a0;
         a1 = obj.a1;
         a2 = obj.a2;
         
         fi_func_1 = obj.fi_func_1;
         fi_func_2 = obj.fi_func_2;
         psi_func_1 = obj.fi_func_1;
         psi_func_2 = obj.fi_func_2;
      end
      
      function f = FindF(obj)
          
          obj.f = zeros(1, 2 * obj.m);

          for ii = 1:2 * obj.m

              f_i = obj.a0 + obj.a1 * obj.x2(ii) + obj.a2 * obj.y2(ii);

              for jj = 1:2 * obj.m

                  f_i = f_i + obj.fi_func_1(jj) * obj.calc.H1(2, 1, ii, jj) / (2 * obj.m) +...
                              obj.fi_func_2(jj) * (obj.calc.H1_1(2, ii, jj) * obj.calc.R(abs(jj - ii), obj.m) + obj.calc.H1_2(2, ii, jj) / (2 * obj.m)) +...
                              obj.psi_func_1(jj) * obj.calc.H2(2, 1, ii, jj) / (2 * obj.m) +...
                              obj.psi_func_2(jj) * (obj.calc.H2_1(2, ii, jj) * obj.calc.R(abs(jj - ii), obj.m) + obj.calc.H2_2(2, ii, jj) / (2 * obj.m));

              end

              obj.f(ii) = f_i;
          end
          
          f = obj.f;
          
          obj.file_manager.WriteRightVectorToFile(obj.f_file, obj.f);
          
         % disp([newline 'f(x) (x on Г2) = ' num2str(obj.f) newline]);
          
       end
      
   end
   
   methods (Access = private)
 
      function A = InitializeMatrix(obj)

         obj.A = zeros(8 * obj.m + 3, 8 * obj.m + 3);
         
         for ii = 1:2 * obj.m   
    
            %---------------------------%

            % 1st utility equation
            obj.A(1, ii) = obj.h;
            obj.A(1, 2 * obj.m + ii) = obj.h;

            % 2nd utility equation
            obj.A(2, ii) = obj.h * obj.x1(ii);
            obj.A(2, 2 * obj.m + ii) = obj.h * obj.x2(ii);
            obj.A(2, 4 * obj.m + ii) = obj.h * obj.n1_on_1(ii);
            obj.A(2, 6 * obj.m + ii) = obj.h * obj.n1_on_2(ii);

            % 3d utility equation
            obj.A(3, ii) = obj.h * obj.y1(ii); 
            obj.A(3, 2 * obj.m + ii) = obj.h * obj.y2(ii);
            obj.A(3, 4 * obj.m + ii) = obj.h * obj.n2_on_1(ii);
            obj.A(3, 6 * obj.m + ii) = obj.h * obj.n2_on_2(ii);

            %---------------------------%

            % 1st equation linear function
            obj.A(3 + ii, 8 * obj.m + 1) = 1;
            obj.A(3 + ii, 8 * obj.m + 2) = obj.x1(ii);
            obj.A(3 + ii, 8 * obj.m + 3) = obj.y1(ii);

            % 2nd equation linear function
            obj.A(3 + 2 * obj.m + ii, 8 * obj.m + 1) = 0;
            obj.A(3 + 2 * obj.m + ii, 8 * obj.m + 2) = obj.n1_on_1(ii);
            obj.A(3 + 2 * obj.m + ii, 8 * obj.m + 3) = obj.n2_on_1(ii);

            % 3d equation linear function
            obj.A(3 + 4 * obj.m + ii, 8 * obj.m + 1) = 0;
            obj.A(3 + 4 * obj.m + ii, 8 * obj.m + 2) = obj.n1_on_2(ii);
            obj.A(3 + 4 * obj.m + ii, 8 * obj.m + 3) = obj.n2_on_2(ii);

            % 4th equation linear function
            obj.A(3 + 6 * obj.m + ii, 8 * obj.m + 1) = 0;
            obj.A(3 + 6 * obj.m + ii, 8 * obj.m + 2) = 0;
            obj.A(3 + 6 * obj.m + ii, 8 * obj.m + 3) = 0;   

            %---------------------------%

            for jj = 1:2 * obj.m
                % 1st equation kernels
                obj.A(3 + ii, jj) = obj.calc.H1_1(1, ii, jj) * obj.calc.R(abs(jj - ii), obj.m) + obj.calc.H1_2(1, ii, jj) / (2 * obj.m);
                obj.A(3 + ii, jj + 2 * obj.m) = obj.calc.H1(1, 2, ii, jj) / (2 * obj.m);
                obj.A(3 + ii, jj + 4 * obj.m) = obj.calc.H2_1(1, ii, jj) * obj.calc.R(abs(jj - ii), obj.m) + obj.calc.H2_2(1, ii, jj) / (2 * obj.m);
                obj.A(3 + ii, jj + 6 * obj.m) = obj.calc.H2(1, 2, ii, jj) / (2 * obj.m);

                % 2nd equation kernels
                obj.A(3 + ii + 2 * obj.m, jj) = obj.calc.H3_1(1, ii, jj) * obj.calc.R(abs(jj - ii), obj.m) + obj.calc.H3_2(1, ii, jj) / (2 * obj.m);
                obj.A(3 + ii + 2 * obj.m, jj + 2 * obj.m) = obj.calc.H3(1, 2, ii, jj) / (2 * obj.m);
                obj.A(3 + ii + 2 * obj.m, jj + 4 * obj.m) = obj.calc.H4_1(1, ii, jj) * obj.calc.R(abs(jj - ii), obj.m) + obj.calc.H4_2(1, ii, jj) / (2 * obj.m);
                obj.A(3 + ii + 2 * obj.m, jj + 6 * obj.m) = obj.calc.H4(1, 2, ii, jj) / (2 * obj.m);

                % 3d equation kernels
                obj.A(3 + ii + 4 * obj.m, jj) = obj.calc.H3(2, 1, ii, jj) / (2 * obj.m);
                obj.A(3 + ii + 4 * obj.m, jj + 2 * obj.m) = obj.calc.H3_1(2, ii, jj) * obj.calc.R(abs(jj - ii), obj.m) + obj.calc.H3_2(2, ii, jj) / (2 * obj.m);
                obj.A(3 + ii + 4 * obj.m, jj + 4 * obj.m) = obj.calc.H4(2, 1, ii, jj) / (2 * obj.m);
                obj.A(3 + ii + 4 * obj.m, jj + 6 * obj.m) = obj.calc.H4_1(2, ii, jj) * obj.calc.R(abs(jj - ii), obj.m) + obj.calc.H4_2(2, ii, jj) / (2 * obj.m);

                % 4th equation kernels
                obj.A(3 + ii + 6 * obj.m, jj) = obj.calc.H5(2, 1, ii, jj, obj.nu) / (2 * obj.m);
                obj.A(3 + ii + 6 * obj.m, jj + 2 * obj.m) = obj.calc.H5_1(2, ii, jj, obj.nu) * obj.calc.R(abs(jj - ii), obj.m) + obj.calc.H5_2(2, ii, jj, obj.nu) / (2 * obj.m);
                obj.A(3 + ii + 6 * obj.m, jj + 4 * obj.m) = obj.calc.H6(2, 1, ii, jj, obj.nu) / (2 * obj.m);
                obj.A(3 + ii + 6 * obj.m, jj + 6 * obj.m) = obj.calc.H6(2, 2, ii, jj, obj.nu) / (2 * obj.m);
            end

            %---------------------------%
         end
         
         % Write to file
         
         obj.file_manager.WriteMatrixToFile(obj.matrix_file, obj.A);

         A = obj.A;
      end
      
      function y = InitializeRightVector(obj)
         obj.y = zeros(8 * obj.m + 3, 1);
         
         obj.y(1) = obj.A0;
         obj.y(2) = obj.A1;
         obj.y(3) = obj.A2;
         
         for ii = 1:2 * obj.m
            obj.y(ii + 3) = 0;
            obj.y(2 * obj.m + 3 + ii) = 0;
            obj.y(4 * obj.m + 3 + ii) = obj.g_k(ii);
            obj.y(6 * obj.m + 3 + ii) =  obj.q_k(ii);
         end
         
         % Write to file
         
         obj.file_manager.WriteRightVectorToFile(obj.vector_file, obj.y);
         
         y = obj.y;
      end
       
      function SetNormal1(obj)
          syms t;

          x_derivative = matlabFunction(diff(obj.x_vector_1(t)), t);
          y_derivative = matlabFunction(diff(obj.y_vector_1(t)), t);
          
          x_2_derivative = matlabFunction(diff(diff(obj.x_vector_1(t), t), t));
          y_2_derivative = matlabFunction(diff(diff(obj.y_vector_1(t), t), t));
          
          % Ignore j
          [obj.x_derivative_1, j] = x_derivative(obj.s);
          [obj.y_derivative_1, j] = y_derivative(obj.s);

          obj.x_second_derivative_1 = x_2_derivative(obj.s);
          obj.y_second_derivative_1 = y_2_derivative(obj.s);

          obj.n1_on_1 = obj.y_derivative_1 ./ sqrt(obj.x_derivative_1.^2 + obj.y_derivative_1.^2);
          obj.n2_on_1 = - obj.x_derivative_1 ./ sqrt(obj.x_derivative_1.^2 + obj.y_derivative_1.^2);
          
          obj.calc.SetNormal1(obj.n1_on_1, obj.n2_on_1);
          obj.calc.SetDerivative1(obj.x_derivative_1, obj.y_derivative_1);
          obj.calc.SetSecondDerivative1(obj.x_second_derivative_1, obj.y_second_derivative_1);
      end
      
      function SetNormal2(obj)
          syms t;

          x_derivative = matlabFunction(diff(obj.x_vector_2(t)), t);
          y_derivative = matlabFunction(diff(obj.y_vector_2(t)), t);
          
          x_2_derivative = matlabFunction(diff(diff(obj.x_vector_2(t), t), t));
          y_2_derivative = matlabFunction(diff(diff(obj.y_vector_2(t), t), t));
          
          % Ignore j
          [obj.x_derivative_2, j] = x_derivative(obj.s);
          [obj.y_derivative_2, j] = y_derivative(obj.s);

          obj.x_second_derivative_2 = x_2_derivative(obj.s);
          obj.y_second_derivative_2 = y_2_derivative(obj.s);
          
          obj.n1_on_2 = obj.y_derivative_2 ./ sqrt(obj.x_derivative_2.^2 + obj.y_derivative_2.^2);
          obj.n2_on_2 = - obj.x_derivative_2 ./ sqrt(obj.x_derivative_2.^2 + obj.y_derivative_2.^2);
          
          obj.calc.SetNormal2(obj.n1_on_2, obj.n2_on_2);
          obj.calc.SetDerivative2(obj.x_derivative_2, obj.y_derivative_2);
          obj.calc.SetSecondDerivative2(obj.x_second_derivative_2, obj.y_second_derivative_2);
      end
   end
  
end