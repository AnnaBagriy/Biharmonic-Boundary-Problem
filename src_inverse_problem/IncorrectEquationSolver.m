classdef IncorrectEquationSolver < handle
     
    properties (SetAccess = public)
        
        n {mustBeNumeric}
        m {mustBeNumeric}
        lambda {mustBeNumeric}
        
    end
    
    properties (SetAccess = private)
        
        b
        A
        q
        
        fi_func_1
        fi_func_2
        psi_func_1
        psi_func_2
        
        a0
        a1
        a2
        
        s
        
        x1
        y1
        x2
        y2
        
        n1_on_1
        n2_on_1
        n1_on_2
        n2_on_2
        
    end
    
    properties (Access = private)
    
        calc
        f
        
    end
        
    methods
        
        function obj = IncorrectEquationSolver(n, m, lambda, f)
            obj.n = n;
            obj.m = m;
            obj.lambda = lambda;
            
            obj.f = f;
        end
        
        function InitializeDensitiesAndConstants(obj, fi_func_1, fi_func_2, psi_func_1, psi_func_2, a0, a1, a2)
            obj.fi_func_1 = fi_func_1;
            obj.fi_func_2 = fi_func_2;
            obj.psi_func_1 = psi_func_1;
            obj.psi_func_2 = psi_func_2;
            
            obj.a0 = a0;
            obj.a1 = a1;
            obj.a2 = a2;
        end
        
        function SetBoundaryAndNormal(obj, s, x1, y1, x2, y2, n11, n21, n12, n22, q_m)
            obj.s = s;
            
            obj.x1 = x1;
            obj.y1 = y1;
            obj.x2 = x2;
            obj.y2 = y2;
            
            obj.n1_on_1 = n11;
            obj.n2_on_1 = n21;
            obj.n1_on_2 = n12;
            obj.n2_on_2 = n22;
            
            obj.calc = Calculator(s);
            
            obj.calc.SetBoundary1(x1, y1, q_m);
            obj.calc.SetBoundary2(x2, y2);
            obj.calc.SetNormal1(n11, n21);
            obj.calc.SetNormal2(n12, n22);
        end
        
        function b = InitializeRightVector(obj)
            
            obj.b = zeros(2 * obj.m, 1);
            
            for ii = 1:2 * obj.m

                b_i = 0;

                for k = 1:2 * obj.m
                    b_i = b_i + obj.calc.H1(2, 1, ii, k) .* obj.fi_func_1(k) + obj.calc.H2(2, 1, ii, k) .* obj.psi_func_1(k) +...
                        (2 * obj.m * obj.calc.R(abs(ii - k), obj.n) .* obj.calc.H1_1(2, ii, k) + obj.calc.H1_2(2, ii, k)) .* obj.fi_func_2(k) +...
                        (2 * obj.m * obj.calc.R(abs(ii - k), obj.n) .* obj.calc.H2_1(2, ii, k) + obj.calc.H2_2(2, ii, k)) .* obj.psi_func_2(k);        
                end

                obj.b(ii) = obj.f(ii) - (obj.a0 + obj.a1 * obj.x2(ii) + obj.a2 * obj.y2(ii)) - b_i / (2 * obj.m);
                
            end
    
            b = obj.b; 

        end
        
        function A = InitializeMatrix(obj)
            
            obj.A = zeros(2 * obj.m, 2 * obj.n + 1);
 
            for ii = 1:2 * obj.m
                for jj = 1:2 * obj.n + 1

                    a_ij = 0;

                    for k = 1:2 * obj.m
                        a_ij = a_ij + obj.calc.l(jj, obj.s(k), obj.n) .* (obj.calc.N_r_1(ii, jj) .* obj.fi_func_1(k)) +...
                                                                          obj.calc.N_r_2(ii, jj) .* obj.psi_func_1(k);
                    end

                    obj.A(ii, jj) = a_ij / (8 * obj.m);

                end
            end
            
            A = obj.A;
            
        end
        
        function q = SolveIncorrectSystem(obj)
           
            obj.InitializeMatrix();
            obj.InitializeRightVector();
            
            AT_A = transpose(obj.A) * obj.A;
            obj.q = (AT_A + obj.lambda * eye(size(AT_A)))^(-1) * transpose(obj.A) * obj.b;
            
            q = obj.q;
            
        end
        
        function q_m = Getq_m(obj)
            obj.SolveIncorrectSystem();
            
            q_m = zeros(1, 2 * obj.m);

            for ii = 1:2 * obj.m

                q_m_k = 0;

                for jj = 1:2 * obj.n + 1
                    q_m_k = q_m_k + obj.q(jj) * obj.calc.l(jj, obj.s(ii), obj.n);
                end

                q_m(ii) = q_m_k;
                
            end
    
        end
        
    end
end

