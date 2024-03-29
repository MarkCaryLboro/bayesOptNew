classdef ei < acqFcn
    % Implement the expected improvement acquisition function

    properties ( Constant = true )
        FcnName         acqFcnType   = acqFcnType( "ei" )
    end % Abstract constant properties
    
    properties ( SetAccess = protected )
        Beta  (1,1) double  { mustBeGreaterThanOrEqual( Beta, 0) }          = 0.01
    end

    methods
        function obj = ei( ModelObj, Beta )
            %--------------------------------------------------------------
            % Expected Improvement (EI) class constructor
            %
            % obj = ei( ModelObj, Beta );
            %
            % Input Arguments:
            %
            % ModelObj  --> Surrogate model for system
            % Beta        --> (double) Hyperparameter {0.01}
            %--------------------------------------------------------------
            arguments
                ModelObj (1,1)          { mustBeNonempty( ModelObj ) }
                Beta     (1,1) double   { mustBeGreaterThanOrEqual( Beta, 0),...
                                          mustBeLessThanOrEqual( Beta, 1)} = 0.01
            end
            obj.Beta = Beta;
            obj.ModelObj = ModelObj;
        end % constructor

        function Fcn = evalFcn( obj, X, Coded, Beta )
            %--------------------------------------------------------------
            % Evaluate the EI acquisition function at the location
            % specified
            %
            % Fcn = obj.evalFcn( X, Beta, Coded );
            %
            % Input Arguments:
            %
            % X     --> Points to evaluate the acquisition function
            % Beta  --> Hyperparameter value
            % Coded --> True (False) if inputs are coded (decoded)
            %--------------------------------------------------------------
            arguments
                obj    (1,1)   ei        {mustBeNonempty( obj )}
                X      (:,:)   double    {mustBeNonempty( X )}
                Coded  (1,1)   logical                                      = false
                Beta   (1,1)   double    { mustBeGreaterThanOrEqual( Beta, 0),...
                                           mustBeLessThanOrEqual( Beta, 1 )} = obj.Beta
            end
            if Coded
                X = obj.ModelObj.decode( X );
            end
            obj = obj.setBeta( Beta );
            [ Z, Mu, Sigma ] = obj.calcZscore( X );
            Zpdf = normpdf( Z );
            Zcdf = normcdf( Z );
            Explore = Sigma .* Zpdf;
            Exploit = ( Mu - obj.Fmax - obj.Beta ) .* Zcdf;
            if obj.Problem
                %----------------------------------------------------------
                % Maximisation problem
                %----------------------------------------------------------
                Fcn = Exploit + Explore;                                    % add exploration bonus
                Fcn = -Fcn;                                                 % necessary as fmincon only minimises
            else
                %----------------------------------------------------------
                % Minimisation problem
                %----------------------------------------------------------
                Fcn = Exploit - Explore;                                    % subtract exploration bonus
            end 
        end % evalFcn

        function obj = setBeta( obj, Beta )
            %--------------------------------------------------------------
            % Set the hyperparameter
            %
            % obj = obj.setBeta( Beta )
            %
            % Input Arguments:
            %
            % Beta    --> (double) hyper-parameter value ( 0 <= Beta <= 1 )
            %--------------------------------------------------------------
            arguments
                obj (1,1)             {mustBeNonempty( obj )}
                Beta  (1,1)   double  { mustBeGreaterThanOrEqual( Beta, 0),...
                    mustBeLessThanOrEqual( Beta, 1 )} = 0.01
            end
            obj.Beta = Beta;
        end % setBeta

        function [ Z, Mu, Sigma ] = calcZscore( obj, X )
            %--------------------------------------------------------------
            % Calculate the adjusted standard normal variable
            %
            % [ Z, Mu, Sigma ] = obj.calcZscore( X );
            %
            % Input Arguments:
            %
            % X         --> (double) Input matrix
            %
            % Output Arguments:
            %
            % Z     - Z-scores
            % Mu    - Surrogate model evaluations g(X).
            % Sigma - Surrogate model standard deviations
            %--------------------------------------------------------------
            arguments
                obj (1,1)   ei      { mustBeNonempty( obj ) }
                X   (:,:)   double  { mustBeNonempty( X ) }                
            end
            [ Mu, Sigma ] = obj.ModelObj.predict( X );
            try
                Z = ( Mu - obj.ModelObj.Fmax - obj.Beta ) ./ Sigma;
            catch 
                Z = zeros( size( Mu ) );
            end
        end % calcZscore
    end % ordinary method signatures
end % classdef