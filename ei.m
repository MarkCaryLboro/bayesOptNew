classdef ei < acqFcn
    % Implement the expected improvement acquisition function

    properties ( Constant = true )
        FcnName         acqFcnType   = acqFcnType( "ei" )
    end % Abstract constant properties
    
    properties ( SetAccess = protected )
        Xi  (1,1)   double { mustBePositive( Xi ), mustBeReal( Xi ) } = 0.01
    end

    methods
        function obj = ei( ModelObj, Xi )
            %--------------------------------------------------------------
            % Expected Improvement (EI) class constructor
            %
            % obj = ei( ModelObj, Xi );
            %
            % Input Arguments:
            %
            % ModelObj  --> Surrogate model for system
            % Xi        --> (double) Hyperparameter {0.01}
            %--------------------------------------------------------------
            arguments
                ModelObj (1,1)          { mustBeNonempty( ModelObj ) }
                Xi       (1,1) double   { mustBeGreaterThanOrEqual( Xi, 0),...
                                          mustBeLessThanOrEqual( Xi, 1)} = 0.01
            end
            if ( nargin > 1 )
                obj.Xi = Xi;
            end
            obj.ModelObj = ModelObj;
        end % constructor

        function Fcn = evalFcn( obj, X, Xi )
            %--------------------------------------------------------------
            % Evaluate the EI acquisition function at the location
            % specified
            %
            % Fcn = obj.evalFcn( X, Xi );
            %
            % Input Arguments:
            %
            % X     --> Points to evaluate the acquisition function
            % Xi    --> Hyperparameter value
            %--------------------------------------------------------------
            arguments
                obj (1,1)   ei          {mustBeNonempty( obj )}
                X   (:,:)   double      {mustBeNonempty( X )}
                Xi  (1,1)   double      { mustBeGreaterThanOrEqual( Xi, 0),...
                                          mustBeLessThanOrEqual( Xi, 1 )} = obj.Xi
            end
            obj = obj.setXi( Xi );
            [ Z, Mu, Sigma ] = obj.calcZscore( X );
            Zpdf = normpdf( Z );
            Zcdf = normcdf( Z );
            Fcn = -( Mu - obj.Fmax - obj.Xi ) .* Zcdf - Sigma .* Zpdf;
            Fcn = ( Sigma > 0 ) .* Fcn;
        end % evalFcn

        function obj = setXi( obj, Xi )
            %--------------------------------------------------------------
            % Set the hyperparameter
            %
            % obj = obj.setXi( Xi )
            %
            % Input Arguments:
            %
            % Xi    --> (double) hyper-parameter value ( 0 <= Xi <= 1 )
            %--------------------------------------------------------------
            arguments
                obj (1,1)   ei          {mustBeNonempty( obj )}
                Xi  (1,1)   double      { mustBeGreaterThanOrEqual( Xi, 0),...
                    mustBeLessThanOrEqual( Xi, 1 )} = 0.01
            end
            obj.Xi = Xi;
        end % setXi

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
                Z = ( Mu - obj.ModelObj.Fmax - obj.Xi ) ./ Sigma;
            catch 
                Z = zeros( size( Mu ) );
            end
        end % calcZscore
    end % ordinary method signatures
end % classdef